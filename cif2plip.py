#!/usr/bin/env python3
"""
cif_to_plip.py
==============
Convert a (Model)CIF structure to a strictly standards-compliant PDB file and
run a full PLIP protein-ligand interaction analysis on it.

Pipeline
--------
1. Parse the CIF with gemmi (the wwPDB's own mmCIF/PDB library).
2. Remap multi-character chain IDs to single-character PDB chain IDs
   (e.g. 'JH1' -> 'A', 'J031' -> 'B') and record the mapping.
3. Rename any hetero-compound whose residue name exceeds the 3-character PDB
   limit (e.g. 'LIG1' -> 'LIG') so the output obeys columns 18-20.
4. Write a strict PDB file (element symbols in cols 77-78, TER records, etc.).
5. Tidy with pdb_tidy, then RENUMBER all serials contiguously. This is critical:
   pdb_tidy leaves a one-number gap at the chain/TER boundary, and OpenBabel
   maps CONECT serials positionally, so any gap makes it silently discard the
   CONECT records past it. The result is a hybrid distance/CONECT perception
   that garbles the ligand (wrong SMILES, fictitious bonds in the PyMOL session,
   even invented interactions). Contiguous serials remove the gap.
6. Generate intramolecular CONECT records for every ligand from covalent-radius
   geometry (the input CIF carries no bond table), against the final serials, so
   PLIP/OpenBabel perceive ligand bond orders, rings and donors/acceptors right.
7. Validate the PDB with pdb_validate.
8. Run PLIP, emitting a PyMOL .pse session (-y), PNG images (-p), XML (-x) and
   TXT (-t) reports for every detected binding site.
9. Flatten the PLIP XML into a tidy interactions CSV + a ligand summary CSV.

The output PDB and report basenames preserve the input filename stem.

Dependencies: gemmi, pdb-tools, plip, openbabel, pymol (open-source).
Author: pipeline scaffolded for Marc Deller.
"""
from __future__ import annotations
import argparse, csv, os, string, subprocess, sys
from collections import defaultdict
import gemmi

# --- Covalent radii (Angstrom), Cordero et al. 2008; enough for organics + common heteroatoms.
COVALENT_RADII = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "P": 1.07, "S": 1.05, "CL": 1.02, "BR": 1.20, "I": 1.39, "SI": 1.11,
    "SE": 1.20, "AS": 1.19, "NA": 1.66, "K": 2.03, "MG": 1.41,
    "CA": 1.76, "FE": 1.32, "ZN": 1.22, "MN": 1.39, "CU": 1.32, "CO": 1.26,
    "NI": 1.24,
}
BOND_TOLERANCE = 0.45          # Angstrom slack added to summed covalent radii
DEFAULT_RADIUS = 0.77          # fallback if an element is missing from the table

PDB_CHAIN_ALPHABET = (string.ascii_uppercase + string.ascii_lowercase
                      + string.digits)


def remap_structure(cif_path: str):
    """Read the CIF and return (structure, chain_map, resname_map, ligand_resnames)."""
    st = gemmi.read_structure(cif_path)
    st.setup_entities()                      # ensure subchains/het flags are sane
    if len(st) == 0:
        sys.exit("ERROR: no models found in the input file.")
    if len(st) > 1:
        # PLIP analyses a single model; keep model 1 only and warn.
        print(f"[warn] {len(st)} models present; keeping model 1 only.")
        while len(st) > 1:
            del st[1]

    model = st[0]

    # --- 1. Single-character chain IDs, assigned in order of appearance.
    chain_map: dict[str, str] = {}
    used = set()
    alphabet = iter(PDB_CHAIN_ALPHABET)
    for ch in model:
        if ch.name in chain_map:
            continue
        nxt = next(a for a in alphabet if a not in used)
        chain_map[ch.name] = nxt
        used.add(nxt)
    for ch in model:
        ch.name = chain_map[ch.name]

    # --- 2. Shorten over-length hetero residue names (PDB cols 18-20 = 3 chars).
    resname_map: dict[str, str] = {}
    ligand_resnames: set[str] = set()
    used_codes = {r.name for ch in model for r in ch if len(r.name) <= 3}
    for ch in model:
        for res in ch:
            is_het = res.het_flag == "H"
            if is_het:
                ligand_resnames.add(res.name)
            if len(res.name) > 3 and res.name not in resname_map:
                base = res.name[:3].upper()
                code, n = base, 0
                while code in used_codes and code != resname_map.get(res.name):
                    n += 1
                    code = f"L{n:02d}"
                resname_map[res.name] = code
                used_codes.add(code)
    if resname_map:
        for ch in model:
            for res in ch:
                if res.name in resname_map:
                    res.name = resname_map[res.name]
        # update ligand set to the post-rename codes
        ligand_resnames = {resname_map.get(x, x) for x in ligand_resnames}

    st.setup_entities()
    return st, chain_map, resname_map, ligand_resnames


def write_pdb(st: gemmi.Structure, out_pdb: str):
    """Write a strict PDB; gemmi handles element cols, TER, occupancy/B padding."""
    doc_opts = gemmi.PdbWriteOptions()
    st.write_pdb(out_pdb, doc_opts)


def _element_of(atom_name: str, raw_elem: str) -> str:
    e = (raw_elem or "").strip().upper()
    if e:
        return e
    # crude fallback from atom name
    a = atom_name.strip().lstrip("0123456789")
    return (a[:2] if a[:2].upper() in COVALENT_RADII else a[:1]).upper()


def renumber_contiguous(pdb_path: str):
    """Renumber every ATOM/HETATM/TER serial sequentially with NO gaps.

    pdb_tidy assigns the inter-chain TER its own serial but then restarts the
    next chain one number too high, leaving a gap (e.g. ...2239, TER 2240,
    HETATM 2242). OpenBabel maps CONECT serials positionally, so a gap shifts
    the map and CONECT records pointing past the gap are silently DISCARDED --
    which is what corrupts ligand bond perception (garbled SMILES, the PyMOL
    distance-bond tangle, even fictitious interactions). Contiguous serials
    remove the gap entirely. Any existing CONECT records are dropped here; they
    are regenerated against the new serials by add_ligand_conect().
    """
    out = []
    new = 0
    for ln in open(pdb_path):
        rec = ln[:6]
        if rec in ("ATOM  ", "HETATM", "TER   "):
            new += 1
            ln = ln[:6] + f"{new:>5}" + ln[11:]
            out.append(ln)
        elif ln.startswith("CONECT"):
            continue                      # stale; regenerated downstream
        else:
            out.append(ln)
    with open(pdb_path, "w") as fh:
        fh.writelines(out)


def add_ligand_conect(pdb_path: str, ligand_resnames: set[str]):
    """Append CONECT records for intramolecular ligand bonds (covalent-radius rule)."""
    if not ligand_resnames:
        return 0
    atoms = []  # (serial, element, x, y, z)
    with open(pdb_path) as fh:
        lines = fh.readlines()
    for ln in lines:
        if ln.startswith("HETATM"):
            resname = ln[17:20].strip()
            if resname in ligand_resnames:
                serial = int(ln[6:11])
                x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
                elem = _element_of(ln[12:16], ln[76:78])
                atoms.append((serial, elem, x, y, z))
    bonds = defaultdict(set)
    for i in range(len(atoms)):
        si, ei, xi, yi, zi = atoms[i]
        ri = COVALENT_RADII.get(ei, DEFAULT_RADIUS)
        for j in range(i + 1, len(atoms)):
            sj, ej, xj, yj, zj = atoms[j]
            rj = COVALENT_RADII.get(ej, DEFAULT_RADIUS)
            d2 = (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2
            cutoff = ri + rj + BOND_TOLERANCE
            if ei == "H" and ej == "H":
                continue
            if 0.16 < d2 <= cutoff * cutoff:   # >0.4 A apart, within covalent range
                bonds[si].add(sj)
                bonds[sj].add(si)
    conect = []
    for serial in sorted(bonds):
        partners = sorted(bonds[serial])
        for k in range(0, len(partners), 4):       # max 4 partners per CONECT line
            chunk = partners[k:k + 4]
            rec = "CONECT" + f"{serial:>5}" + "".join(f"{p:>5}" for p in chunk)
            conect.append(rec.ljust(80) + "\n")   # strict 80-column width
    # insert CONECT block immediately before END / MASTER
    out, inserted = [], False
    for ln in lines:
        if ln.startswith(("END", "MASTER")) and not inserted:
            out.extend(conect)
            inserted = True
        out.append(ln)
    if not inserted:
        out.extend(conect)
        out.append("END\n")
    with open(pdb_path, "w") as fh:
        fh.writelines(out)
    n_bonds = sum(len(v) for v in bonds.values()) // 2
    return n_bonds


def run_cli(cmd: list[str], **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def plip_xml_to_csv(xml_path: str, inter_csv: str, summary_csv: str):
    """Flatten a PLIP XML report into a tidy interactions CSV + a ligand summary CSV."""
    import xml.etree.ElementTree as ET
    root = ET.parse(xml_path).getroot()

    # headline-distance tag, in priority order, per interaction class
    DIST_KEYS = ["dist", "dist_d-a", "dist_a-w", "centdist", "dist_h-a"]
    # only unambiguous protein-side / ligand-side atom indices go in the idx columns;
    # donor/acceptor indices (H-bonds, water bridges) are ambiguous so they stay in geometry
    LIG_IDX = ["ligcarbonidx", "lig_idx_list", "lig_idx", "metal_idx"]
    PROT_IDX = ["protcarbonidx", "prot_idx_list", "target_idx"]
    SKIP = {"resnr", "restype", "reschain", "resnr_lig", "restype_lig",
            "reschain_lig", "ligcoo", "protcoo"} | set(LIG_IDX) | set(PROT_IDX)

    rows = []
    summary = []
    for site in root.findall(".//bindingsite"):
        ident = site.find("identifiers")
        gett = lambda tag: (ident.findtext(tag) or "").strip() if ident is not None else ""
        lig_id = f"{gett('hetid')}:{gett('chain')}:{gett('position')}"
        counts = defaultdict(int)
        interactions = site.find("interactions")
        for category in (list(interactions) if interactions is not None else []):
            itype = category.tag.replace("_interactions", "").replace("_", " ").strip()
            for entry in category:
                fields = {c.tag: (c.text or "").strip() for c in entry}
                dist = next((fields[k] for k in DIST_KEYS if fields.get(k)), "")
                lig_idx = next((fields[k] for k in LIG_IDX if fields.get(k)), "")
                prot_idx = next((fields[k] for k in PROT_IDX if fields.get(k)), "")
                geom = "; ".join(f"{k}={v}" for k, v in fields.items()
                                 if k not in SKIP and v not in ("", None))
                rows.append({
                    "ligand": lig_id,
                    "interaction_type": itype,
                    "prot_resnr": fields.get("resnr", ""),
                    "prot_restype": fields.get("restype", ""),
                    "prot_chain": fields.get("reschain", ""),
                    "lig_restype": fields.get("restype_lig", ""),
                    "lig_chain": fields.get("reschain_lig", ""),
                    "distance_A": dist,
                    "prot_atom_idx": prot_idx,
                    "lig_atom_idx": lig_idx,
                    "geometry": geom,
                })
                counts[itype] += 1
        summary.append({
            "ligand": lig_id,
            "ligand_type": gett("ligtype"),
            "smiles": gett("smiles"),
            "inchikey": gett("inchikey"),
            "total_interactions": sum(counts.values()),
            **{f"n_{k.replace(' ', '_')}": v for k, v in counts.items()},
        })

    cols = ["ligand", "interaction_type", "prot_resnr", "prot_restype",
            "prot_chain", "lig_restype", "lig_chain", "distance_A",
            "prot_atom_idx", "lig_atom_idx", "geometry"]
    with open(inter_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    scols = sorted({k for s in summary for k in s},
                   key=lambda c: (c not in ("ligand", "ligand_type",
                                            "total_interactions"), c))
    with open(summary_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=scols)
        w.writeheader()
        for s in summary:
            w.writerow({k: s.get(k, 0) for k in scols})
    return len(rows), len(summary)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="CIF -> strict PDB -> PLIP analysis")
    ap.add_argument("cif", help="input .cif / .mmcif file")
    ap.add_argument("-o", "--outdir", default=None,
                    help="output directory (default: alongside input)")
    ap.add_argument("--no-plip", action="store_true", help="convert only")
    args = ap.parse_args()

    cif_path = os.path.abspath(args.cif)
    stem = os.path.splitext(os.path.basename(cif_path))[0]
    outdir = os.path.abspath(args.outdir or os.path.dirname(cif_path))
    os.makedirs(outdir, exist_ok=True)
    pdb_path = os.path.join(outdir, f"{stem}.pdb")

    print(f"[1/6] parsing + remapping  {cif_path}")
    st, chain_map, resname_map, ligs = remap_structure(cif_path)
    print(f"      chain map     : {chain_map}")
    print(f"      resname map   : {resname_map or '(none)'}")
    print(f"      ligand codes  : {sorted(ligs) or '(none)'}")

    print(f"[2/6] writing PDB          {pdb_path}")
    write_pdb(st, pdb_path)

    print(f"[3/6] tidying PDB (TER records, atom renumbering, column padding)")
    tidy = run_cli(["pdb_tidy", pdb_path])
    if tidy.returncode == 0 and tidy.stdout:
        with open(pdb_path, "w") as fh:
            fh.write(tidy.stdout)
    # Remove the serial gap pdb_tidy leaves at chain/TER boundaries, otherwise
    # OpenBabel discards CONECT records past the gap and mis-perceives the ligand.
    renumber_contiguous(pdb_path)

    print(f"[4/6] generating ligand CONECT records (contiguous post-tidy serials)")
    nb = add_ligand_conect(pdb_path, ligs)
    print(f"      {nb} intramolecular ligand bonds written")

    print(f"[5/6] validating PDB against format rules")
    val = run_cli(["pdb_validate", pdb_path])
    warn = (val.stdout + val.stderr).strip()
    print("      " + (warn.replace("\n", "\n      ") if warn else "no validation warnings"))

    if args.no_plip:
        print("done (conversion only).")
        sys.exit(0)

    print(f"[6/6] running PLIP")
    plip_dir = os.path.join(outdir, f"{stem}_plip")
    os.makedirs(plip_dir, exist_ok=True)
    # PLIP 3.x flags:  -p pictures(PNG)  -y PyMOL .pse session  -x XML  -t TXT
    plip_cmd = [sys.executable, "-m", "plip.plipcmd",
                "-f", pdb_path, "-o", plip_dir, "-pyxt"]
    res = run_cli(plip_cmd, env={**os.environ})
    xml_path = os.path.join(plip_dir, f"{stem}_report.xml")
    if not os.path.isfile(xml_path):
        print("      PLIP did not produce an XML report. stderr tail:")
        print("      " + "\n      ".join(res.stderr.splitlines()[-8:]))
        sys.exit(1)

    inter_csv = os.path.join(outdir, f"{stem}_interactions.csv")
    summ_csv = os.path.join(outdir, f"{stem}_ligand_summary.csv")
    n_rows, n_sites = plip_xml_to_csv(xml_path, inter_csv, summ_csv)
    print(f"      {n_sites} binding site(s), {n_rows} interactions")

    print("\nDeliverables")
    print("-" * 60)
    print(f"  PDB (strict + CONECT)   : {pdb_path}")
    for f in sorted(os.listdir(plip_dir)):
        tag = {".pse": "PyMOL session", ".png": "ray-traced image",
               ".xml": "PLIP XML report", ".txt": "PLIP text report"}.get(
                   os.path.splitext(f)[1], "")
        if tag:
            print(f"  {tag:<23} : {os.path.join(plip_dir, f)}")
    print(f"  interactions CSV        : {inter_csv}")
    print(f"  ligand summary CSV      : {summ_csv}")
    print("done.")
