# cif_to_plip

**CIF → strict PDB → PLIP protein-ligand interaction analysis**

A single-file Python pipeline that converts a ModelCIF or mmCIF co-folding output into a format-compliant PDB file and feeds it directly into [PLIP](https://github.com/pharmai/plip) for full protein-ligand interaction profiling. Designed for downstream analysis of [Boltz-2](https://github.com/jwohlwend/boltz) and similar co-folding outputs where the input CIF carries no bond table and often contains multi-character chain IDs that break the PDB format specification.

---

## Why this script exists

Modern co-folding tools (Boltz-2, AlphaFold 3, Chai-1) emit ModelCIF files with:

- **Multi-character chain IDs** (`JH1`, `J031`, …) that exceed the single-character PDB column 22 limit
- **Residue names longer than 3 characters** (`LIG1`, `SMOL`, …) that overflow PDB columns 18–20
- **No CONECT records** — bond topology is absent from the file

Naïvely converting such files with `gemmi` or `pdb_tidy` alone introduces a subtle but critical bug: `pdb_tidy` leaves a one-serial gap at each chain/TER boundary. Because OpenBabel maps CONECT records **positionally** by serial number, any gap causes records past that point to be silently discarded. The result is corrupted ligand perception — wrong SMILES, fictitious bonds in the PyMOL session, and invented PLIP interactions. `cif_to_plip` diagnoses and resolves all of these issues before PLIP ever sees the structure.

---

## Pipeline

```
CIF input
   │
   ├─ [1] gemmi parse + sanity check (single model enforced)
   ├─ [2] Remap multi-character chain IDs → single-character (A, B, C …)
   ├─ [3] Rename >3-char residue names → 3-char codes
   ├─ [4] Write strict PDB (element cols 77-78, TER records, occupancy/B padding)
   ├─ [5] pdb_tidy → renumber all serials contiguously (removes TER-boundary gap)
   ├─ [6] Generate intramolecular CONECT records from covalent-radius geometry
   ├─ [7] pdb_validate — format compliance check
   └─ [8] PLIP → XML + TXT + PNG + PyMOL .pse
              └─ [9] Flatten XML → interactions CSV + ligand summary CSV
```

---

## Outputs

| File | Description |
|---|---|
| `<stem>.pdb` | Standards-compliant PDB with contiguous serials and CONECT records |
| `<stem>_plip/<stem>_report.xml` | Full PLIP XML report |
| `<stem>_plip/<stem>_report.txt` | Human-readable PLIP text report |
| `<stem>_plip/*.png` | Ray-traced binding-site images (one per ligand) |
| `<stem>_plip/*.pse` | PyMOL session with PLIP interaction visualisation |
| `<stem>_interactions.csv` | Tidy, flattened interactions table (one row per contact) |
| `<stem>_ligand_summary.csv` | Per-ligand summary with SMILES, InChIKey, interaction counts |

---

## Installation

Conda is recommended to manage the mixed Python/OpenBabel/PyMOL stack:

```bash
conda create -n cif2plip python=3.11 -y
conda activate cif2plip

# Core dependencies
pip install gemmi pdb-tools plip

# OpenBabel (required by PLIP for ligand perception)
conda install -c conda-forge openbabel -y

# Open-source PyMOL (required by PLIP for .pse session generation)
conda install -c conda-forge pymol-open-source -y
```

Verify the install:

```bash
python -m plip.plipcmd --help
pdb_tidy --help
```

---

## Usage

```
python cif_to_plip.py <input.cif> [-o <output_dir>] [--no-plip]
```

### Arguments

| Argument | Description |
|---|---|
| `cif` | Input `.cif` / `.mmcif` file (required) |
| `-o`, `--outdir` | Output directory (default: same directory as input) |
| `--no-plip` | Convert to PDB only; skip PLIP analysis |

### Examples

Run the full pipeline:

```bash
python cif_to_plip.py boltz_output/jak1_lig.cif -o results/jak1/
```

Convert only (useful for format validation or feeding into other tools):

```bash
python cif_to_plip.py model.cif --no-plip
```

---

## Interaction types reported

PLIP detects and reports the following contact classes, all of which appear as distinct `interaction_type` values in the output CSV:

- Hydrophobic contacts
- Hydrogen bonds (donor/acceptor, with geometry)
- Water bridges
- Salt bridges
- π–stacking (face-to-face and edge-to-face)
- π–cation interactions
- Halogen bonds
- Metal coordination

---

## Technical notes

### Chain ID remapping
Chain IDs are assigned in order of appearance using the 62-character alphabet `A–Z`, `a–z`, `0–9`. The mapping is printed to stdout at runtime.

### Residue name truncation
Names exceeding 3 characters are truncated to their first 3 characters. Collisions are resolved by assigning sequential codes `L01`, `L02`, … The mapping is printed to stdout.

### CONECT generation
Bonds are perceived from covalent radii (Cordero et al., 2008) with a 0.45 Å tolerance. H–H pairs are excluded. The bond table is written using strictly 80-column CONECT records and inserted immediately before the `END` record, after all serials have been made contiguous.

### Serial renumbering
All `ATOM`, `HETATM`, and `TER` records are renumbered sequentially from 1 with no gaps. This step is the fix for the OpenBabel CONECT-parsing bug described above.

---

## Dependencies

| Package | Role |
|---|---|
| [gemmi](https://gemmi.readthedocs.io) | CIF parsing and PDB writing |
| [pdb-tools](https://www.bonvinlab.org/pdb-tools/) | `pdb_tidy`, `pdb_validate` |
| [PLIP](https://github.com/pharmai/plip) | Protein-ligand interaction profiling |
| [OpenBabel](https://openbabel.org) | Ligand bond-order perception (PLIP dependency) |
| [PyMOL (OSS)](https://pymol.org) | `.pse` session generation (PLIP dependency) |

---

## Licence

MIT

---

## Citation

If you use this pipeline in published work, please cite PLIP:

> Adasme MF et al. (2021). PLIP 2021: expanding the scope of the protein-ligand interaction profiler to DNA and RNA. *Nucleic Acids Research*, 49(W1), W530–W534. https://doi.org/10.1093/nar/gkab294

And gemmi if the CIF conversion step is central to your workflow:

> Wojdyr M (2022). GEMMI: A library for structural biology. *Journal of Open Source Software*, 7(73), 4200. https://doi.org/10.21105/joss.04200
