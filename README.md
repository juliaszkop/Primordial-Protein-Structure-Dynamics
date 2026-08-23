## Input directory structure

The pipeline operates on a per-fragment basis. Before any processing, the
following structure must be prepared under `fragment_generator/`:

```text
fragment_generator/
├── PDB_sequences/
│   └── pdb_seqres.txt                 # FASTA file with all PDB sequences (shared across fragments)
├── fragment_1/
│   └── ancestral_data_no_gaps.dat     # source sequences and PDB IDs for this fragment
├── fragment_2/
│   └── ancestral_data_no_gaps.dat
├── fragment_3/
│   └── ancestral_data_no_gaps.dat
└── ...
```

### File descriptions

- `fragment_N/ancestral_data_no_gaps.dat` — whitespace-separated table with
  columns `pdb_id`, `sequence`, `start_end`, `scope`, derived from the SI of
  Alva et al. (2015); one file per fragment
- `PDB_sequences/pdb_seqres.txt` — PDB sequence database in FASTA format,
  available from the RCSB PDB website

### Additional requirements

- A local PDB mirror managed by the `localpdb` package (expected at
  `/data/LOCAL_PDB`), set up with the `PDBChain` plugin

  ### `0_fragment_reps_generator.py`

Extracts experimental structures of a given ancestral fragment from the PDB and
allocates MD replicas across them.

For each source sequence listed in `ancestral_data_no_gaps.dat`, the script
searches `pdb_seqres.txt` for exact substring matches, resolves the matching
chains against the local PDB mirror, and locates the fragment within each
structure by scanning the sequence of CA-containing residues. Every hit is
written as a separate PDB file restricted to a single model and chain, with
alternate locations reduced to conformer A.

Because NMR entries contribute multiple models of the same chain, the MD budget
is distributed across hits rather than assigned per file: each source sequence
receives a fixed total of 25 replicas, allocated so that every hit is sampled at
least once and the remainder is drawn with replacement. When a sequence yields
more than 25 hits, 25 are drawn without replacement.

All extracted fragments are then superimposed on main-chain atoms (N, CA, C, O)
using the Kabsch algorithm. The first structure serves as the initial reference,
an average structure is computed from the superimposed coordinates, and every
fragment is re-aligned onto that average.

**Usage**

```bash
python 0_fragment_reps_generator.py fragment_1 
```

**Outputs**

- `fragment_1_done/` — one PDB file per extracted hit, aligned to the
  average structure, plus `average.pdb`
- `fragment_1/fragment_templates.csv` — one row per extracted structure:

  | column | description |
  | --- | --- |
  | `pdb_file` | file name of the extracted fragment |
  | `counts` | number of MD replicas allocated to this structure |
  | `method` | experimental method of the parent entry |
  | `resolution` | resolution in Å (empty for NMR) |
  | `sequence` | fragment sequence used as the search query |
  | `scope` | SCOP classification of the source domain |
  | `parent` | source PDB ID and residue range |

  ### `1_reps_folder_maker.py`

Builds the simulation systems for every MD replica of a given ancestral
fragment.

For each structure listed in `fragment_templates.csv`, the extracted PDB file is
passed to `make_system_caps` from the `ff19_system_maker` toolkit. The peptide
termini are capped with acetyl and N-methyl groups, the topology is generated
with the ff19SB force field, and the peptide is placed in a rhombic dodecahedral
box, solvated with OPC water, neutralised and brought to 0.15 M NaCl. A
structure that fails to parameterise is reported and skipped.

The resulting directory is then replicated as many times as the replica count
assigned to that structure in the previous step, so that every copy can be
minimised and equilibrated independently with its own randomly drawn velocities.
Systems are grouped by parent structure, meaning each source sequence of a
fragment keeps a separate subtree and the replicas of one sequence are never
mixed with those of another.

**Usage**

```bash
python 1_reps_folder_maker.py fragment_1
```

**Outputs**

- `fragment_1/reps_reference/<parent>/<structure>_<k>/` — one directory per MD
  replica, containing the solvated system ready for energy minimisation
- `fragment_1/all_fragment_reps.csv` — one row per replica:

  | column | description |
  | --- | --- |
  | `dir` | replica directory name, `<structure>_<k>` |
  | `method` | experimental method of the parent entry |
  | `resolution` | resolution in Å (empty for NMR) |
  | `sequence` | fragment sequence |
  | `parent` | source PDB ID and residue range |
