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
conda activate localpdb
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
passed to `make_system_caps` from the `ff19_system_maker` toolkit. Hydrogens are
stripped, the peptide termini are capped with acetyl (ACE) and N-methyl (NME)
groups, and the topology is generated with `tleap` using the ff19SB force field.
The Amber topology is converted to GROMACS format, after which the peptide is
placed in a rhombic dodecahedral box with a minimum solute-boundary distance of
1.25 nm, solvated with OPC water, and neutralised to a NaCl concentration of
0.15 M. Position restraints are generated for the protein heavy atoms, the
system is energy-minimised, and the minimised coordinates are made whole and
centred in the box. A structure that fails at any of these stages is reported
and skipped.

Each system directory therefore contains everything needed to start
equilibration: `system.gro`, `topology.top`, `posre.itp`, the minimised
structure `mini_c.gro`, and the ion and water parameter files.

> **Requires AmberTools and GROMACS.** `tleap` and the Amber force field files
> are not bundled with this repository. The pipeline was run with AmberTools 23;
> later releases are expected to work, but the ff19SB parameters must match those
> reported in the thesis. AmberTools is distributed free of charge, mostly under
> the GNU General Public License, and can be obtained from
> [ambermd.org](https://ambermd.org/GetAmber.php#ambertools). Install it into a
> dedicated conda environment and activate that environment before running this
> script, so that `tleap` and `gmx` are both on `PATH`.

The resulting directory is then replicated as many times as the replica count
assigned to that structure in the previous step, so that every copy can be
equilibrated independently with its own randomly drawn velocities. Systems are
grouped by parent structure, meaning each source sequence of a fragment keeps a
separate subtree and the replicas of one sequence are never mixed with those of
another.

The script does not overwrite existing systems: an interrupted run leaves the
partially built tree in place and must be restarted after removing
`reps_reference`. Because the path to the toolkit is resolved relative to the
working directory, the script has to be run from the directory containing
`ff19_system_maker`.

**Usage**

```bash
conda activate AmberTools23
python 1_reps_folder_maker.py fragment_1
```

**Outputs**

- `fragment_1/reps_reference/<parent>/<structure>_<k>/` — one directory per MD
  replica, containing the minimised, solvated system
- `fragment_1/all_fragment_reps.csv` — one row per replica:

  | column | description |
  | --- | --- |
  | `dir` | replica directory name, `<structure>_<k>` |
  | `method` | experimental method of the parent entry |
  | `resolution` | resolution in Å (empty for NMR) |
  | `sequence` | fragment sequence |
  | `parent` | source PDB ID and residue range |

---

### `run_selected_eq.py`

Runs the restrained equilibration of every system listed in
`all_fragment_reps_sampled.csv`.

Each replica passes through three stages of decreasing restraint strength, all
referenced to the minimised structure. Restraints are first applied to the
protein heavy atoms at 1000 kJ mol⁻¹ nm⁻² for 50 ps at constant volume, with
velocities drawn at random, so replicas of the same structure diverge from this
point onwards. Restraints are then reduced to the main-chain atoms at
500 kJ mol⁻¹ nm⁻² and finally to the backbone at 250 kJ mol⁻¹ nm⁻², both stages
at constant pressure. The final coordinates and the trajectory of the last stage
are made whole and centred in the box.

The `.mdp` files defining the three stages are read from `md_skelet/`, which is
shared by all fragments. Runs are executed sequentially on GPU 0; the loop stops
on the first failing system.

**Usage**

Run from inside the `reps_reference` directory of the fragment, with
`md_skelet` two levels above it:

```bash
cd fragment_1/reps_reference
python ../../run_selected_eq.py
```

**Outputs**

Written into each replica directory:

- `equ_2_c.gro` — final equilibrated structure, used to define the conformational
  basin and to compute the descriptors
- `equ_2.xtc` — trajectory of the last equilibration stage
- `nvt.*`, `equ_1.*`, `equ_2.*` — GROMACS output of the individual stages
