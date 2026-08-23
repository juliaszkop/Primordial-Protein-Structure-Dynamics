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
