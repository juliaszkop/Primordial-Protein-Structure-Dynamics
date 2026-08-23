#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import SeqIO
import pandas as pd
from localpdb import PDB

import MDAnalysis as mda

import numpy as np
import os
import gzip
import argparse
from typing import List, Tuple

from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select, is_aa
from Bio.SeqUtils import seq1

import random
from collections import Counter


# In[ ]:


parser = argparse.ArgumentParser(
    description="Extract experimental structures of an ancestral fragment from the PDB "
                "and allocate MD replicas across them."
)
parser.add_argument(
    "fragment_dir",
    help="directory of the fragment to process, e.g. fragment_37",
)
parser.add_argument(
    "-o", "--out-dir",
    default=None,
    help="directory for the extracted structures "
         "(default: <fragment_dir>_done, next to the fragment directory)",
)
parser.add_argument(
    "--seqres",
    default="PDB_sequences/pdb_seqres.txt",
    help="FASTA file with PDB sequences (default: %(default)s)",
)
args = parser.parse_args()

fragment_dir = os.path.normpath(args.fragment_dir)
if args.out_dir is not None:
    fragment_templates_dir = os.path.normpath(args.out_dir)
else:
    fragment_templates_dir = fragment_dir + "_done"


#load LOCAL PDB into dataframe, including clustering30 data 
lpdb = PDB(db_path='/data/LOCAL_PDB', plugins=['PDBChain'])
# lpdb.load_clustering_data(redundancy=50)
# lpdb.load_clustering_data(redundancy=100)


# In[ ]:


_seqres_cache = {}


def load_seqres(pdb_fasta_file):
    """
    read a fasta file once and keep it in memory as a list of (id, sequence)
    """
    records = _seqres_cache.get(pdb_fasta_file)
    if records is None:
        records = [
            (record.id, str(record.seq))
            for record in SeqIO.parse(pdb_fasta_file, "fasta")
        ]
        _seqres_cache[pdb_fasta_file] = records
    return records


def get_matches(pdb_fasta_file,query):
    """
    search fasta file for matching sequence string
    """
    for rid, seq in load_seqres(pdb_fasta_file):
        if query not in seq:
            continue
    
        # Handle common formats
        if "|" in rid:
            parts = rid.split("|")
            if len(parts) >= 3:
                pdb_id = parts[1].lower()
                chain = parts[2]
            else:
                continue
        elif "_" in rid:
            pdb_id, chain = rid.split("_", 1)
            pdb_id = pdb_id.lower()
        else:
            # fallback: first 4 chars = pdb id
            pdb_id = rid[:4].lower()
            chain = None
    
        matches.append((pdb_id+"_"+chain))
    return matches        


# In[ ]:


def md_counts_from_nmr_models(models, total_md=25, seed=None):
    """
    models: list of model identifiers (e.g. model indices or names)
    total_md: total number of MD relaxations to allocate
    seed: optional RNG seed for reproducibility

    returns: dict {model: number_of_md_runs}
    """
    if seed is not None:
        random.seed(seed)

    M = len(models)
    if M == 0:
        raise ValueError("Model list is empty")

    # Start by including each model once (if possible)
    counts = Counter(models)

    if M >= total_md:
        # If there are more models than MD budget,
        # randomly select total_md models without replacement
        selected = random.sample(models, total_md)
        return Counter(selected)

    # Otherwise, sample remaining with replacement
    remaining = total_md - M
    extra_samples = random.choices(models, k=remaining)
    counts.update(extra_samples)

    return dict(counts)


# In[ ]:


# Below set functions to extract fragments from pdb files based on sequence

def load_structure(path: str, struct_id: str = "X"):
    lower = path.lower()
    base = lower[:-3] if lower.endswith(".gz") else lower
    parser = MMCIFParser(QUIET=True) if base.endswith((".cif", ".mmcif")) else PDBParser(QUIET=True)
    if lower.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
            return parser.get_structure(struct_id, handle)
    return parser.get_structure(struct_id, path)


def chain_residues_with_coords(chain):
    residues, letters = [], []
    for res in chain.get_residues():
        if not is_aa(res, standard=False):
            continue
        if "CA" not in res:
            continue
        try:
            letters.append(seq1(res.get_resname()))
        except Exception:
            continue
        residues.append(res)
    return residues, "".join(letters)


def find_fragment_hits(structure, query_seq: str):
    query_seq = query_seq.strip().upper()
    hits = []
    for model in structure:
        for chain in model:
            residues, seq = chain_residues_with_coords(chain)
            start = 0
            while True:
                i = seq.find(query_seq, start)
                if i == -1:
                    break
                frag_res = residues[i:i + len(query_seq)]
                hits.append((model.id, chain.id, frag_res, i))
                start = i + 1
    return hits


class FragmentSelectByModelChain(Select):
    """
    Select fragment residues, but also restrict to a given NMR model_id and chain_id.
    """
    def __init__(self, model_id, chain_id, allowed_residues, keep_altloc="A_or_blank"):
        super().__init__()
        self.model_id = model_id
        self.chain_id = chain_id
        self.allowed = set(allowed_residues)
        self.keep_altloc = keep_altloc

    def accept_model(self, model):
        return 1 if model.id == self.model_id else 0

    def accept_chain(self, chain):
        return 1 if chain.id == self.chain_id else 0

    def accept_residue(self, residue):
        return 1 if residue in self.allowed else 0

    def accept_atom(self, atom):
        if self.keep_altloc == "all":
            return 1
        alt = atom.get_altloc()
        if self.keep_altloc == "A_or_blank":
            return 1 if alt in (" ", "A") else 0
        if self.keep_altloc == "blank_only":
            return 1 if alt == " " else 0
        raise ValueError("keep_altloc must be one of: 'all', 'A_or_blank', 'blank_only'")


def write_fragment_pdbs(
    input_path: str,
    query_seq: str,
    out_dir: str,
    out_prefix: str,
    keep_altloc: str = "A_or_blank",
):
    os.makedirs(out_dir, exist_ok=True)
    structure = load_structure(input_path, struct_id="IN")
    hits = find_fragment_hits(structure, query_seq)

    io = PDBIO()
    io.set_structure(structure)   # <-- IMPORTANT: keep full structure as root

    out_files = []
    _scores = []
    for k, (model_id, chain_id, frag_res, start_idx) in enumerate(hits, start=1):
        selector = FragmentSelectByModelChain(
            model_id=model_id,
            chain_id=chain_id,
            allowed_residues=frag_res,
            keep_altloc=keep_altloc,
        )

        out_name = f"{out_prefix}_{start_idx+1}_{start_idx+len(query_seq)}_{k}.pdb"
        out_path = os.path.join(out_dir, out_name)

        io.save(out_path, select=selector)
        # print(out_path)

        # keep only non-empty outputs (END-only files are tiny)
        if os.path.getsize(out_path) < 50:
            os.remove(out_path)
            continue

        out_files.append(out_name)

    if len(out_files) > 0:
        out_data =  md_counts_from_nmr_models(out_files,total_md=25)
        return out_data
    else:
        return 0


# In[ ]:


def convert_start_end(s):
    return (
        s.strip("()")
         .replace(":", "_")
         .replace("-", "_")
    )


# In[ ]:


ref_file = "ancestral_data_no_gaps.dat" #sequences and pdb-ids from SI of e-life paper for the given fragment
ref_fragment_file = os.path.join(fragment_dir,ref_file)
ref_fragments = pd.read_csv(ref_fragment_file,sep= r"\s+",index_col="pdb_id")
ref_fragments


# In[ ]:


# search for matching sequences & extact fragments
os.makedirs(fragment_templates_dir, exist_ok=True)
fragment_reps = pd.DataFrame()
for pdb_id,ref_fragment in ref_fragments.iterrows():
    sequence = ref_fragment["sequence"]
    res_range = convert_start_end(ref_fragment["start_end"])
    matches = []
    matches = get_matches(args.seqres,sequence)
    # print(matches)
    structures = lpdb.chains.loc[(lpdb.chains.index.isin(matches)) & (lpdb.chains["pdb_fn"].notna())]
    # print("SEQ: ",sequence)
    for chain,row in structures.iterrows():
        # print(chain,sequence)
        pdb_names = write_fragment_pdbs(
            input_path=row["pdb_fn"],
            query_seq=sequence,
            out_dir=fragment_templates_dir,
            keep_altloc="A_or_blank", 
            out_prefix = chain
        )
        if pdb_names == 0:
            continue
        fragment_reps_tmp = pd.DataFrame(pdb_names.items(),columns=["pdb_file","counts"])
        fragment_reps_tmp["method"] = row["method"]
        fragment_reps_tmp["resolution"] = row["resolution"]
        fragment_reps_tmp["sequence"] = sequence
        fragment_reps_tmp["scope"] = ref_fragment["scope"]
        fragment_reps_tmp["parent"] = pdb_id+"_"+res_range
        fragment_reps = pd.concat([fragment_reps,fragment_reps_tmp])    

fragment_reps_file_name = "fragment_templates.csv"
fragment_reps_file = os.path.join(fragment_dir,fragment_reps_file_name)
fragment_reps.to_csv(fragment_reps_file, index=False)


# In[ ]:


# two functions for RMSD-based structure alignment
def centroid(X):
    return X.mean(axis=0)
def kabsch_rotate(A, B):
    """
    Compute rotation matrix R that rotates A onto B (both shape (N,3)).
    Returns R (3x3), centroid_A, centroid_B.
    """
    # centroids
    cA = centroid(A)
    cB = centroid(B)
    # center coords
    A0 = A - cA
    B0 = B - cB
    # covariance
    C = A0.T @ B0
    # SVD
    U, S, Vt = np.linalg.svd(C)
    # Correct for reflection
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    M = np.diag([1.0, 1.0, d])
    R = Vt.T @ M @ U.T
    return R, cA, cB


# In[ ]:


# produce average structure and align all fragments based on main-chain atoms
first = True
selection = "name CA C O N"
average_name = "average.pdb"
average_file = os.path.join(fragment_templates_dir,average_name)
n = 0

# IGNORE = {"2bcw_B_13_37_1.pdb"}

for i,row in fragment_reps.iterrows():
    fragment_rep_name = row["pdb_file"]
    # if fragment_rep_name in IGNORE:
    #     continue
    fragment_rep_file = os.path.join(fragment_templates_dir,fragment_rep_name)
    u = mda.Universe(fragment_rep_file)
    sel_atoms = u.select_atoms(selection)
    if first:
        all_pos = np.zeros((len(fragment_reps),len(sel_atoms),3))
        all_pos[0,:,:] = sel_atoms.positions.copy()
        first = False
    else:
        sel_pos = sel_atoms.positions.copy()
        R, cA, cB = kabsch_rotate(sel_pos, all_pos[0])    
        all_pos[n,:,:] = (sel_pos - cA) @ R.T + cB
    n += 1
avg_pos = np.mean(all_pos,axis=0)
avg_pos = avg_pos - np.mean(avg_pos,axis=0)
sel_atoms.positions= avg_pos
sel_atoms.write(average_file)

for i,row in fragment_reps.iterrows():
    fragment_rep_name = row["pdb_file"]
    # if fragment_rep_name in IGNORE:
    #     continue
    fragment_rep_file = os.path.join(fragment_templates_dir,fragment_rep_name)
    u = mda.Universe(fragment_rep_file)
    sel_atoms = u.select_atoms(selection)
    sel_pos = sel_atoms.positions.copy()
    all_pos = u.atoms.positions.copy()
    R, cA, cB = kabsch_rotate(sel_pos, avg_pos)
    all_pos = (all_pos - cA) @ R.T + cB
    u.atoms.positions = all_pos
    u.atoms.write(fragment_rep_file)




# In[ ]:
