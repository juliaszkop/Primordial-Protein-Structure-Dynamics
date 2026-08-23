#!/usr/bin/env python
# coding: utf-8

import MDAnalysis as mda
import pandas as pd
import numpy as np
import os
import subprocess
import shutil
import argparse
from pathlib import Path
# note: the script assumes that AmberTools are installed!


parser = argparse.ArgumentParser(
    description="Build one simulation directory per MD replica from the extracted "
                "fragment structures."
)
parser.add_argument(
    "fragment_dir",
    help="directory of the fragment to process, e.g. fragment_27",
)
parser.add_argument(
    "-t", "--templates-dir",
    default=None,
    help="directory holding the extracted structures "
         "(default: <fragment_dir>_done, next to the fragment directory)",
)
parser.add_argument(
    "--tool-folder",
    default="ff19_system_maker",
    help="folder with the ff19 files and the make_system_caps script "
         "(default: %(default)s)",
)
args = parser.parse_args()

fragment_dir = os.path.normpath(args.fragment_dir)
if args.templates_dir is not None:
    fragment_templates_dir = os.path.normpath(args.templates_dir)
else:
    fragment_templates_dir = fragment_dir + "_done"


#main dir and necessary files
fragment_reps_file_name = "fragment_templates.csv"
fragment_reps_dir_name = "reps_reference"
fragment_reps_file = os.path.join(fragment_dir,fragment_reps_file_name)
fragment_reps_dir = os.path.join(fragment_dir,fragment_reps_dir_name)
fragment_reps = pd.read_csv(fragment_reps_file,index_col=None)


# procedure to make files for MD runs
# make sure to specify proper location of tool_folder (with ff19 files and scripts)
def make_sim_dir(pdb_file,target_path,tool_folder=args.tool_folder):
    env = os.environ.copy()
    tool_script = "./make_system_caps"
    target_pdb = "peptide.pdb"

    tool_path = Path(tool_folder)
    pdb_tmp = tool_path / target_pdb
    local_dir = tool_path / target_pdb.replace(".pdb","")
    shutil.copy2(pdb_file,pdb_tmp)

    result = subprocess.run(
        [tool_script, target_pdb],
        cwd=tool_path,
        env=env,
        capture_output=True,
        text=True
    )
    output = result.stdout
    if output[-3:-1] == "OK":
        shutil.copytree(local_dir,target_path)
        return 0
    else:
        print("problem processing ",pdb_file)
        return 1


#create empty df for storing the resulting folders
all_reps = pd.DataFrame()
os.makedirs(fragment_reps_dir, exist_ok=True)
parents = fragment_reps["parent"].unique()
#iterate over parent structures (original structures associated with ancestral sequences)
for parent in parents:
    parent_dir = os.path.join(fragment_reps_dir,parent)
    os.makedirs(parent_dir, exist_ok=True)
    structures = fragment_reps.loc[fragment_reps["parent"] == parent]
    #iterate over all structures belonging to given parent
    for i,row in structures.iterrows():
        max_n_dirs = row["counts"]
        source_pdb_name = row["pdb_file"]
        source_pdb = source_pdb_name.replace(".pdb","")
        source_pdb_file = os.path.join(fragment_templates_dir,source_pdb_name)
        target_dir = os.path.join(parent_dir,source_pdb)
        # print(source_pdb_file,target_dir)
        #make simulation dir for MD
        result = make_sim_dir(source_pdb_file,target_dir)

        # if dir made ok, multiply it for all desired copies
        if result == 0:
            first_dir_name = source_pdb+"_1"
            first_dir = os.path.join(parent_dir,first_dir_name)
            shutil.move(target_dir,first_dir)
            all_reps_tmp = pd.DataFrame([[first_dir_name,row["method"],row["resolution"],row["sequence"],parent]],columns=["dir","method","resolution","sequence","parent"])
            all_reps = pd.concat([all_reps,all_reps_tmp])
            for n_dir in range(2,max_n_dirs+1):
                dir_name = source_pdb+"_"+str(n_dir)
                next_dir = os.path.join(parent_dir,dir_name)
                shutil.copytree(first_dir,next_dir)
                all_reps_tmp = pd.DataFrame([[dir_name,row["method"],row["resolution"],row["sequence"],parent]],columns=["dir","method","resolution","sequence","parent"])
                all_reps = pd.concat([all_reps,all_reps_tmp])


#save processed dirs to csv
all_reps_file_name = "all_fragment_reps.csv"
all_reps_file = os.path.join(fragment_dir,all_reps_file_name)
all_reps.to_csv(all_reps_file, index=False)
