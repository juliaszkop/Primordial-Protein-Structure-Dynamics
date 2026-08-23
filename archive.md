
# **Stability and Dynamics of Primordial Protein Supersecondary Structure Elements**

## **Overview**

This project investigates how **ancient peptides**, which existed before modern folded proteins, could have formed **stable structural motifs**.
The central hypothesis is that the earliest proteins evolved from **short peptides acting as cofactors** for RNA-based catalysis during the RNA world. Initially dependent on RNA scaffolds for structure and activity, these peptides gradually evolved hydrophobic interactions with one another, eventually allowing them to fold **independently of RNA**.

A foundational resource for the work is:
**“A Vocabulary of Ancient Peptides at the Origin of Folded Proteins” – Alva et al., 2015**

## **Ancient Peptide Dictionary**

Researchers previously assembled a “dictionary” of ancient peptides — **40 consensus fragments**, 9–38 residues long.
These fragments:

* were standardized to consensus lengths
* are enriched for functions relevant to the RNA world (binding RNA, nucleotides, metal cofactors)
* often adopt unusual, open conformations, consistent with early RNA-assisted folding

<img width="715" height="573" alt="image" src="https://github.com/user-attachments/assets/24cc3c5c-1251-4167-afa2-16fd30c13d45" />

---

# **Goals of the Project**

A major objective was to identify **structural decoys** — proteins that *look structurally similar* to ancient peptide fragments but **lack evolutionary relatedness**.
These decoys allow testing how special or unique the conformations of ancient peptides are.

**Selection criteria for decoys:**

* **Low RMSD** → structural similarity
* **Low sequence identity** → evolutionary independence

---

# **Methods**

## **1. Local Protein Database Construction**

To enable efficient large-scale structural searches:

* A local PDB database was built using **localpdb**
* Annotated with:

  * **ECOD** classification
  * **cluster30**
  * **DSSP** secondary structure
* These annotations allow grouping by structure and evolutionary family

---

## **2. Structural Matching Using MASTER**

Each ancient peptide fragment was compared against all protein structures using the **MASTER algorithm**:

* Identifies 3D-similar regions regardless of sequence
* All hits below an RMSD threshold were extracted as potential decoys

---

## **3. Filtering Pipeline**

A series of filters was applied to ensure decoys are **truly non-homologous**:

* **ECOD** + **cluster30** → evolutionary filtering
* **BLOSUM62 scoring** → sequence dissimilarity
* Duplicate removal
* Ranking by RMSD

For every peptide, the **5 best structural decoys** (lowest RMSD) were selected for molecular dynamics.

*→ Tutaj wstawiany jest skrypt związany z filtrowaniem.*

---

## **4. Resolution-Based Quality Control**

To guarantee structural accuracy:

* Each candidate PDB structure was checked for **crystallographic resolution**
* If worse than **2.0 Å**, the system searched for a higher-resolution structure of the same sequence
* Ensures minimal noise in MD simulations

*→ Tutaj wstawiany jest drugi skrypt dotyczący kontroli jakości PDB.*

---

## **5. Molecular Dynamics Simulations**

All-atom MD simulations were performed in **GROMACS** under physiological conditions:

* Temperature: **310 K**
* Pressure: **1 atm**
* Salt concentration: **150 mM NaCl**
* Box: **dodecahedral**, periodic
* For each structure: **5 × 10 ns** runs

This produced robust statistics for comparing peptides vs. decoys.

<img width="1307" height="577" alt="image" src="https://github.com/user-attachments/assets/5ba6b259-41d1-43f3-b6bc-baa7f120c247" />


---

# **Descriptors Used in Analysis**

## **Geometry & Dynamics**

* **RMSD per residue**
* **RMS distance**
* **Radius of gyration** (overall and per-axis)

## **Secondary Structure & Composition**

* **DSSP assignments** (helix, sheet, coil)
* **Amino acid frequency**
* **Peptide character** (polar / acidic / non-polar / basic)

## **Intermolecular & Environmental Features**

* **Heavy-atom contacts**
* **Side-chain contacts**
* **Hydrophobic contacts**
* **Hydration field:**

  * local water density
  * free-energy water distribution

*→ Tutaj wstawiany jest skrypt do analizy descriptorów.*

---

# **Results**

## **Descriptor Heatmaps**

The heatmaps show:

* Peptide fragments cluster tightly
* Decoy structures are more spread out
* Indicates conserved structural tendencies among ancient peptides

---

## **Descriptor Correlations**

Pairwise correlation histograms reveal:

* Most descriptor pairs are near **zero correlation**
* Strong positive/negative correlations highlight meaningful relationships

  * e.g., expected: RMSD ↔ heavy-atom contacts (negative)
  * unexpected: hydration ↔ residue composition

These relationships offer clues about ancient structural constraints.

---

## **Principal Component Analysis (PCA)**

* First two components capture **~38%** of variance
* Ancient peptide fragments form a **compact, homogeneous cluster**
* Decoy structures show much higher dispersion

**→ Suggests ancient peptides shared a limited number of stable structural configurations.**

---

# **Repository Structure**

You can structure the repository like this (suggestion):

```
/data
  /pdb_local
  /fragments
  /decoys

/scripts
  filtering_pipeline.py
  resolution_checker.py
  descriptor_analysis.py
  md_preparation.py

/results
  /heatmaps
  /pca
  /correlations

/docs
  figures/
```

---


# **License**

MIT
