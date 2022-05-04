This describes how to run the *align-molecules* job from the *comp chem* category in the *rdkit* collection.

## What the job does

This job takes a template molecule and aligns molecules to it. The matching atoms are determined using Maximum Common Substructure (MCS).

## Implementation details

* Job implementation: [/align_mol.py]()
* Job definition: `jobs.align-molecules` in [/data-manager/rdkit.yaml]()

## How to run the job

### Inputs

* **Reference molecule**: the molecule to align to (.mol or .sdf)
* **Molecules to align**: the molecules to align (.sdf)

### Options

* **Output file name**: name for the output SD-file
* **Reference atoms to align**: atoms in the reference molecule to use for the alignment (1)
* **Keep all alignments**: keep all alignments rather than just the best (e.g. in the case of symmetry related MCS matches)
* **MCS: maximize bonds**:  (2)
* **MCS: match valences**: (2)
* **MCS: ring atom only matches ring atom**:  (2)
* **MCS: match complete rings**:  (2)
* **MCS: match chiral tag**:  (2)
* **MCS: atom compare method**: (2)
* **MCS: bond compare method**: (2)
* **MCS: ring compare method**: (2)


Notes:
(1) Atom numbers as integers in the order of the atoms in the CTAB section of the reference molecule. Specified as a comma separated list of integers. The atom numbers are zero indexed, with the first atom having an index of zero.

(2) The *MCS* parameters correspond to RDKit's options for how to determine a MCS match. See the 
[RDKit docs](http://rdkit.org/docs/source/rdkit.Chem.rdFMCS.html?highlight=findmcs#rdkit.Chem.rdFMCS.FindMCS) for more details.

### Outputs

A SD-file is output with the input molecules aligned to the reference molecule. A *RMSD* property is added that has the RMSD for the alignment.