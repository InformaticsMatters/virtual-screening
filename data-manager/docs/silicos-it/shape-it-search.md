# Job: shape-it

This topic describes how to run the `shape-it` job from the `virtual screening` category in the `silicos-it` collection.

This performs shape based alignment of a set of molecules against a target molecule.
For more info [look here](http://www.silicos-it.be/software.html#shape-it)


## How to run the job

### Inputs

* **Molecule to align to**:  The molecule to align the candidate to (.mol or .sdf).
* **Molecules to align**: The candidates to align (.sdf)

### Options

* **Output file base name**: Base name for output files (without extension).
* **Score cutoff**: Similarity score cuttoff (value between 0 and 1). Only output molecule more similar than this value.
* **Rank scores by**: Sort the outputs by this metric. Options are TANIMOTO, TVERSKY_DB and TVERSKY_REF.
* **Keep only this many best scores**: output only this many molecules.
* **Do not include reference mol in output**: Do not write reference molecule as the first record of the output.
* **Inputs don't need aligning**: Input molecules are already aligned and just need scoring.
* **Additional simulated annealing steps**: Number of additional simulated annealing steps.


### Outputs

Two files are created based on the *Output file base name* option. A SD-file (.sdf extension) contains the molecules
and their properties and a tab separated values file (.tab extension) contains just the scores.

## Related topics

* [align-it job](align-it-search.md)
* [open3dAlign job](../rdkit/open3dalign.md)