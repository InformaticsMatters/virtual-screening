This describes how to run the *open3dalign* job from the *comp chem* category in the *rdkit* collection.

## What the job does

Open3DAlign is a molecular alignment tool that can perform a rigid alignment of conformers to a reference structure. The original tool can be found [here](http://open3dalign.sourceforge.net/) and is described [here](https://doi.org/10.1007/s10822-011-9462-9), though the version used here is the one found in [RDKit](http://rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html?highlight=open3d). TODO - explain the differences.

It takes as inputs a set of 3D conformers to be aligned and the reference molecule to align them to. The output is the conformers aligned to the reference molecule. The outputs can be filtered according to the quality of the alignment.

Typically this can be used as part of a ligand based virtual screening campaign to identify molecules with similar 3D shape and molecular properties to a known active molecule.

## Implementation details

Job implementation: [/open3dalign.py]()
Job definition: `jobs.open3dalign` in [/data-manager/rdkit.yaml]()

## How to run the job

### Inputs

**Molecules to screen**: A SD file containing the 3D conformers to screen.
**Query Molecule**: The reference molecule to which the conformers should be aligned.

### Options

**Output file name**: The name of the output SD file.
**Use Crippen contributions**: Include Crippen lipophilicity criteria
**Remove Hs from outputs**: Remove explicit hydrogens from the output
**Filter threshold**: Only include molecules above this threshold (1)

Notes:
(1) Filtering uses the `o3da_score_rel` field which is the relative alignment score. The alignment of the reference molecule to itself would give a score of 1. Any other alignment will give a score less than that. Values should be between 0 and 1. If no value is given then all molecules will be output.

### Outputs

A SD file containing the aligned inputs, optionally filtered by the *Filter threshold* option.

## Related topics

* [Description of the sharded molecule system](https://discourse.squonk.it/t/the-sharded-molecule-system/88)
* [prep-enum-conf-lists job](../im-virtual-screening/prep-enum-conf-lists.md)
* [generate-low-energy-conformers job](generate-low-energy-conformers.md)
* [assemble-conformers job](../im-virtual-screening/job-assemble-conformers.md)
* [align-it job](../silicos-it/align-it.md)
* [shape-it job](../silicos-it/shape-it.md)
* [ultrafast-shape-recognition](../im-virtual-screening/ultrafast-shape-recognition.md)