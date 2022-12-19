# Job: ph4-align-to-fragments

This describes how to run the `ph4-align-to-fragments-smiles`, `ph4-align-to-fragments-simple` and
`ph4-align-to-fragments-workflow` jobs from the `ligand based virtual screening` category in the `im-virtual-screening`
collection.

## What the job does

These job runs the PLANTS in align mode to perform flexible alignment of ligands to one or more fragment molecules.
PLANTS is primarily a docking algorithm, but with the `--mode align` option it can use the same Ant Colony Optimistaion
methodology to perform flexible alignment of a molecule to a fixed or flexible target.

This job performs flexible alignment to the fixed 3D structures of one or more fragment molecules.
The prime aim it to try to identify larger molecules that combine features of (typically) two observed fragment poses.

The 3 variants are:

* **ph4-align-to-fragments-smiles**: molecules to be aligned specified as an option in SMILES format. Suitable for a handful on inputs.
* **ph4-align-to-fragments-simple**: reads molecules to be aligned from a file in SDF or SMILES format. Suitable for 100's or 1000's of inputs.
* **ph4-align-to-fragments-workflow**: parallelised workflow that also re-scores using Open3DAlign. Suitable for 100's of thousands of inputs.

Whilst PLANTS does quite a good job of generating interesting poses, the score it generates for these poses does not 
necessarily rank poses in what would seem to be the optimum order. To address this, the workflow job adds an additional
alignment step using Open3DAlign (implemented as part of RDKit), which does a rigid re-alignment of the molecule to the 
fragments and generates scores that appear to be more meaningful.

The complete workflow job does this:

1. Split the SD-file with the candidate molecules in chunks of 5000 molecules so that each chunk can be processed in parallel.
2. Perform flexible alignment using PLANTS. The fragments are treated as fixed and the candidate molecule as flexible.
   The number of conformers generated is determined by the *count* option.
3. Performed rigid alignment of the generated conformer to the fragments using Open3D align.
4. Combine the results for all chunks back into a single SD-file.
5. As long as a *Field to group molecules* option is specified the best aligned molecule in each group (using the
   *o3da_score_rel* field) is determined and then those are sorted using that score.

When handling alignment to multiple fragments, the PLANTS stage handles this as alignment to 2 separate fixed molecules,
whilst Open3DAlign can only align to a single molecule so the multiple fragments are merged into a single molecule
(potentially with overlapping atoms).

## Implementation details

The workflow job is implemented as a [Nextflow](https://www.nextflow.io/) workflow.

* Nextflow workflow: [frag-merge-pharmacophore.nf](/frag-merge-pharmacophore.nf)
* Job definition: `jobs.ph4-align-to-fragments-*` in [im-virtual-screening.yaml](/data-manager/im-virtual-screening.yaml)
* About PLANTS: https://pubs.acs.org/doi/10.1021/ci1000218
* About Open3DAlign: https://open3dalign.sourceforge.net/

## How to run the job

### Inputs

* **Candidates to process**: the molecules to dock in SDF or SMILES format.
* **Fragments to align to**: the fragment molecules, 3D molecules in SDF format. Typically 2 molecules, but could be 1 or 3.

### Options

For all 3 jobs:

* **Number of mols to generate**: the number of 3D structures PLANTS should try to generate. Should be no need for more than 10.
* **PLANTS Torsion weight setting**: lower values allows to generate better aligned molecules but with less "correct"
    torsion angles (1).
* **PLANTS RMSD**: RMSD for discarding similar conformers (2)
* **PLANTS cuttoff score**: Only keep molecules with score better (e.g. lower numbers) than this cutoff. e.g. -50

Additional options for the *smiles* and *simple* jobs:

* **Filename for results**: the file name for the results (can include a directory path)

Additional options for the *workflow* job:

* **Directory for results**: directory for the results. If not defined the top level directory is used, which is not a good idea.
* **Filename for results**: the file name for the results (do not include a directory path).
* **Field to group molecules**: the field name that identifies the consecutive groups of molecules (3).
* **Open3DAlign use Crippen**: Open3DAlign should use Crippen logP atom contributions instead of MMFF atom types and charges
* **Open3DAlign cuttoff score**: Only keep molecules with score better (e.g. higher numbers) than this cutoff. e.g. 0.4

Additional options for the *workflow* and *simple* jobs:

* **Generate 3D coords**: generate 3D coordinates for the input molecules (using OpenBabel). This should be specified if
  the inputs are not already 3D structures e.g. when using SMILES as input (coordinates are always generated for the 
  *smiles* job)

Notes:

1. PLANTS's default value is 20, but this job uses a default of 5 to allow to generate better aligned structures at the cost
   of slightly less ideal torsion angles.
2. A default value of 2 is used.
3. This should be a field in the inputs that uniquely identifies each molecule. If PLANTS is told to generate 10 molecules
   (the **count** option) then there will be 10 consecutive molecules in the PLANTS output with that value. Whn the value
   changes a new group commences. If it's the title line (first line in the record) that identifies the molecule then use
   the value *_Name*.

### Outputs

### Files

Two files can be output. You always get the file specified by the *Filename for results* option containing the entire
set of generated molecules and scores. If you specify a value for the *Field to group molecules* then you also get a 
file that has the best aligned structure in the group (determined from the *o3da_score_rel* field) and then have those
remaining structures sorted so that the first molecule in the file has the best alignment score. If the output file is
specified as `results.sdf` then this second file will be called `results-best.sdf`.

### Fields

The following fields are added by this job (all fields in the input will aslo be retained):

* **PH4_SCORE**: The PLANTS alignment score, a pseudo delta G value (the more negative the better)
* **o3da_score**: The Open3DAlign score. The numbers are largely meaningless by the bigger the value the better the alignment.
* **o3da_score_rel**: The Open3DAlign score relative to the score for aligning the fragments to themselves. Value between
  0 and 1, with one being a perfect alignment.
* **o3da_align**: Open3DAlign alignment score (unclear what this is)