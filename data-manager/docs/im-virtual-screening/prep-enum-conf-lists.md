This describes how to run the `prep-enum-conf-lists` job from the `virtual screening` category in the `im-virtual-screening` collection.

## What the job does

This job takes a list of candidate molecules and works out which ones need enumeration and 3D conformer generation. Typically the input for this job is the output of the [filter](https://discourse.squonk.it/t/job-filter/68) or [max-min-picker](https://discourse.squonk.it/t/job-max-min-picker/69) jobs which both generate a tab separated file containing the smiles, the UUID and the sha256 digest. The sha256 is a pointer into the sharded data structure (see the [shard](https://discourse.squonk.it/t/job-shard/67) job). In there will be a `digest.json` file (where digest if the full sha256 digest). When the molecule is enumerated a `digest.smi` is created, and when those enumerated molecules have a 3D conformer generated a `digest.sdf` file is created. This job finds if the `digest.smi` and `digest.sdf` files are present for each member of the input file and, if not, adds that molecule to the list for enumeration or conformer generation.

The output of the job is two files, one for enumeration (specified by the `outputFileEnum` option) and one for conformer generation (specified by the `outputFileConf` option). Those files are inputs for the [enumerate-candidates](https://discourse.squonk.it/t/job-enumerate-candidates/71) and [generate-conformer](https://discourse.squonk.it/t/job-generate-conformer/72) jobs. Note that the `enumerate-candidates` job must be run before the `generate-conformer` job.

## Implementation details

* Job implementation: [/prepare_enum_conf_lists.py]()
* Job definition: `jobs.prep-enum-conf-lists` in [/data-manager/virtual-screening.yaml]()

## How to run the job

### Inputs
* **Molecules to evaluate**: the input molecules to examine
* **Directory with sharded data**: the directory with the sharded data (typically `molecules/sha256`)

### Options
* **Filename for molecules needing enumeration**: the name of the output file for the molecules needing enumeration
* **Filename for molecules needing 3D conformers**: the name of the output file for the molecules needing conformer generation

## Related topics

* [Description of the sharded molecule system](https://discourse.squonk.it/t/the-sharded-molecule-system/88)
* [assemble-conformers job](assemble-conformers.md)