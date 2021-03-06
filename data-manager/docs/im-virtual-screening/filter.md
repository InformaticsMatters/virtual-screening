# Job: filter

This describes how to run the `filter` job from the `virtual screening` category in the `im-virtual-screening` collection.

## What the job does

This job uses data generated by the [shard job](shard.md) and allows to select a suitable subset of molecules using calculated
molecular properties. Typically used to select candidate molecules for virtual screening,

## Implementation details

* Job implementation: [filter.py](/filter.py)
* Job definition: `jobs.filter` in [im-virtual-screening.yaml](/data-manager/im-virtual-screening.yaml)

## How to run the job

### Inputs

* **Directories with molecules to filter**: the directories containing the sharded data generated by the `shard` job e.g. `molecules/chemspace_feb_2021`.

### Options
Required:
* **Output file name**: name of the output file e.g. `filtered.smi`
* **Minimum Heavy Atom Count**: Minimum number of heavy atoms
* **Maximum Heavy Atom Count**: Maximum number of heavy atoms

Optional:
* **Minimum number of rotatable bonds**, **Maximum number of rotatable bonds**, **Minimum number of rings**, **Maximum number of rings**, **Minimum number of aromatic rings**, **Maximum number of aromatic rings**, **Minimum number of chiral centres**, **Maximum number of chiral centres**, **Minimum number of undefined chiral centres**, **Maximum number of undefined chiral centres**, **Minimum sp3**, **Maximum sp3**: all with integer values of zero or more.

## Related topics

* [Description of the sharded molecule system](https://discourse.squonk.it/t/the-sharded-molecule-system/88)
* [shard](shard.md) job
* [max-min-picker](../rdkit/max-min-picker.md) job