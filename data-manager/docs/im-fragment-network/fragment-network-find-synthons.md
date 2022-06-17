# Job: fragment-network-find-synthons

This describes how to run the `fragment-network-find-synthons` job from the `virtual screening` category in the `im-fragment-network` collection.

This uses the fragment network to generate "synthons" of a query molecule.
A "synthon" is a fragment of the molecule that could potentially be encountered in other molecules.

The resulting synthons can be used in the [fragment-network-synthon-expansion](fragment-network-synthon-expansion.md) job.

## Implementation details

* Job implementation: [fn_find_synthons.py](/fn_find_synthons.py)
* Job definition: `jobs.fragment-network-find-synthons` in [fragnet-search.yaml](/data-manager/fragnet-search.yaml)

## How to run the job

### Inputs

This job has no file inputs.

### Options

* **Query SMILES**: The SMILES of the molecule to start from.
* **Output file**: The name of the output file e.g. `synthons.smi`.

## Related topics

* [About the fragment network](https://squonk.it/fragnet-search-ui)
* [fragment-network-expand](fragment-network-expand.md) job
* [fragment-network-synthon-expansion](fragment-network-synthon-expansion.md) job
