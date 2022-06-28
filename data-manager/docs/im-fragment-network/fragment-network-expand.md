# Job: fragment-network-expand

This describes how to run the `fragment-network-expand` job from the `virtual screening` category in the `im-fragment-network` collection.

This uses the fragment network to find molecules related to a query molecule. The results are those molecules within a
specified number of hops from the query and some optional molecular filters such as the change in heavy atom or ring atom
count.


## Implementation details

* Job implementation: [fn_synthon_expansion.py](/fn_synthon_expansion.py)
* Job definition: `jobs.fragment-network-synthon-expansion` in [fragnet-search.yaml](/data-manager/fragnet-search.yaml)

## How to run the job

### Inputs

This job has no file inputs.

### Options

* **Query SMILES**: The SMILES of the molecule to start from.
* **Output file**: The name of the output file e.g. `expansions.smi`.
* **Number of fragment network hops**: The number of edges to traverse in the fragment network. Max allowed is 3.
* **Min heavy atom count**: Permitted reduction in heavy atom count for resulting molecules.
* **Max heavy atom count**: Permitted increase in heavy atom count for resulting molecules.
* **Min ring atom count**: Permitted reduction in ring atom count for resulting molecules.
* **Max ring atom count**: Permitted increase in ring atom count for resulting molecules.

## Related topics

* [About the fragment network](https://squonk.it/fragnet-search-ui)
* [fragment-network-find-synthons](fragment-network-find-synthons.md) job
* [fragment-network-synthon-expansion](fragment-network-synthon-expansion.md) job