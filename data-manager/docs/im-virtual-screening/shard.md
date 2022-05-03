This describes how to run the `shard` job from the `virtual screening` category in the `im-virtual-screening` collection.

## What the job does

The job takes molecules in smiles format as input and:
* standardizes the molecules
* calculates a series of molecular properties of each molecule
* generates a sha256 digest of the standardized smiles
* shards them into a series of shared directories based on the first 4 characters of the sha256 digest (1)
* shards them into series of files defined by the heavy atom count of the molecule, also writing the molecular properties to the file (2)

Note (1).  The resulting directory structure will be up to 256 top level directories based on the first 2 characters of the digest, and in each of those top level directories up to 256 second level directories based on characters 3 and 4 of the digest. In those 2nd level directories are files with a base name of the digest  and a `.json` extension. Those files contain the standardized smiles of the molecule, and UUID that is used to identify the molecule and information about the suppliers of that molecule. When a second dataset (e.g. from a different supplier) is sharded into the same directory system the information for that supplier is appended to the existing information.  An example file looks like this:
```
{
  "smiles": "CN(CCNC(=O)OC(C)(C)C)C(=O)C1CC1(Cl)Cl",
  "uuid": "8f4579c7-b7ad-4975-9a31-b0f871cd0f35",
  "suppliers": {
    "chemspace": [
      "CSMB00114229931"
    ]
  }
}
```

Note (2). A directory is created for that dataset (supplier) and in that directory will be files named `*.smi` where the base name is the number of heavy atoms for the molecule. e.g. all molecules with 22 heavy atoms will be in a  file named 22.smi. Molecular properties and other information is written to the tile (tab separated text) which looks like this:
```
smiles	uuid	id	orig_smiles	sha256	hac	rot_bonds	ring_count	aromatic_ring_count	chiral_centres	undefined_chiral_centres	num_sp3
c1ccc(COc2ccc(OCc3ccccc3)cc2)cc1	47110a97-5e79-4715-983a-22e7eaa1f873	CSMB00000000337	C(Oc1ccc(OCc2ccccc2)cc1)c3ccccc3	7609dafe9ffe71cbb76d4061fac621e09091466708904d7dc0db9735037f8058	22	6	3	3	0	0	2
O=C(c1ccncc1)N(Cc1ccccc1)c1ccccc1	97c64582-06dc-483d-9344-dec058e0ce24	CSMB00000001310	O=C(N(Cc1ccccc1)c2ccccc2)c3ccncc3	5f9dba93286c4aba8e0fce67f5253934eae4e4a46cb6eec9cc473b4f77c7b9a1	22	4	3	3	0	0	1
CCCNC(=O)C(NC(=O)OCc1ccccc1)C(C)CC	b3126958-1952-4e6f-acfc-52464053997b	CSMB00000001328	CCCNC(=O)C(NC(=O)OCc1ccccc1)C(C)CC	397b8859e77c7c3a9f067ff261552aed3802890e24c63a5ad4708b7a6d110c38	22	8	1	1	2	2	9
CCC(C)C(NC(=O)OCc1ccccc1)C(=O)NC(C)C	9ff36477-18e3-422b-8715-d8de78a60df0	CSMB00000001412	CCC(C)C(NC(=O)OCc1ccccc1)C(=O)NC(C)C	b837e9b87fcd19b8f7c534aa090a46b20749de0b1aba6ad5cc7b26a71df41a19	22	7	1	1	2	2	9
```

The main purpose of running the shard job is to generate data that can easily filtered to select molecules of interest. See the `filter` job that can be used to do this.

## Implementation details

Job implementation: [/shard.py]()
Job definition: `jobs.shard` in [/data-manager/virtual-screening.yaml]()

## How to run the job

### Inputs
**Molecules to shard**: one or more files from a vendor in SMILES format. An example file with a small subset of ChemSpace can be found [here](https://github.com/InformaticsMatters/virtual-screening/blob/main/data/100000.smi).
### Options
**Source of the molecules**: a name that describes the source of the molecules e.g. `ChemSpace`
**Version of the molecule source**: something that identifies the version of the molecules such as a release number or date e.g.  `feb_2021`
**Index of vendor code field**: the zero based index of the column in the input that contains the vendor's code for that molecule e.g. `1`
**Output directory**: the location where the sharded data is written. By convention use the value `molecules` unless you are just testing and do not want to modify that location.

## Limitations
1. It is currently assumed that the input is tab separated. The python module allows the delimited to be specified, but this option is not yet exposed in the job.
2. The same applies to being able to skip a header line which is not currently possible.

## Related topics

* [Description of the sharded molecule system](https://discourse.squonk.it/t/the-sharded-molecule-system/88)
* [filter job](filter.md)