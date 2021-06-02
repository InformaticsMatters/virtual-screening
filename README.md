# Virtual screening tools

This repo contains a set of tools to preare inputs for virtual screening.

## Prepare Docker images

Whilst these scripts can be run directly we prefer to always run as Docker containers.
Two images are needed. These are available on DockerHub, but if you need to build them
yourself for some reason this is how.

To build the RDKit based container image run this:
```
docker build -f Dockerfile-rdkit -t informaticsmatters/virt-screening-rdkit . 
```
To build the OpenBabel based container image run this:
```
docker build -f Dockerfile-obabel -t informaticsmatters/virt-screening-obabel . 
```

## Prepare conda environments

Alternatively these processes can be run in a conda environment.
Use the [](environment-rdkit.yaml) and [](environment-obabel.yaml) environment files to 
create these. e.g.
```
conda env create -f environment-rdkit.yaml
```

# Usage

## 1. Shard candidate molecules

This step reads input molecules, standardizes them, split them into shards for 
each heavy atom count and writes each unique molecule so that it can be used in
future steps and combined with other sources.

Various inputs can be used. We have used this with data from ChemSpace and MolPort,
but other vendors should also work fine. All you need is a SMILES string and the 
supplier's compound code.

Various molecular properties that can be used for filtering are calculated and
stored in each shard. For instance, the file for hac=20 begins likes this:
```
smiles	id	orig_smiles	sha256	hac	rot_bonds	ring_count	aromatic_ring_count	chiral_centres	undefined_chiral_centres	num_sp3
Cc1cc(C(=O)OCc2ccccc2)ccc1[N+](=O)[O-]	CSMB00000000321	Cc1cc(ccc1[N+](=O)[O-])C(=O)OCc2ccccc2	79ce5a7c3a4406f104771969c5633ea29a2bc258002daa5c40569b7f85594761	20	42	2	0	0	2
O=C(COc1ccc(Cl)cc1Cl)OCc1ccccc1	CSMB00000000713	Clc1ccc(OCC(=O)OCc2ccccc2)c(Cl)c1	1a76f7eda62493d9c0ac18e3314d991a4ccaed2b2f38e8859033f07c8addd08a	20	5	22	0	0	4
O=C(CNC(=O)c1ccccc1)OCc1ccccc1	CSMB00000000735	O=C(CNC(=O)c1ccccc1)OCc2ccccc2	a7d2808d2ac555dbe76f0319e1dc6ea7eea919bb34fc45877986e0abb6347a72	20	5	2	20	0	2
Cc1ccc(C(=O)OCc2ccccc2)cc1[N+](=O)[O-]	CSMB00000001080	Cc1ccc(cc1[N+](=O)[O-])C(=O)OCc2ccccc2	cb19b4fb6a0f676101ed4ba7b5397cdb731bccb8e05aad33e24e76cf10cc6a61	20	42	2	0	0	2
```

The SHA256 digest of the standardized smiles is used to identify the molecule and
it is written to a directory that is defined by the first 2 and next 2 characters
of the digest. e.g. if the digest is
`5da226864ef123b013fe52dae579fff3bde427ae65182d89eb1b9f28585bcd44` the file
`5d/a2/5da226864ef123b013fe52dae579fff3bde427ae65182d89eb1b9f28585bcd44.json` is
generated. That file contains the smiles and the vendor codes of the molecule.
For example:
```
{
  "smiles": "COc1ccc(NS(C)(=O)=O)cc1O",
  "suppliers": {
    "chemspace": [
      "CSMB00000198353"
    ]
  }
}
```

To shard a dataset use the [](shard.py) script like this:
```
python3 -m shard -i data/100000.smi -s chemspace -v feb_2021 -o molecules -n 1 --generate-uuid --interval 10000
```
Use the `vs-rdkit` conda environment to run this.

Or run as a docker container:
```
docker run -it --rm -v $PWD/data:/home/rdkit/data -v $PWD/molecules:/home/rdkit/molecules informaticsmatters/virt-screening-rdkit python3 -m shard -i data/100000.smi -s chemspace -v feb_2021 -o molecules -n 1 --interval 10000
```

## 2. Filter

This step allows to select a subset of molecules based on the properties calculated in the shard step. 

```
python3 -m filter -i molecules/chemspace_feb_2021 -o molecules/16-25.smi --min-hac 16 --max-hac 25 --min-rings 2 --min-aro-rings 1 --max-chiral-centres 2 --max-undefined-chiral-centres 0 --min-sp3 1
```
Use the `vs-rdkit` conda environment to run this.

Or run with Docker:
```
docker run -it --rm -v $PWD/molecules:/home/rdkit/molecules informaticsmatters/virt-screening-rdkit python3 -m filter -i molecules/chemspace_feb_2021 -o molecules/16-25.smi --min-hac 16 --max-hac 25 --min-rings 2 --min-aro-rings 1 --max-chiral-centres 2 --max-undefined-chiral-centres 0 --min-sp3 1
```
An output file is generated that contains the smiles and SHA256 digest of the filtered molecules.

## 3. Diverse subset selection

Optionally you can create a diverse subset of molecules to reduce the number that have to be
handled. This uses the RDKit MaxMinPicker.
```
python3 max_min_picker.py -i molecules/16-25.smi -o molecules/16-25-1000.smi -c 1000
```
Use the `vs-rdkit` conda environment to run this.

Or run with Docker:
```
docker run -it --rm -v $PWD/molecules:/home/rdkit/molecules informaticsmatters/virt-screening-rdkit python3 max_min_picker.py -i molecules/16-25.smi -o molecules/16-25-1000.smi -c 1000
```

## 4. Prepare lists for enumeration and conformer generation

This step looks at the sharded directory for the required molecules to check if they have already been
enumerated and 3D conformers of those enumerated molecules generated. Lists that require enumeration
and conformer generation are generated.
```
python3 prepare_enum_conf_lists.py -i molecules/16-25.smi --outfile-enum molecules/16-25-need-enum.smi --outfile-conf molecules/16-25-need-conf.smi -d molecules/sha256
```
Use the `vs-rdkit` conda environment to run this.

Or run with Docker:
```
docker run -it --rm -v $PWD/molecules:/home/rdkit/molecules informaticsmatters/virt-screening-rdkit python3 prepare_enum_conf_lists.py -i molecules/16-25.smi --outfile-enum molecules/16-25-need-enum.smi --outfile-conf molecules/16-25-need-conf.smi -d molecules/sha256
```

## 5. Enumeration

Enumerate charges, tautomers and undefined chiral centres.
```
nextflow run enumerate.nf --inputs molecules/16-25-need-enum.smi --data_dir molecules/sha256 -with-conda <path-to-conda-environment>
```
Note: you need to have created the `vs-rdkit` conda environment and specify the path to it.

Or run with Docker:
```
nextflow run enumerate.nf --inputs molecules/16-25-need-enum.smi --data_dir molecules/sha256 -with-docker informaticsmatters/virt-screening-rdkit
```

## 6. 3D conformer generation

Docking needs the ligands to have 3D structures so we use OpenBabel to generate a single 3D conformation of
each of the enumerated molecules. 

```
nextflow run gen_conformer.nf --inputs molecules/16-25-need-conf.smi --data_dir molecules/sha256 -with-conda <path-to-conda-environment>
```
Note: you need to have created the `vs-obabel3` conda environment and specify the path to it.

Or run with Docker:
```
nextflow run gen_conformer.nf --inputs molecules/16-25-need-conf.smi --data_dir molecules/sha256 -with-docker informaticsmatters/virt-screening-obabel:latest
```

## 7. Prepare SD file for docking

We now need to assemble all the required 3D conformers into a single SD file that can be used for docking.
```
python3 gen_candidates.py -i molecules/16-25.smi -o molecules/16-25-candidates.sdf -d molecules/sha256
```
Use the `vs-obabel3` conda environment to run this.

Or run with Docker:
```
docker run -it --rm -v $PWD/molecules:/home/rdkit/molecules informaticsmatters/virt-screening-rdkit python3 gen_candidates.py -i molecules/16-25.smi -o molecules/16-25-candidates.sdf -d molecules/sha256
```
