# Virtual screening tools

This repo contains a set of tools to prepare inputs for virtual screening.

## Prepare Docker images

Whilst these scripts can be run directly we prefer to always run as Docker containers.
Two images are needed. These are available on DockerHub, but if you need to build them
yourself for some reason this is how...

To build the container images run this:

    $ IMAGE_TAG=1.0.0 docker-compose build 

Or, to build using the `latest` tag: -

    $ docker-compose build 

## Prepare conda environments

Alternatively these processes can be run in a conda environment.
Use the [](environment.yaml) environment files to 
create a conda environment named `im-vs-prep`. e.g.
```
conda env create -f environment.yaml
conda activate im-vs-prep
```
This environment contains all the tools needed to run these processes.

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
./shard.py -i data/100000.smi -s chemspace -v feb_2021 -o molecules -n 1 --interval 10000
```
Use the `im-vs-prep` conda environment to run this.

Or run as a docker container:
```
docker run -it --rm -v $PWD:$PWD -w $PWD -u 1000:1000\
  informaticsmatters/vs-prep:$IMAGE_TAG \
  /code/shard.py -i data/100000.smi -s chemspace -v feb_2021 -o molecules -n 1 --interval 10000
```

This runs the container as your user ID rather than the default of `root` so that the created files
have the correct ownnership.
Update the `-u 1000:1000` bit with your actual user and group ID).
You will need to create the `molecules` dir before running this so that it has the right ownership.

## 2. Filter

This step allows to select a subset of molecules based on the properties calculated in the shard step.

```
./filter.py -i molecules/chemspace_feb_2021 -o 16-25.smi --min-hac 16 --max-hac 25 --min-rings 2 --min-aro-rings 1 --max-chiral-centres 2 --max-undefined-chiral-centres 0 --min-sp3 1
```
Use the `im-vs-prep` conda environment to run this.

Or run with Docker:
```
docker run -it --rm -v $PWD:$PWD -w $PWD -u 1000:1000 informaticsmatters/vs-prep:$IMAGE_TAG\
  /code/filter.py -i molecules/chemspace_feb_2021 -o 16-25.smi --min-hac 16 --max-hac 25 --min-rings 2\
  --min-aro-rings 1 --max-chiral-centres 2 --max-undefined-chiral-centres 0 --min-sp3  1
```
An output file is generated that contains the smiles and SHA256 digest of the filtered molecules.

## 3. Diverse subset selection

Optionally you can create a diverse subset of molecules to reduce the number that have to be
handled. This uses the RDKit MaxMinPicker.
```
./max_min_picker.py -i 16-25.smi -o 16-25-1000.smi -c 1000
```
Use the `vs-rdkit` conda environment to run this.

Or run with Docker:
```
docker run -it --rm -v $PWD:$PWD -w $PWD -u 1000:1000 informaticsmatters/vs-prep:$IMAGE_TAG\
  /code/max_min_picker.py -i 16-25.smi -o 16-25-1000.smi -c 1000
```

## 4. Prepare lists for enumeration and conformer generation

This step looks at the sharded directory for the required molecules to check if they have already been
enumerated and 3D conformers of those enumerated molecules generated. Lists that require enumeration
and conformer generation are generated.
```
./prepare_enum_conf_lists.py -i 16-25.smi --outfile-enum need-enum.smi --outfile-conf need-conf.smi -d molecules/sha256
```
Use the `im-vs-prep` conda environment to run this.

Or run with Docker:
```
docker run -it --rm -v $PWD:$PWD -w $PWD -u 1000:1000 informaticsmatters/vs-prep:$IMAGE_TAG\
  /code/prepare_enum_conf_lists.py -i 16-25.smi\
  --outfile-enum need-enum.smi --outfile-conf need-conf.smi -d molecules/sha256
```

## 5. Enumeration

Enumerate charges, tautomers and undefined chiral centres.
```
nextflow run enumerate.nf --inputs need-enum.smi --data_dir molecules/sha256 -with-conda <path-to-conda-environment>
```
Note: you need to have created the `vs-rdkit` conda environment and specify the path to it.

Or run with Docker:
```
nextflow run enumerate.nf --inputs need-enum.smi --data_dir molecules/sha256 -with-docker informaticsmatters/vs-prep:$IMAGE_TAG
```

## 6. 3D conformer generation

Docking needs the ligands to have 3D structures so we use OpenBabel to generate a single 3D conformation of
each of the enumerated molecules. 

```
nextflow run gen_conformer.nf --inputs need-conf.smi --data_dir molecules/sha256 -with-conda <path-to-conda-environment>
```
Note: you need to have created the `vs-obabel3` conda environment and specify the path to it.

Or run with Docker:
```
nextflow run gen_conformer.nf --inputs need-conf.smi --data_dir molecules/sha256 -with-docker informaticsmatters/vs-prep:$IMAGE_TAG
```

## 7. Prepare SD file for docking

We now need to assemble all the required 3D conformers into a single SD file that can be used for docking.
```
./gen_candidates.py -i 16-25-1000.smi -o 16-25-candidates.sdf -d molecules/sha256
```
Use the `vs-obabel3` conda environment to run this.

Or run with Docker:
```
docker run -it --rm -v $PWD:$PWD -w $PWD -u 1000:1000 informaticsmatters/vs-prep:$IMAGE_TAG\
  /code/gen_candidates.py -i 16-25-1000.smi -o 16-25-candidates.sdf -d molecules/sha256
```

## 8. Prepare PDB file for docking using OpenBabel

Docking tools need the protein to be properly prepared, wich is quite a complex topic.
rDock needs the protein to be input in MOL2 format, whereas smina needs it in PDBQT format,
but the typical starting point is a PDB file.
We'll use OpenBabel to do the conversion, and to prototnate the protein at a suitable pH.
Thisis a relatively simple and crude approach.
rDock needs polar hydrogens to be present, but doesn't care if non-polar ones are present or not.
```
obabel data/dhfr-receptor.pdb -O dhfr-receptor-ph7.mol2 -p 7
```
Use the `vs-obabel3` conda environment to run this.
Similar commands can be used to generate PDBQT format for input to smina.

Or run with Docker:
```
docker run -it --rm -v $PWD:$PWD -w $PWD -u 1000:1000 informaticsmatters/vs-prep:$IMAGE_TAG\
  obabel data/dhfr-receptor.pdb -O dhfr-receptor-ph7.mol2 -p 7
```

## 9. Prepare PDB file for docking using pdb2pqr

OpenBabel only does a very basic preparation of the protein, protonating at a specific pH but
not taking the precise environment of the titratable groups into consideration.
A more sophisticated approach can be performed using [pdb2pqr](https://github.com/Electrostatics/pdb2pqr)
which can correct for missing heavy atoms and, using [propka](https://github.com/jensengroup/propka)
can predict ionisation states including local information such as hydrogen bonding networks.
```
pdb2pqr30 --with-ph 7.0 --titration-state-method propka --pdb-output dhfr-receptor-ph7.pdb\
  data/dhfr-receptor.pdb dhfr-receptor-ph7.pqr
```
Use the `im-vs-prep` conda environment to run this.

Or run with Docker:
```
docker run -it --rm -v $PWD:$PWD -w $PWD -u 1000:1000 informaticsmatters/vs-prep:$IMAGE_TAG\
  pdb2pqr30 --with-ph 7.0 --titration-state-method propka\
  --pdb-output dhfr-receptor-ph7.pdb\
  data/dhfr-receptor.pdb dhfr-receptor-ph7.pqr
```

Note: if using pdb2pqr you still need to convert to MOL2 and/or PDBQT formats for docking
using OpenBabel as above, but don't specify the `-p 7` option so that OpenBabel does not
protonate the protein.

## 10. Docking with rDock

We use the [](rdock-docking.nf) Nextflow workflow for performing the docking.
This workflow splits the input SDF into multiple chunks that can be docked in parallel, and then
collates the reulst into a single output and does some very basic analysis.

Before we can do this we must prepare a rDock configuration file (.prm file) and generate the 
cavity defintion (.as file). The [](prepare_rdock.py) script can be used to do this. Currently 
we don't have a conda environment for this so it must be done with Docker. Berfore running copy the 
`data/dhfr-ligand.mol` file to this directory:
```
cp data/dhfr-ligand.mol .
```
Now create the rDock configuration:
```
docker run -it --rm -v $PWD/:$PWD -w $PWD -u 1000:1000 informaticsmatters/vs-rdock:$IMAGE_TAG\
  /code/prepare_rdock.py --receptor dhfr-receptor-ph7.mol2 --ligand dhfr-ligand.mol --output docking
```
This creates a basic `docking.prm` file with the rDock configuration and the `docking.as` file with
the active site definition. You can edit the `docking.prm` file if you need a different configuration.

Now we can run the docking.
```
nextflow run rdock-docking.nf\
  --ligands 16-25-candidates.sdf\
  --protein dhfr-receptor-ph7.mol2\
  --prmfile docking.prm\
  --asfile docking.as\
  --num_dockings 5
```
We use `--num_dockings 5` to speed things up. A real run would typically use 25-50 dockings for each candidate.


## 11. Docking with smina

Alternatively we can use [smina](https://sourceforge.net/projects/smina/), an
[AutoDock VINA](http://vina.scripps.edu/) derivative that is faster and easier to use.

Like rDock it can handle inputs in SDF format, but needs the receptor in PDBQT fromat.
It can generate the box configuration (partly equivalent of the rdock cavity definition)
on the fly if you specify a sample ligand.

Prepare the PDBQT file as we did for the MOL2 file needed for rDock:
```
docker run -it --rm -v $PWD:$PWD -w $PWD -u 1000:1000 informaticsmatters/vs-prep:$IMAGE_TAG\
  obabel data/dhfr-receptor.pdb -O dhfr-receptor-ph7.pdbqt -p 7
```
And prepare the candidate ligand:
```
docker run -it --rm -v $PWD:$PWD -w $PWD -u 1000:1000 informaticsmatters/vs-prep:$IMAGE_TAG\
  obabel data/dhfr-ligand.mol -h -O dhfr-ligand.pdbqt
```
Now we can run the docking:
```
nextflow run smina-docking.nf\
      --ligands 16-25-candidates.sdf\
      --protein dhfr-receptor-ph7.pdbqt\
      --ligand dhfr-ligand.pdbqt\
      --padding 4\
      --exhaustiveness 8\
      --scoring_function default
```
We only expose a few of the key smina parameters. For the `scoring_function` you could specify 
ad4_scoring, dkoes_fast, dkoes_scoring, dkoes_scoring_old, vina or vinardo instead of default.

Note that the workflow is quite smart in that you can pass in the receptor and ligand in other
formats and they will be converted to PDBQT format automatically.


## 12. Rescoring and interactions using ODDT
[ODDT](https://github.com/oddt/oddt) provides a wide range of tools that are useful in virtual screening.
We use it to run some more sophisticated scoring functions on the docked poses and to calculate the
interactions the poses make with the receptor to help in prioritising the docked poses.

To run this on the rDock results:
```
nextflow run oddt.nf --ligands results/results_rdock.sdf --protein data/dhfr-receptor.pdb
```
This generates the `results/results_oddt.sdf` output. This workflow can also be run on the
smina poses.



