# Fragment merging

## Preparation

Molport and Chemspace had already been sharded in the docking workflow. See that for details of how to prepare these.
This allows us to start by selecting the molecules.

Create `frag-merge` dir for the work:
```
mkdir frag-merge
```

## Select the molecules to select from

We expect to find molecules around 18 heavy atoms, so we choose candidate with 16 - 20 heavy atoms, plus some other filters
tht mean they will be reasonably 'lead like'.

```
./filter.py -i molecules/chemspace_feb_2021 molecules/molport_jan_2021 -o frag-merge/mols.smi --min-hac 16 --max-hac 20 --min-rings 2 --min-aro-rings 1 --max-chiral-centres 2 --min-sp3 1
2021-11-06T09:45:45+00:00 # INFO -EVENT- filter:  Namespace(inputs=['molecules/chemspace_feb_2021', 'molecules/molport_jan_2021'], outfile='frag-merge/mols.smi', min_hac=16, max_hac=20, min_rotb=None, max_rotb=None, min_rings=2, max_rings=None, min_aro_rings=1, max_aro_rings=None, min_chiral_centres=None, max_chiral_centres=2, min_undefined_chiral_centres=None, max_undefined_chiral_centres=None, min_sp3=1, max_sp3=None)
2021-11-06T09:45:45+00:00 # INFO -EVENT- Inspecting files molecules/chemspace_feb_2021/[0-9]*.smi
2021-11-06T09:45:45+00:00 # INFO -EVENT- Inspecting files molecules/molport_jan_2021/[0-9]*.smi
2021-11-06T09:45:45+00:00 # INFO -EVENT- Files to process: ['molecules/chemspace_feb_2021/16.smi', 'molecules/chemspace_feb_2021/17.smi', 'molecules/chemspace_feb_2021/18.smi', 'molecules/chemspace_feb_2021/19.smi', 'molecules/chemspace_feb_2021/20.smi', 'molecules/molport_jan_2021/16.smi', 'molecules/molport_jan_2021/17.smi', 'molecules/molport_jan_2021/18.smi', 'molecules/molport_jan_2021/19.smi', 'molecules/molport_jan_2021/20.smi']
2021-11-06T09:45:45+00:00 # INFO -EVENT- Processing molecules/chemspace_feb_2021/16.smi
2021-11-06T09:45:51+00:00 # INFO -EVENT- Processing molecules/chemspace_feb_2021/17.smi
2021-11-06T09:45:55+00:00 # INFO -EVENT- Processing molecules/chemspace_feb_2021/18.smi
2021-11-06T09:45:59+00:00 # INFO -EVENT- Processing molecules/chemspace_feb_2021/19.smi
2021-11-06T09:46:02+00:00 # INFO -EVENT- Processing molecules/chemspace_feb_2021/20.smi
2021-11-06T09:46:03+00:00 # INFO -EVENT- Processing molecules/molport_jan_2021/16.smi
2021-11-06T09:46:04+00:00 # INFO -EVENT- Processing molecules/molport_jan_2021/17.smi
2021-11-06T09:46:04+00:00 # INFO -EVENT- Processing molecules/molport_jan_2021/18.smi
2021-11-06T09:46:04+00:00 # INFO -EVENT- Processing molecules/molport_jan_2021/19.smi
2021-11-06T09:46:05+00:00 # INFO -EVENT- Processing molecules/molport_jan_2021/20.smi
2021-11-06T09:46:07+00:00 # INFO -EVENT- Matched 4566242 out of 7082729 records. 53437 duplicates
```


## Find similar molecules

```
./screen.py --smiles 'CC(=O)NC1=CN=CC=C1C' 'CC(NC(=O)C)c1cccc(Cl)c1' --input frag-merge/mols.smi --output frag-merge/screened.smi --metric tversky --descriptor morgan2 --interval 10000 --sim-index 5 --threshold 0.3
...
...
2021-11-06T13:47:09+00:00 # INFO -EVENT- Processed 4560000 records, 4362 hits
2021-11-06T13:47:10+00:00 # INFO -EVENT- Inputs: 4566242 Hits: 4386 Errors: 17 Time (s): 632.3182473182678
```

## Prepare conformers

```
./prepare_enum_conf_lists.py -i frag-merge/screened.smi --outfile-enum frag-merge/need-enum.smi --outfile-confs frag-merge/need-confs.smi
prepare_enum_conf_lists.py:  Namespace(input='frag-merge/screened.smi', outfile_enum='frag-merge/need-enum.smi', outfile_confs='frag-merge/need-confs.smi', data_dir='molecules/sha256')
2021-11-06T17:06:48+00:00 # INFO -EVENT- Processing file frag-merge/screened.smi
2021-11-06T17:06:48+00:00 # INFO -EVENT- Processed 4386 records. 0 duplicates, 0 errors. 1396 already enumerated, 1396 already have low energy conformers, 2990 need enumeration, 2990 need low energy conformers generated
```

```
nextflow -log frag-merge/.nextflow.log run enumerate.nf --inputs frag-merge/need-enum.smi -with-docker informaticsmatters/vs-prep:latest -with-trace frag-merge/trace-enum.txt -with-report frag-merge/report-enum.html
N E X T F L O W  ~  version 20.10.0
Launching `enumerate.nf` [mad_fermat] - revision: 638c75eb4d
executor >  local (4)
[9a/6511da] process > splitter      [100%] 1 of 1 ✔
[2a/3ca63b] process > enumerate (3) [100%] 3 of 3 ✔
Completed at: 06-Nov-2021 17:36:47
Duration    : 25m
CPU hours   : 0.9
Succeeded   : 4

```


```
nextflow -log frag-merge/.nextflow.log run le_conformers.nf --inputs frag-merge/need-confs.smi -with-docker informaticsmatters/vs-prep:latest -with-trace frag-merge/trace-confs.txt -with-report frag-merge/report-confs.html
N E X T F L O W  ~  version 20.10.0
Launching `le_conformers.nf` [festering_kowalevski] - revision: 502d68b8aa
executor >  local (31)
[ab/4b21f9] process > splitter            [100%] 1 of 1 ✔
[07/250c1e] process > gen_conformers (27) [100%] 30 of 30 ✔
Completed at: 06-Nov-2021 21:39:28
Duration    : 3h 59m 58s
CPU hours   : 68.0
Succeeded   : 31
```


```
./assemble_conformers.py -i frag-merge/screened.smi -m low-energy -o frag-merge/conformers-merge.sdf.gz --interval 10000
2021-11-07T07:53:44+00:00 # INFO -EVENT- assemble_conformers.py:  Namespace(input='frag-merge/screened.smi', output='frag-merge/conformers-merge.sdf.gz', data_dir='molecules/sha256', mode='low-energy', exclude_base=False, exclude_tautomers=False, exclude_microstates=False, interval=10000)
2021-11-07T07:59:36+00:00 # INFO -EVENT- Processed 4386 inputs with 4386 unique mols in 351.3297746181488s. 0 errors, 0 duplicates
```

## Shape screening

USRCAT using the full merged fragments

```
./usr.py -i $D/conformers-merge.sdf.gz -q notebooks/merged1.mol -o $D/results-merge1-usrcat.sdf -t 0.6 -m usrcat -g std_smi --interval 10000
2021-11-07T08:08:30+00:00 # INFO -EVENT- usr.py:  Namespace(group_by_field='std_smi', inputs='frag-merge/conformers-merge.sdf.gz', interval=10000, method='usrcat', outfile='frag-merge/results-merge-usrcat.sdf', query='merged.sdf', threshold=0.6)
2021-11-07T08:08:30+00:00 # INFO -EVENT- Opening frag-merge/results-merge1-usrcat.sdf as output
...
...
...
2021-11-07T09:45:38+00:00 # INFO -EVENT- Processed 1331721 conformers. Generated 76 outputs. 0 errors. Average similarity is 0.32917631516010165
```

Electroshape using the full merged fragments
```
./usr.py -i $D/conformers-merge.sdf.gz -q notebooks/merged1.mol -o $D/results-merge1-electroshape.sdf -t 0.9 -m electroshape -g std_smi --interval 10000
2021-11-07T11:34:17+00:00 # INFO -EVENT- usr.py:  Namespace(group_by_field='std_smi', inputs='frag-merge/conformers-merge.sdf.gz', interval=10000, method='electroshape', outfile='frag-merge/results-merge-electroshape.sdf', query='merged.sdf', threshold=0.9)
...
...
...
2021-11-07T12:42:45+00:00 # INFO -EVENT- Processed 1331721 conformers. Generated 23 outputs. 0 errors. Average similarity is 0.6508460694033357
```

USRCAT using the truncated Mpro-x1382 fragment
```
./usr.py -i $D/conformers-merge.sdf.gz -q notebooks/merged2.mol -o $D/results-merge2-usrcat.sdf -t 0.6 -m usrcat -g std_smi --interval 10000
2021-11-07T13:29:28+00:00 # INFO -EVENT- usr.py:  Namespace(group_by_field='std_smi', inputs='frag-merge/conformers-merge.sdf.gz', interval=10000, method='usrcat', outfile='frag-merge/results-merge2-usrcat.sdf', query='notebooks/merged2.sdf', threshold=0.6)
...
...
...
2021-11-07T15:07:21+00:00 # INFO -EVENT- Processed 1331721 conformers. Generated 235 outputs. 0 errors. Average similarity is 0.32401849441262187

```

Electroshape using the truncated Mpro-x1382 fragment
```
./usr.py -i $D/conformers-merge.sdf.gz -q notebooks/merged2.mol -o $D/results-merge2-electroshape.sdf -t 0.9 -m electroshape -g std_smi --interval 10000
...
...
...
2021-11-07T16:15:49+00:00 # INFO -EVENT- Processed 1331721 conformers. Generated 62 outputs. 0 errors. Average similarity is 0.6521101257849915
```
