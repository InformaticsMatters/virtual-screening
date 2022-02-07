params.interval = 1000
params.num_charges = 2
params.try_embedding = true
params.add_hydrogens = true
params.max_tautomers = 25

process enumerate {

    container 'informaticsmatters/vs-prep:latest'

    input:
    file inputs
    file data_dir

    """
    /code/enumerate.py -i $inputs --data-dir $data_dir --interval $params.interval\
      --enumerate-tautomers --enumerate-chirals --enumerate-charges\
      ${params.num_charges ? '--num-charges ' + params.num_charges : ''}\
      ${params.try_embedding ? '--try-embedding' : ''}\
      ${params.add_hydrogens ? '--add-hydrogens' : ''}\
      --max-tautomers $params.max_tautomers
    """
}