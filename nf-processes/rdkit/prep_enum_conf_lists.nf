
params.enum_out = 'need-enum.smi'
params.conf_out = 'need-confs.smi'
params.interval = 10000

process prep_lists {

    container 'informaticsmatters/vs-prep:latest'

    input:
    path inputs // .smi
    path data_dir  // e.g. molecules

    output:
    path params.enum_out
    path params.conf_out

    """
    /code/prepare_enum_conf_lists.py\
      --input $inputs\
      --data-dir $data_dir\
      --outfile-enum $params.enum_out\
      --outfile-confs $params.conf_out\
      --interval $params.interval
    """
}