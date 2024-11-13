params.publish_dir = ''
params.publish_dir_mode = 'copy'
params.interval = 10000
params.delimiter = 'tab'
params.name_col = 1
params.skip_lines = 0


process standardize {

    container 'informaticsmatters/vs-moldb:latest'
    if (params.publish_dir) { publishDir params.publish_dir, mode: params.publish_dir_mode }

    input:
    file input

    output:
    file 'std_*'

    """
    python -m moldb.standardize -i $input -o std_$input.name\
      --delimiter $params.delimiter\
      --name-column $params.name_col\
      --skip-lines $params.skip_lines\
      --interval $params.interval
    """
}