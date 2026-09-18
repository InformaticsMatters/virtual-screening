params.interval = 1000
params.purge = true

params.library_name = 'no_name'

process load_standardized {

    container 'informaticsmatters/vs-moldb:2.0.0'
    maxForks 1
    errorStrategy { sleep(Math.pow(2, task.attempt) * 1000 as long); return 'retry' }
    maxRetries 3

    input:
    path inputs

    script:
    """
    python -m moldb.load_standardized --input '$inputs' --library-name '$params.library_name'
    """
}

process load_molprops {

    container 'informaticsmatters/vs-moldb:2.0.0'
    maxForks 1
    errorStrategy { sleep(Math.pow(2, task.attempt) * 1000 as long); return 'retry' }
    maxRetries 3

    input:
    path inputs

    output:
    path 'count.txt' // the number of outputs, used for the cost

    script:
    """
    python -m moldb.load_molprops --input '$inputs'

    # count the number of lines
    COUNT=\$(wc -l < '$inputs')
    echo \$COUNT > count.txt
    """
}

process load_enum {

    container 'informaticsmatters/vs-moldb:2.0.0'
    maxForks 1
    errorStrategy { sleep(Math.pow(2, task.attempt) * 1000 as long); return 'retry' }
    maxRetries 3

    input:
    path inputs

    output:
    path 'count.txt' // the number of outputs, used for the cost

    script:
    """
    python -m moldb.load_enums --input '$inputs' --interval $params.interval\
      ${params.purge ? '--purge' : ''}

    # count the number of lines
    COUNT=\$(wc -l < '$inputs')
    echo \$COUNT > count.txt
    """
}

process load_conf {

    container 'informaticsmatters/vs-moldb:2.0.0'
    maxForks 1
    errorStrategy { sleep(Math.pow(2, task.attempt) * 1000 as long); return 'retry' }
    maxRetries 3

    input:
    path inputs

    output:
    path 'count.txt' // the number of outputs, used for the cost

    script:
    """
    python -m moldb.load_confs --input '$inputs' --interval $params.interval\
      ${params.purge ? '--purge' : ''}

    # count the number of lines
    COUNT=\$(wc -l < '$inputs')
    echo \$COUNT > count.txt
    """
}