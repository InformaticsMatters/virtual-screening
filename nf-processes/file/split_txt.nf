params.digits = 6
params.prefix = 'x'
params.suffix = '.txt'
params.header = false // if true then skip the first line
params.chunk_size = 1000

process split_txt {

    container 'informaticsmatters/vs-prep:latest'

    input:
    file inputs

    output:
    file params.prefix + '*' + params.suffix

    """
    ${params.header ? 'tail -n +2' : 'cat'} $inputs |
    split -l $params.chunk_size -d -a $params.digits --additional-suffix $params.suffix - $params.prefix
    """
}