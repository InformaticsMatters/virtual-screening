params.digits = 6
params.prefix = 'x'
params.suffix = '.txt'
params.chunk_size = 1000

process split_txt {

    container 'informaticsmatters/vs-prep:latest'

    input:
    file inputs

    output:
    file params.prefix + '*' + params.suffix

    """
    split -l $params.chunk_size -d -a $params.digits --additional-suffix $params.suffix $inputs $params.prefix
    """
}