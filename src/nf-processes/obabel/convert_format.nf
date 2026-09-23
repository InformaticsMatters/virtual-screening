/** Convert molecule file formats using OpenBabel.
The input and output formats are determined using the file extensions
*/
params.scratch = false

process convert_format {

    container 'informaticsmatters/vs-prep:3.1.0'
    scratch params.scratch

    input:
    path input
    val input_extensions // e.g. ['.pdb', '.mol2']
    val output_extension // e.g. '.pdbqt'
    val output_file      // the output file name without the extension

    output:
    path "${output_file}${output_extension}"

    script:
    def found = input_extensions.any { ext -> input.name.endsWith(ext) }
    if (found)
        """
        echo 'Converting ${input.name} to ${output_extension} format'
        obabel $input -O '${output_file}${output_extension}'
        """

    else if (input.name.endsWith(output_extension))
        """
        cp $input '${output_file}${output_extension}'
        """

    else
        """
        echo 'Input ${input.name} must be in one of ${input_extensions.join(' ')} formats'
        exit 1
        """
}
