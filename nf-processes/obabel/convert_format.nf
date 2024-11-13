/** Convert molecule file formats using OpenBabel.
The input and output formats are determined using the file extensions
*/
params.scratch = false
params.input_extensions = ['.pdb', '.mol2']
params.output_extension = '.pdbqt'
params.output_file = 'reformatted'

process convert_format {

    container 'informaticsmatters/vs-prep:latest'
    scratch params.scratch

    input:
    file input

    output:
    file params.output_file + params.output_extension

    script:
    found = false
    params.input_extensions.each() {
      if (input.name.endsWith(it))
         found = true
    }
    if ( found )
       """
       echo 'Converting ${input.name} to ${params.output_extension} format'
       obabel $input -O '${params.output_file + params.output_extension}'
       """

    else if ( input.name.endsWith(params.output_extension) )
      """
      cp $input '${params.output_file + params.output_extension}'
      """

    else
      """
      echo 'Input ${input.name} must be in one of ${params.input_extensions.join(' ')} formats'
      exit 1
      """
}