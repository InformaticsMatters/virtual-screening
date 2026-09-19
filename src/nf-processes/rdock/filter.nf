params.top = 1

// empty string means nothing will be published. Change this to where you want your outputs.
params.publish_dir = ''
params.publish_dir_mode = 'copy'

/** Sort within the groups defined by group_by_field to find the best n results (sort_field, sort_descending, top).
*/
process sd_sort_and_top {

    container 'informaticsmatters/vs-rdock:3.0.0'
    publishDir params.publish_dir ?: '.', mode: params.publish_dir_mode, enabled: params.publish_dir as boolean

    input:
    path inputs
    val sort_field      // e.g. 'SCORE'
    val sort_descending // true or false
    val group_by_field  // e.g. '_Name'
    val outputfile      // just the filename (no path). It will be placed in publish_dir

    output:
    path "${outputfile}"

    script:
    """
    sdsort -n -s -f${sort_field} -id${group_by_field} ${sort_descending ? '-r' : ''} $inputs |\
      sdfilter -f'\$_COUNT == ${params.top}' -s${group_by_field} > ${outputfile}
    """
}

/** Sort within the groups defined by group_by_field to find the best result (sort_field, sort_descending)
and sort those best results.
*/
process sd_best_sorted {

    container 'informaticsmatters/vs-rdock:3.0.0'
    publishDir params.publish_dir ?: '.', mode: params.publish_dir_mode, enabled: params.publish_dir as boolean

    input:
    path inputs
    val sort_field      // e.g. 'SCORE'
    val sort_descending // true or false
    val group_by_field  // e.g. '_Name'
    val outputfile      // just the filename (no path). It will be placed in publish_dir

    output:
    path "${outputfile}"

    script:
    """
    sdsort -n -s -f${sort_field} -id${group_by_field} ${sort_descending ? '-r' : ''} $inputs |\
      sdfilter -f'\$_COUNT == 1' -s${group_by_field} |\
      sdsort -n -f${sort_field} ${sort_descending ? '-r' : ''} > ${outputfile}
    """
}


/** Sort within the groups defined by group_by_field to find the best result (sort_field, sort_descending)
and sort those best results and keep the top n (params.top)
*/
process sd_best_sorted_top {

    container 'informaticsmatters/vs-rdock:3.0.0'
    publishDir params.publish_dir ?: '.', mode: params.publish_dir_mode, enabled: params.publish_dir as boolean

    input:
    path inputs
    val sort_field      // e.g. 'SCORE'
    val sort_descending // true or false
    val group_by_field  // e.g. '_Name'
    val outputfile      // just the filename (no path). It will be placed in publish_dir

    output:
    path "${outputfile}"

    script:
    """
    sdsort -n -s -f${sort_field} -id${group_by_field} ${sort_descending ? '-r' : ''} $inputs |\
      sdfilter -f'\$_COUNT == 1' -s${group_by_field} |\
      sdsort -n -f${sort_field} ${sort_descending ? '-r' : ''} |\
      sdfilter -f'\$_REC <= ${params.top}' > ${outputfile}
    """
}