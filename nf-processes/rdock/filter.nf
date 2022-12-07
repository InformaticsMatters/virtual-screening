params.sort_field = 'SCORE'
params.sort_descending = false
params.group_by_field = "_Name"
params.top = 1

// outputfile should be just the filename (no path). It will be placed in publish_dir
params.outputfile = 'results.sdf'
// empty string means nothing will be published. Change this to where you want your outputs.
params.publish_dir = ''
params.publish_dir_mode = 'copy'

/** Sort within the groups defined by group_by_field to find the best n results (sort_field, sort_descending, top).
*/
process sd_sort_and_top {

    container 'informaticsmatters/vs-rdock:latest'
    if (params.publish_dir) { publishDir params.publish_dir, mode: params.publish_dir_mode }

    input:
    path inputs

    output:
    path params.outputfile

    """
    sdsort -n -s -f${params.sort_field} -id${params.group_by_field} ${params.sort_descending ? '-r' : ''} $inputs |\
      sdfilter -f'\$_COUNT == ${params.top}' -s${params.group_by_field} > ${params.outputfile}
    """
}

/** Sort within the groups defined by group_by_field to find the best result (sort_field, sort_descending)
and sort those best results.
*/
process sd_best_sorted {

    container 'informaticsmatters/vs-rdock:latest'
    if (params.publish_dir) { publishDir params.publish_dir, mode: params.publish_dir_mode }

    input:
    path inputs

    output:
    path params.outputfile

    """
    sdsort -n -s -f${params.sort_field} -id${params.group_by_field} ${params.sort_descending ? '-r' : ''} $inputs |\
      sdfilter -f'\$_COUNT == 1' -s${params.group_by_field} |\
      sdsort -n -f${params.sort_field} ${params.sort_descending ? '-r' : ''} > ${params.outputfile}
    """
}


/** Sort within the groups defined by group_by_field to find the best result (sort_field, sort_descending)
and sort those best results and keep the top n (params.top)
*/
process sd_best_sorted_top {

    container 'informaticsmatters/vs-rdock:latest'
    if (params.publish_dir) { publishDir params.publish_dir, mode: params.publish_dir_mode }

    input:
    path inputs

    output:
    path params.outputfile

    """
    sdsort -n -s -f${params.sort_field} -id${params.group_by_field} ${params.sort_descending ? '-r' : ''} $inputs |\
      sdfilter -f'\$_COUNT == 1' -s${params.group_by_field} |\
      sdsort -n -f${params.sort_field} ${params.sort_descending ? '-r' : ''} |\
      sdfilter -f'\$_REC <= ${params.top}' > ${params.outputfile}
    """
}