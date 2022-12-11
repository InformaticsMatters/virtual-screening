# Job: filter-sdf

This describes how to run the `sort-sdf` and `filter-sdf` jobs from the `virtual screening` category in the 
`im-virtual-screening` collection.

## What the job does

The `sort-sdf` job sorts a SD-file's contents and optionally keeps the best *n* records.

The `filter-sdf` job sorts the records of a SD-file within consecutive groups keeping the best, then sorts those records
and optionally keeps the best *n* records. The records within a group must be identified by the value of a particular
field (the `groupBy` field), and must be consecutive.

## Implementation details

The work in this job is done by the `sdsort` and `sdfilter` programs from rDock. Consult the 
[rDock documentation](http://rdock.sourceforge.net/wp-content/uploads/2015/08/rDock_User_Guide.pdf) for details.

* Job definition: `jobs.filter-sdf` in [im-virtual-screening.yaml](/data-manager/im-virtual-screening.yaml)
* Scripts: [rdock_filter_sdf_sort_top.sh](rdock_filter_sdf_sort_top.sh) and 
  [rdock_filter_sdf_group.sh](rdock_filter_sdf_group.sh)

## How to run the job

### Inputs

* **Molecules to filter**: the input SD-file

### Options

* **Output file name**: the name for the output SD file
* **Sort field**: the name of the field to sort by (assumed to have numeric values)
* * **Sort descending**: Whether to sort descending or ascending (default false, so sort is ascending)
* **Group by field**: the name of the field used to identify the consecutive groups within the input. If not defined the 
title line (first line of the SDF record) is used. This is only used in the `filter-sdf` variant.
* **Keep best n records**: Keep this many top records after sorting.