This describes how to run the *filter-sdf* job from the *virtual screening* category in the *im-virtual-screening* collection.

## What the job does

This job filters an SD-file's contents. It has two primary options, a field name to use for sorting the records and a 
field name to group consecutive records. It is assumed that your input SD file has groups of consecutive records 
identified by the `groupBy` field. Within each group the records are sorted according to the value of the `sortField`,
then the best one is kept and finally all those best records are again sorted so that the first record in the output SD 
file is the best record.

## Implementation details

The work in this job is done by the `sdsort` and `sdfilter` programs from rDock. Consult the 
[rDock documentation](http://rdock.sourceforge.net/wp-content/uploads/2015/08/rDock_User_Guide.pdf) for details.

Job definition: `jobs.filter-sdf` in [/data-manager/virtual-screening.yaml]()

## How to run the job

### Inputs

**Molecules to filter**: the input SD file

### Options

**Output file name**: the name for the output SD file
**Sort field**: the name of the field to sort by (assumed to have numeric values)
**Group by field**: the name of the field used to identify the consecutive groups within the input. If not defined the 
title line (first line of the SDF record) is used.
**Sort descending**: Whether to sort descending or ascending (default false, so sort is ascending)