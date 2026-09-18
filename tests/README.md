# Running tests

There are a few unit tests written. More need writing.
These are located under the tests directory, in a directory which has the name of the
conda environment in which they should be run. e.g. the tests in `tests/im-vs-prep`
should be run from the `im-vs-prep` conda environment which is defined in the 
`environment-im-vs-prep.yaml` file.

To run the tests you need to set the `PYTHONPATH` environment variable to the `src`
directory of this repo and then run the tests from the top level directory. e.g.

```commandline
$ export PYTHONPATH=$PWD/src
$ pytest tests/im-vs-prep/
```