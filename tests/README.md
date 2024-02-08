# Running tests

There are a few unit tests written. More need writing.
These are located under the tests directory, in a directory which has the name of the
conda environment in which they should be run. e.g. the tests in `tests/im-vs-prep`
should be run from the `im-vs-prep` conda environment which is defined in the 
`environment-im-vs-prep.yaml` file.

To run the tests you need to set the `PYTHONPATH` environment variable to the top
level directory of this repo and then run the tests from that dir. e.g.

```commandline
$ export PYTHONPATH=$PWD
$ pytest tests/im-vs-prep/
```