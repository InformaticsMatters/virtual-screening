# A brief guide to developing tools and jobs

This guide describes the process for creating new tools, and how to allow them to be run in Squonk2 as *jobs*.

These tools are command line tools that typically read some input and write some output. That same tool can be
packaged up in a Docker container image and, given a *job definition*, can be executed in Squonk2 as a job.

## 1. Creating a repository

We try to create a set of jobs that fall into a particular category. The ones in this GitHub repository are related
to virtual screening. If the tool you want to create is related to virtual screening then it might be appropriate to 
add it to this repository (fork the repository and give us a pull request with your changes).

However, if your tool is not related to virtual screening then you should create a new repository. Squonk2 can load jobs 
from multiple repositories, so it helps to keep things separate where appropriate.

If you do create your own repository, then we suggest to follow the patterns used here, though if you have good reason
to do things in other ways then that's probably OK too.

## 2. Creating the conda environment

We assume here that you are creating your tool with Python.
You don't have to use Python. See [here](https://github.com/InformaticsMatters/squonk2-cdk) for a Squonk2
job that is written in Java.

The easiest way to run your tools is to use an existing conda environment, or if none is appropriate, then create a
new one. See the various `environment-*.yaml` files for the ones that already exist. For example the
[environment-im-prep.yaml]() file defines a conda environment that contains RDKit, OpenBabel and some other tools.

Create your conda environment with something like this:
```
conda env create -f environment-im-prep.yaml
```
Then activate your environment with something like this:
```
conda activate im-vs-prep
```

## 3. Creating the Python module

If you are using Python then we strongly suggest you use Python3 if possible.

We suggest you follow the patterns used in the modules in this repository:

* Use [argparse](https://docs.python.org/3/library/argparse.html) for handling the command line options.
* Follow the conventions used in these modules for the naming of the command line options.
* Provide a Python function in your module that allows the function of the tool be used from another module. e.g.
  have the main entrypoint parse the commandline arguments, prepare anything that is necessary and then call a function
  that does the real work. That same function should be callable from another Python module or script.
* Keep the module relatively simple and concise. If your tool needs to do one thing followed by another thing then
  probably you should create two separate tools.
* Consider using the utility functions provided e.g. in [rdkit_utils.py]().
* Log information to STDOUT. See next section for more.

## 4. Logging

If you are wanting your tool to be a Squonk2 job then you should pay attention to logging *Data Manager* event and cost 
messages.

Event messages are lines in STDOUT conforming to a particular pattern that the Data Manager job executor looks for and
reports as *events* that are shown in the job execution UI. For instance an event message looks like this:
```
2022-02-03T16:39:27+00:00 # INFO -EVENT- Hello World!
```
The bit after the `-EVENT-` section is the message that is reported.
Not everything that is logged needs to be a DataManger Event. You might want to log some more verbose log messages that
won't get reported through the Squonk2 Data Manager. Typically you only want a small number (5-10) of events to be 
reported.

Cost messages are lines in STDOUT conforming to a particular pattern that the Data Manager job executor looks for and
uses to monitor job usage. Ultimately this could result in the user incurring charges and you get paid for them using
your tool, but you should write cost messages even if you are not wanting to monitise your work (a zero cost can be
assigned to your job) as it allows usage of your job to be recorded. A cost message looks like this:
```
2022-02-03T16:40:16+00:00 # INFO -COST- 5.7 1
```
The two numbers after the `-COST-` bit are:

1. The cost in some arbitrary units e.g. number of molecules processed, or whatever is an appropriate measure of the 'cost'
2. The sequence of the cost event, this message being the first (1).

With cost messages it is usually best to log the cost event once processing is finished. In this case just one message 
needs to be logged. However, if your tool is likely to run for a long time it is better to log cost messages at regular 
intervals e.g. every 5 minutes, or every 100,000 molecules processed.

The cost value (the first number) can either be incremental or absolute. Consider these two sets of messages:

```
2022-02-03T16:40:16+00:00 # INFO -COST- 5.7 1
2022-02-03T16:40:16+00:00 # INFO -COST- 8.3 2
```

```
2022-02-03T16:40:16+00:00 # INFO -COST- +5.7 1
2022-02-03T16:40:16+00:00 # INFO -COST- +8.3 2
```

The first pair is absolute, the second incremental. The difference is in the second the cost is prefixed with a `+`.
In the first pair the final cost value is 8.3, in the second it is 14.0. Using absolute costs is probably easier and
better, but there are times when you might want to use incremental costs.

If this all sounds a bit complex and you're thinking of ignoring this then don't! We have created a 
[simple library](https://github.com/InformaticsMatters/data-manager-job-utilities)
that makes logging these messages very simple. It's on [PyPi](https://pypi.org/project/im-data-manager-job-utilities/)
and easy to install and use. The modules in this repository contain lots of examples for how to use it. For instance,
the [minimize.py]() module emits event messages like this:
```
DmLog.emit_event("Force field could not be set up for molecule", count)
```
and emits cost messages at regular intervals like this:
```
if success % 10000 == 0:
    DmLog.emit_cost(success)
```

## 5. Creating the Dockerfile

To be runnable as a Squonk2 job then you must build a Docker container image that allows your tool to be run using 
the command that you will define in the job definition (see the next section).

If you are using one of our conda environments then you can use the corresponding container image, but if there is 
no container image suitable for your tool then create a Dockerfile that allows one to be created. Look at the
`Dockerfile-*` files in this repo as examples. You can even use a container image that contains conda as your base
image and conda install the same packages as your conda environment (see [Dockerfile-prep]() as an example).

Make sure you push the Docker image to a container repository if it is going to be used in Squonk2.
Generally you should use a public container registry such as [Dockerhub](https://hub.docker.com/) but you can also use
a private prository if necessary. If so contact us about how to specify the appropriate pull secrets in Squonk2.

## 6. Creating the job definition

To get your tool to run in Squonk2 as a *job* you need to write a *job definition* (a YAML file) and add that file to a
*job manifest* (either an existing one, or create a new one). Both of these files must exist in the `data-manager`
directory at the top level of your repository.

The job manifest is simple. e.g. [data-manager/manifest-virtual-screening.yaml]() looks like this:
```yaml
---
kind: DataManagerManifest
kind-version: '2021.1'

job-definition-files:
- virtual-screening.yaml
- rdkit.yaml
- xchem.yaml
```

It simply defines 3 job defintion YAML files.

The *job definition* YAML file is more complex. It's best understood by looking at an example, for instance
[data-manager/rdkit.yaml]() which defines a number of jobs that use RDKit.

The file can define one or more jobs. The file and job have these key sections:

```yaml
collection: rdkit
```
This is a top level property and all jobs in this file belong to this collection.

All following sections are job specific.

```yaml
    category: comp chem
```
This defines the category of the job.

The `collection` and the `category` can be used for filtering jobs in the Squonk2 UI. 

```yaml
    image:
      name: informaticsmatters/vs-prep
      tag: 'latest'
      project-directory: /data
      working-directory: /data
```
This defines the container image you built and pushed in the previous section.

```yaml
    command: >-
      /code/max_min_picker.py --input '{{ inputFile }}'
      {% if seeds is defined %}--seeds{% for file in seeds %} '{{ file }}'{% endfor %}{% endif %}
      --output '{{ outputFile }}'
      --count {{ count }}
      {% if threshold is defined %}--threshold {{ threshold }}{% endif %}
      --interval 10000
```
This defines the command that is executed. It uses Jinja2 templating to fill in the values of the inputs and options
that we'll see next. The filled in template is used as the command that is executed when the job is run in the Squonk2
Data Manager in the Kubernetes cluster. Think of it being the `<command>` bit when running with docker:
```commandline
docker run -it yourorg/yourcontainer <command>
```
Be careful to prevent hacking attacks by putting substituted strings in single quotes.

```yaml
      inputs:
        type: object
        required:
        - inputFile
        properties:
          inputFile:
            title: Molecules to pick from
            mime-types:
            - squonk/x-smiles
            type: file
          seeds:
            title: Molecules that are already picked
            mime-types:
            - squonk/x-smiles
            type: file
            multiple: true
```
This defines the files that are the inputs to your tool. In this case there are two inputs, with `seeds` being optional.

```yaml
      outputs:
        type: object
        properties:
          outputFile:
            title: Output file
            mime-types:
            - chemical/x-csv
            creates: '{{ outputFile }}'
            type: file
```
This defines the outputs of your tool - the files that are created.

```yaml
      options:
        type: object
        required:
        - count
        properties:
          outputFile:
            title: Output file name
            type: string
            pattern: "^[A-Za-z0-9_/\\.\\-]+$"
            default: diverse.smi
          count:
            title: Number of molecules to pick
            type: integer
            minimum: 1
          threshold:
            title: Similarity threshold
            type: number
            minimum: 0
            maximum: 1
```
This defines the user definable options for your job. This uses [JSON schema](https://json-schema.org/) notation with 
the user interface for the job executor in Squonk2 being automatically generated from this. Try to include as much
validation as possible (especially with string options) to prevent hacking attempts.  

This is not an exhaustable list of the sections, but covers the key aspects.
Look at the existing job definitions for full details, or contact us if you need more info.

## 7. Creating the job tester

Writing the job definition is tricky and can be subject to silly typo or formatting errors.
To assist with this we have created [jote](https://github.com/InformaticsMatters/data-manager-job-tester),
the job tester.

Jote lets you test a job definition. It can:
- validate the YAML
- perform some basic checks
- execute the job in Docker in a way that is very similar to how it will execute in Kubernetes so that you can test
  that it runs correctly

As an example, here is a test definition. It lives inside the job definition in the job definition YAML file. It comes from the same 
RDKit job example used above which includes additional comments explaining the meaning of the elements:
```yaml
    tests:
      simple-execution:
        inputs:
          inputFile: data/mols.smi
        options:
          outputFile: diverse.smi
          count: 100
        checks:
          exitCode: 0
          outputs:
          - name: diverse.smi
            checks:
            - exists: true
            - lineCount: 100
```
You will notice that each test defines some inputs and options that are needed by the job and defines some basic checks
on the outputs that should be created by running the job. The inputs and options are used to generate the *command*
(see above) and then that command is run in docker to generate the outputs which are then checked for validity.

There is a [environment-im-jote.yaml]() conda environment file that you can use to create a jote conda environment.
Jote is run against a particular manifest file, and you can restrict it to jobs from a particular collection or even a 
particular job (use `jote --help` to see all the options).

For instance, to run the test we have been looking at:

```commandline
(im-jote) $ jote -m manifest-virtual-screening.yaml -c rdkit -j max-min-picker
# Using manifest "manifest-virtual-screening.yaml"
# Found 10 tests
# Limiting to Collection "rdkit"
# Limiting to Job "max-min-picker"
  ---
+ collection=rdkit job=max-min-picker test=simple-execution
> run-level=Undefined
> image=informaticsmatters/vs-prep:latest
> command="/code/max_min_picker.py --input 'mols.smi'  --output 'diverse.smi' --count 100  --interval 10000"
# Creating test environment...
# docker-compose (1.29.2, build unknown)
# Created
# path=/data/github/im/virtual-screening/data-manager/jote/rdkit.max-min-picker.simple-execution
input_files=['data/mols.smi']
# Copying inputs (from "${PWD}/data")...
# + data/mols.smi
# Copied
# Executing the test ("docker-compose up")...
# Executed (exit code 0)
# Checking...
# - diverse.smi
#   exists (True) [OK]
#   lineCount (100) [OK]
# Checked
# Deleting the test...
# Deleted
  ---
Done (OK) passed=1 skipped=0 ignored=0 failed=0 
```

## 8. Deploying to Squonk2

Jobs definitions are loaded into Squonk2 Data Manager using the manifest file mentioned above. This allows to provide
a level of granularity of types of jobs. An individual manifest can be loaded into one Squonk2 instance but not into
a different instance.

You will need to contact the administrator of your Squonk2 instance who will use the
`/admin/job-manifest` API endpoint to load your manifest. 
Also, if your manifest or the job definitions it includes
changes the administrator will need to use the `/admin/job-manifest/load` endpoint to reload the job manifests.