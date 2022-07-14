# Running the Jaqpot models

Build the docker image
```
docker-compose build jaqpot
```

Enter the container:
```
docker run -it --rm -v $PWD:$PWD -w $PWD -u 1000:1000 informaticsmatters/vs-jaqpot:latest bash
```

Set the environment variable for the API
```
export JAQPOT_API_KEY=...
```

Check it runs:
```
./jaqpot_model_exec.py --help
```

Run a model:
```
./jaqpot_model_exec.py -i data/100.smi -m K1J1JzuZJxjY1q8YCD3h
```

