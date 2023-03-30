# Job: pk-tmax-cmax-sim

This describes how to run the `pk-tmax-cmax-sim` job from the `dmpk` category in the `im-dmpk` collection.

## What the job does

This job generates estimates of absorption from t1/2 and Tmax after po administration (one-compartment)
It is based on original work by Amit Kumar Garg <a.garg@sygnaturediscovery.com>.

## Implementation details

* Job implementation: [pk_tmax_cmax_sim.py](/dmpk/pk_tmax_cmax_sim.py)
* Job definition: `jobs.pk-tmax-cmax-sim` in [dmpk.yaml](/data-manager/dmpk.yaml)

## How to run the job

The job has no file inputs.

### Options

* **Output file base name**: base name for outputs (.png and .txt will be generated)
* **Half life**: half life in hours
* **Half life absorption**: absorption half like in hours
* **Initial dose**: initial dose in mg
* **AUC**: AUC value in mg/L*hr
* **Time**: Time in hours
* **Plot width**: Plot width in inches (default 10)
* **Plot height**: Plot height in inches (default 4)
* **Font size**: Font size in points (default 12)

### Outputs

Two files are generated:

1. an image file in PNG format with the plots of the absorption with the specified parameters
2. a text file listing the inputs and results

Example output would look like this:

![pk_tmax_cmax_sim output](pk_tmax_cmax_sim.png)

The corresponding text file would contain:
```
Inputs:
  half-life = 0.79
  absorption = 0.5
  dose = 0.14
  auc = 0.88
  time = 8.0
Outputs:
  Tmax(hr) = 0.899
  Cmax(mg/L) = 0.351
  Kel(hr-1) = 0.877
  Ka(hr-1) = 1.39
  V/F(L) = 0.181
```