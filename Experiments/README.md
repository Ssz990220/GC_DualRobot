# Experiment
## Description
Two physical experiments are carried out to verify the effectiveness of the proposed algorithm. Results are compared with state-of-the-art algorithms.


https://github.com/Ssz990220/GC_DualRobot/assets/39759334/e73ba843-6a4b-46a1-a7c1-af8639ed3166


## Experiment 1:
*Section V: Cross Validation* in the paper

![Image](../Assets/CrossRef.png "Cross Reference Result")
## Experiment 2:
*Section V: Measurement Experiments* in the paper

![Image](../Assets/Measure%20Region.jpg "Measurement Region")

## Dataset Setup
The dataset is available at [here](https://drive.google.com/file/d/1KBOPB5leS9vUCu4oeyLxHO-R50ZSIwa4/view?usp=sharing). The dataset is a zip file, which contains the following folders:

* result/       # Comes from data.zip
* traj/         # Comes from data.zip
* meshes/       # Comes from data.zip

These folders should be extracted at the same level as `AXBYCZ.jl` and `pre_process.jl`. The structure of the folder should be
```
Experiments
    result/       # Comes from data.zip
    traj/         # Comes from data.zip
    meshes/       # Comes from data.zip
    AXBYCZ.jl
    pre_process.jl
    validation_set_processing.m
```

## Reproduction
The `AXBYCZ.jl` is the main script to run the experiment. The `pre_process.jl` is a helper functions to process the data. The `validation_set_processing.m` is a script to process the data for the measurement experiment.

To reproduce the cross-validation result, you need to
1. Open a Pluto server
2. Open `AXBYCZ.jl` in Pluto
3. activate the deactivated cell
![Image](../Assets/ActivateCell.png "Activate Cell")

This may take a while to run. The result will be plotted in the following cell.

To reproduce the measurement experiment result, you need to
1. Wait for `AXBYCZ.jl` to finish running
2. Run `validation_set_processing.m` in MATLAB
