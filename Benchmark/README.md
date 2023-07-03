# Benchmark

Scripts in this folder carries out all of the simulation experiments.

## Usage
Most of the commands to run the simulation has been given in the README.md in main folder.

Therefore, we only introduce how to tune the parameters here.

`Paper_Becnhmark_ErrorProp.jl` and `Paper_Benchmark_Sim_AIO.jl` runs with multi-processing for fastest execution. However, multi-processing will take up so much RAM (~13GB for 24 processes) and make your machine fry. To reduce the number of process to run the program (definitely slower, but safer and prevent bricking your machine), you can manually set the # of process you run the program.

To do so, you need to change the following lines in the script.
```
# line 7 in Paper_Benchmark_Sim_AIO.jl
# line 9 in Paper_Benchmark_ErrorProp.jl
addprocs(parse(Int,ENV["NUMBER_OF_PROCESSORS"])) -> addprocs(NPROCs)
```
NPROCs should not exceed your thread count.
