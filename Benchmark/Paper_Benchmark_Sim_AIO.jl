## All in one Simulation Benchmark

using Distributed, Printf, LaTeXStrings
using JLD, HDF5

@info "Starting Distribution Pools..."
addprocs(parse(Int,ENV["NUMBER_OF_PROCESSORS"]))
# Loading Packages...
@info "Loading Packages..."
@everywhere include("../project_utils/DualRobotCalibrator.jl")     ## Solver Wrappers Defined here
@everywhere using .DualRobotCalibrator, StaticArrays

@everywhere const N = 3;

function main(args)
    @everywhere const UNIT = "m"; 
    if args[1] == "N" || args[1] == "normal" || args[1] == "Normal" || args[1] == "n"
        @everywhere const NOISE_DISTRIBUTION  = "normal";
    elseif args[1] == "U" || args[1] == "Uniform" || args[1] == "u" || args[1] == "uniform"
        @everywhere const NOISE_DISTRIBUTION  = "uniform";
    end

    @info string(NOISE_DISTRIBUTION, " noise is injected")

    n_times = parse(Int,args[2])
    @spawnat :any (:(n_times = parse(Int,args[2])));

    name = string("Paper_Bench_Sim_AIO_",n_times,"_",NOISE_DISTRIBUTION)
    @info string("Result will be saved at ./images/",name,"/")
    
    # Config
    @everywhere include("./BenchmarkUtils/configAll.jl")

    @info "Start Running Sim for Wang, Wu, Fu, Fu1/3, Proposed..."
    @info "Functions are compiled in the first batch, thus it takes longer than the rest..."
    @time Error = benchmark_p(n_angles, n_dises, n_sample, n_times, get_data, [liao, wang_cic, G3,fu_bayro,fu_bayro13, GA_Full], err_func, DIM_ERROR, [conf,conf,confG3,conf,conf, confₘ₂]);

    ## Saving Result
    @info "Saving result..."
    begin
        mkpath(string("./images/",name))
        save(string("./images/",name,"/Result.jld"),"Error", Error.s,"N_sample",n_sample,"unit",UNIT)
    end

    ## Run Ma's Sim
    @info "Start Running Sim for Ma..."
    @info "Functions are compiled in the first batch, thus it takes longer than the rest..."
    ErrorMa = benchmark_p(n_angles, n_dises, n_sampleMa, n_times, get_dataMa, [Ma], err_func, DIM_ERROR,[confMa]);

    
    ## Run MaBatch's Sim
    # @info "Start Running Sim for MaBatch..."
    # @info "Functions are compiled in the first batch, thus it takes longer than the rest..."
    # ErrorMaOne = benchmark_p(n_angles, n_dises, n_sample, n_times, get_dataMaOne, [MaOne], err_func, DIM_ERROR,[confMa]);

    ## Save Ma's Result
    @info "Saving result..."
    save(string("./images/",name,"/ResultMa.jld"),"Error", ErrorMa.s,"N_sample",n_sampleMa*N,"unit",UNIT)
    # save(string("./images/",name,"/ResultMaOne.jld"),"Error", ErrorMaOne.s,"N_sample",n_sample,"unit",UNIT)
    
    @info string("Result saved at ","./images/",name,"/")
end

main(ARGS)
