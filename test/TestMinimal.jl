include("../project_utils/DualRobotCalibrator.jl")
using .DualRobotCalibrator, StaticArrays, Printf
global N = 3; global UNIT = "m"; NOISE_DISTRIBUTION = "normal"
include("../Benchmark/BenchmarkUtils/configAll.jl")

conf = ((max_iter = 200, μ = [1,1e-5,1,1,1],m=1, η=0.1, err = "NAN", stop = true,τ=0.3,reg=true,svd=false))
confG3 = ((max_iter = 200, m = 3, η=0.3,μ=0.1,scale=1, err = "NAN", stop = true, reg = true, τ=0.3))
Ntime = 100;

function main(args)
    counterWu = 0;
    counterWang = 0;
    counterG3 = 0;
    n = parse(Int, args[1]);
    for i = 1:Ntime
        data = get_data(n,0.0,0.0);
        result_Wu = get_error_XYZ(DualRobotCalibrator.liao(data,conf),Ans)
        result_Wang = get_error_XYZ(DualRobotCalibrator.wang_cic(data,conf),Ans)
        result_G3 = get_error_XYZ(DualRobotCalibrator.G3(data,confG3),Ans)
        counterWu = maximum(result_Wu) > 1e-10 ? counterWu + 1 : counterWu;
        counterWang = maximum(result_Wang) > 1e-10 ? counterWang + 1 : counterWang;
        counterG3 = maximum(result_G3) > 1e-10 ? counterG3 + 1 : counterG3;
    end
    @printf "%d data points are collected...\n" n
    @printf "Wu: %d times failed in %d trails\n" counterWu Ntime
    @printf "Wang: %d times failed in %d trails\n" counterWang Ntime
    @printf "G3: %d times failed in %d trails\n" counterG3 Ntime
end;

main(ARGS)

