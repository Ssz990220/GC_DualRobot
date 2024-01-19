using Pkg # Takes a while to import Plot Packages
Pkg.activate(".")
using Plots, StatsPlots

begin   # import calibrator & configs
    include("./project_utils/DualRobotCalibrator.jl")
    using .DualRobotCalibrator
    
    const N = 3;    # Number of dataset pairs for Ma
    const UNIT = "m";   # UNIT is vital
    const NOISE_DISTRIBUTION  = "normal";   # normal distributed noise as default, can be replaced with "uniform"
    
    include("./Benchmark/BenchmarkUtils/configAll.jl")
end

begin   # data gen config
    n_samples = 72;     # A reasonable-sized data set
    n_samplesNa = n_samples ÷ N;
    conf = ((max_iter = 400, μ = [1,1e-5,1,1,1],m=3, η=0.1, err = "NAN", stop = true,τ=0.3,reg=true,svd=false))
    noise_angles = 0.5/180*π;       # A reasonable noise level
    noise_dis = UNIT == "m" ? 0.5/1000 : 0.5;   # A reasonable noise level
end
begin
    begin # solve
        data = get_data(n_samples,noise_angles,noise_dis);   # random data generator

        # resultG3 = G3(data,confG3)
        resultW = wang_cic(data,conf)
        resultL = liao(data,conf)
        resultFu = fu_bayro13(data,conf)
        # resultGA = GA_Full(data,confG3)

        dataMa = get_dataMa(n_samplesNa,noise_angles,noise_dis);
        resultMa = Ma(data,confMa);
    end;

    begin   #evaluate
        # eG3 = get_error_XYZ(resultG3,Ans)
        eW = get_error_XYZ(resultW,Ans)
        eL = get_error_XYZ(resultL,Ans)
        eFu = get_error_XYZ(resultFu,Ans)
        eMa = get_error_XYZ(resultMa,Ans)
        # eGA = get_error_XYZ(resultGA,Ans)
        # e = [eG3 eW eL eFu eGA eMa]'
        e = [eW eL eFu eMa]'
    end;

    plotMa = false
    l = @layout [a b]
    if plotMa
        f1 = groupedbar(e[1:4,1:3]*1e3,xticks = (1:4,["Wang" "Wu" "Fu" "Ma"]),labels = ["Rx" "Ry" "Rz"],ylabel = "rot error / 1e-3 rad")
        f2 = groupedbar(e[1:4,4:6]*1e3,xticks = (1:4,["Wang" "Wu" "Fu" "Ma"]),labels = ["tx" "ty" "tz"],ylabel = "transl error/mm")
        plot([f1,f2],layout=l)
    else
        f1 = groupedbar(e[1:3,1:3]*1e3,xticks = (1:3,["Wang" "Wu" "Fu"]),labels = ["Rx" "Ry" "Rz"],ylabel = "rot error / 1e-3 rad")
        f2 = groupedbar(e[1:3,4:6]*1e3,xticks = (1:3,["Wang" "Wu" "Fu"]),labels = ["tx" "ty" "tz"],ylabel = "transl error/mm")
        plot(f1,f2,layout=l)
    end
end
mkpath(string("./images/"))
savefig(string("./images/demo.svg"))