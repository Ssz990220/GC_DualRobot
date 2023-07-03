NAME = "Paper_Sim_ErrorProp"
labels = ["Wu" "Wang" "G3" "Fu" "Fu13" "Ma"]
discribe = "Benchmark on different Algorithms under different noise level"

using Distributed, Printf, LaTeXStrings
using JLD, HDF5

@info "Starting Distribution Pools..."
addprocs(parse(Int,ENV["NUMBER_OF_PROCESSORS"]))
# Loading Packages...
@info "Loading Packages..."
@everywhere include("../project_utils/DualRobotCalibrator.jl")     ## Solver Wrappers Defined here
@everywhere using .DualRobotCalibrator

@everywhere const N = 3;        # Number of dataset pairs for Ma


function main()
    @everywhere const UNIT = "m"; 
    @everywhere const NOISE_DISTRIBUTION  = "normal";
    @everywhere const n_times = 1000;

    name = string(NAME,"_",n_times,"_",NOISE_DISTRIBUTION)
    @info string("Result will be saved at ./images/",name,"/")
    
    # Config All
    @everywhere include("./BenchmarkUtils/configAll.jl")
    @everywhere n_sample = [36]
    @everywhere n_sampleMa = [36÷N]
    @everywhere conf = ((max_iter = 400, μ = [1,1e-5,1,1,1],m=3, η=0.1, err = "NAN", stop = true,τ=0.3,reg=true,svd=false))

    ##########################################
    ##          Rot Error Prop Exp          ##
    ##########################################
    @info "Start Testing Rotational Error Propagation..."
    @info "Only Rotational noise is injected..."

    # EXP Config    
    @everywhere n_angles = (0.0:0.3:3.0)./180.0.*π  # in radius
    @everywhere n_dises = [0.0 for i = 1:size(n_angles,1)]

    @info "Start Running Sim for Wang, Wu, Fu, Fu1/3, Proposed..."
    @info "Functions are compiled in the first batch, thus it takes longer than the rest..."
    @time ErrorR = benchmark_p(n_angles, n_dises, n_sample, n_times, get_data, [liao, wang_cic, G3,fu_bayro,fu_bayro13], err_func, DIM_ERROR, [conf,conf,confG3,conf,conf], disp="noise");

    ## Run Ma's Sim
    @info "Start Running Sim for Ma..."
    @info "Functions are compiled in the first batch, thus it takes longer than the rest..."
    @time ErrorMaR = benchmark_p(n_angles, n_dises, n_sampleMa, n_times, get_dataMa, [Ma], err_func, DIM_ERROR,[confMa], disp="noise");

    ## Saving Result
    Error_Rotation = cat(ErrorR.s,ErrorMaR.s;dims=5);
    @info "Saving result For Rotation Error Propagation..."
    mkpath(string("./images/",name))
    save(string("./images/",name,"/ResultR.jld"),"Error", Error_Rotation,"N_sample",n_sample,"unit",UNIT)

    #############################################
    ##          Transl Error Prop Exp          ##
    #############################################
    @info "Start Testing Translation Error Propagation..."
    @info "Only tranlsational noise is injected..."

    @everywhere n_dises2 = (0.0:0.3:3.0)./1000
    @everywhere n_angles2 = [0.0 for i = 1:size(n_dises2,1)]

    @info "Start Running Sim for Wang, Wu, Fu, Fu1/3, Proposed..."
    @info "Functions are compiled in the first batch, thus it takes longer than the rest..."
    @time ErrorT= benchmark_p(n_angles2, n_dises2, n_sample, n_times, get_data, [liao, wang_cic, G3,fu_bayro,fu_bayro13], err_func, DIM_ERROR, [conf,conf,confG3,conf,conf], disp="noise");

    ## Run Ma's Sim
    @info "Start Running Sim for Ma..."
    @info "Functions are compiled in the first batch, thus it takes longer than the rest..."
    @time ErrorMaT = benchmark_p(n_angles2, n_dises2, n_sampleMa, n_times, get_dataMa, [Ma], err_func, DIM_ERROR,[confMa], disp="noise");

    Error_Transl = cat(ErrorT.s,ErrorMaT.s;dims=5);
    ## Saving Result
    @info "Saving result For Translation Error Propagation..."
    save(string("./images/",name,"/ResultT.jld"),"Error", Error_Transl,"N_sample",n_sample,"unit",UNIT)
    
end

main()

begin
    using Statistics, Plots
    folder = string(NAME,"_",n_times,"_",NOISE_DISTRIBUTION);
    Error_Rotation = load(string("./images/",folder,"/ResultR.jld"))["Error"]
    ER = zeros(size(n_angles,1),DIM_ERROR,size(labels,2))
    i = 1
    @inbounds for j ∈ eachindex(n_angles)
        @inbounds for s ∈ eachindex(labels)
            @inbounds for m =1:DIM_ERROR
                ER[j,m,s] = Statistics.mean(abs.(Error_Rotation[j,i,m,:,s]))
                # if mod(m,2)==0
                    ER[j,m,s] = ER[j,m,s]*1000            # Convert to millimeter for better readability 
                # end
            end
        end
    end
end

begin
    Title = [string("R_",s,mode," Err") for s in ["X","Y","Z","AXBYCZ"] for mode in [" Angle"," Dis"]]
    ylabels = [string(s,mode[1]," Err / ",mode[2]) for s in ["X","Y","Z","AXB=YCZ"] for mode in [(" Rot"," 10^{-3} Rad"),(" Dis","mm")]]
    s = 1;
    n_angles_degree = n_angles ./ pi * 180
    
    l = @layout [a b; c d; e f]
    Ps = Vector{Plots.Plot{Plots.GRBackend}}(undef,6)
    for j = 1:6     # 6 errors angle-x, transl-x, angle-y, transl-y, angle-z, transl-z
        Ps[j] = plot(n_angles_degree,ER[s:end,j,2],linestyle=:dash,linecolor=RGB(109/255, 174/255, 209/255),label=false,
            linewidth=3,markershape=:utriangle,markercolor=RGB(109/255, 174/255, 209/255),markersize=5,#xlabel="Rotational Error / °",
            guidefont = (14,"times"), tickfont = (12), legendfontsize = 12, grid=true,
            ylabel = ylabels[(i-1)*DIM_ERROR+j],size=(400,300))
        plot!(n_angles_degree,ER[s:end,j,1],linestyle=:dashdotdot,linecolor=RGB(17/255, 70/255, 127/255),label=false,
            linewidth=3,markershape=:star5,markercolor=RGB(17/255, 70/255, 127/255),markersize=5)
        plot!(n_angles_degree,ER[s:end,j,4],linestyle=:dash,linecolor=RGB(245/255, 201/255, 91/255),label=false,
            linewidth=3,markershape=:hexagon,markercolor=RGB(245/255, 201/255, 91/255),markersize=5)
        plot!(n_angles_degree,ER[s:end,j,5],linestyle=:solid, linecolor=RGB(245/255, 201/255, 91/255),label=false,
            linewidth=3,markershape=:rect,markercolor=RGB(245/255, 201/255, 91/255),markersize=5)
        plot!(n_angles_degree,ER[s:end,j,3],linestyle=:solid,linecolor=RGB(182/255, 35/255, 48/255),label=false,
            linewidth=3,markershape=:circle,markercolor=RGB(182/255, 35/255, 48/255),markersize=5)
        # plot!(n_angles_degree,ER[s:end,j,6],linestyle=:solid,linecolor=RGB(248/255, 178/255, 147/255),label = false,
        #     linewidth=3,markershape=:dtriangle,markercolor=RGB(248/255, 178/255, 147/255),markersize=5)
    end
    plot(Ps...,layout=l,size=(800,800))
    savefig(string("./images/",folder,"/","RotErrProp",".svg"))
end 

begin
    Error_Translation = load(string("./images/",folder,"/ResultT.jld"))["Error"]
    Et = zeros(size(n_dises2,1),DIM_ERROR,size(labels,2))
    i = 1
    @inbounds for j ∈ eachindex(n_dises2)
        @inbounds for s ∈ eachindex(labels)
            @inbounds for m =1:DIM_ERROR
                Et[j,m,s] = Statistics.mean(abs.(Error_Translation[j,i,m,:,s]))
                # if mod(m,2)==0
                    Et[j,m,s] = Et[j,m,s]*1000            # Convert to millimeter for better readability 
                # end
            end
        end
    end
end

begin
    using Plots
    Title2 = [string("T_",s,mode," Err") for s in ["X","Y","Z","AXBYCZ"] for mode in [" Angle"," Dis"]]
    ylabels = [string(s,mode[1]," Err / ",mode[2]) for s in ["X","Y","Z","AXB=YCZ"] for mode in [(" Rot"," 10^{-3} Rad"),(" Dis","mm")]]
    s = 1;
    n_dises2_m = n_dises2 .* 1e3
    
    l = @layout [a b; c d; e f]
    Pst = Vector{Plots.Plot{Plots.GRBackend}}(undef,6)
    for j = 1:6     # 6 errors angle-x, transl-x, angle-y, transl-y, angle-z, transl-z
        Pst[j] = plot(n_dises2_m,Et[s:end,j,2],linestyle=:dash,linecolor=RGB(109/255, 174/255, 209/255),label=false,
            linewidth=3,markershape=:utriangle,markercolor=RGB(109/255, 174/255, 209/255),markersize=5,#xlabel="Translational Error / mm",
            guidefont = (14,"times"), tickfont = (12), legendfontsize = 12, grid=true,
            ylabel = ylabels[(i-1)*DIM_ERROR+j],size=(400,300))
        plot!(n_dises2_m,Et[s:end,j,1],linestyle=:dashdotdot,linecolor=RGB(17/255, 70/255, 127/255),label=false,
            linewidth=3,markershape=:star5,markercolor=RGB(17/255, 70/255, 127/255),markersize=5)
        plot!(n_dises2_m,Et[s:end,j,4],linestyle=:dash,linecolor=RGB(245/255, 201/255, 91/255),label=false,
            linewidth=3,markershape=:hexagon,markercolor=RGB(245/255, 201/255, 91/255),markersize=5)
        plot!(n_dises2_m,Et[s:end,j,5],linestyle=:solid, linecolor=RGB(245/255, 201/255, 91/255),label=false,
            linewidth=3,markershape=:rect,markercolor=RGB(245/255, 201/255, 91/255),markersize=5)
        plot!(n_dises2_m,Et[s:end,j,3],linestyle=:solid,linecolor=RGB(182/255, 35/255, 48/255),label=false,
            linewidth=3,markershape=:circle,markercolor=RGB(182/255, 35/255, 48/255),markersize=5)
        # plot!(n_angles_degree,ER[s:end,j,6],linestyle=:solid,linecolor=RGB(248/255, 178/255, 147/255),label = false,
        #     linewidth=3,markershape=:dtriangle,markercolor=RGB(248/255, 178/255, 147/255),markersize=5)
    end
    plot(Pst...,layout=l,size=(800,800))
    savefig(string("./images/",folder,"/","TranslErrProp",".svg"))
end 