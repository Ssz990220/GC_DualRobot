NAME = "Paper_Bench_Sim_AIO"
# name = "Paper_Bench_Sim_AIO_1000_normal"
labels = ["Wu" "Wang" "Proposed" "Fu" "Fu2/3" "Ma"]
levels = ["low" "mid" "high"]
begin
    using Statistics
    using JLD, HDF5
    using Plots
    using LaTeXStrings
    using Printf
end
function main(args)
    if args[1] == "N" || args[1] == "normal" || args[1] == "Normal" || args[1] == "n"
        NOISE_DISTRIBUTION  = "normal";
    elseif args[1] == "U" || args[1] == "Uniform" || args[1] == "u" || args[1] == "uniform"
        NOISE_DISTRIBUTION  = "uniform";
    end
    N = args[2]
    name = NAME * "_" * N * "_" * NOISE_DISTRIBUTION;

begin
    @info string("Results loaded from ./images/",name)
    Error = load(string("./images/",name,"/Result.jld"))["Error"]
    n_sample = load(string("./images/",name,"/Result.jld"))["N_sample"]
    unit = load(string("./images/",name,"/Result.jld"))["unit"]
    @show size(Error)
    num_angles = size(Error,1); num_sample = size(Error,2);
    DIM_ERROR = size(Error,3); n_times = size(Error,4);
    n_methods = size(Error,5);
    
    E = zeros(num_angles,num_sample,DIM_ERROR,n_methods+1)        # Add 1 dim for Ma
    @inbounds for i = 1:num_sample
        @inbounds for j = 1:num_angles
            @inbounds for s = 1:n_methods
                @inbounds for m =1:DIM_ERROR
                    valid = (.!any(abs.(Error[j,i,[1,3,5],:,s]).>0.2,dims=1)).&& (.!any(abs.(Error[j,i,[2,4,6],:,s]).>0.2,dims=1))      # Remove outliers, all for PGA method
                    Err = Error[j,i,m,valid[:],s]                                                                                       # Other methods works perfectly without outliers
                    E[j,i,m,s] = Statistics.mean(Err)                                                                                   # Replace Err with Error[j,i,m,:,s] to include outliers
                    E[j,i,m,s] = E[j,i,m,s]*1000            # Convert to millimeter & 1e-3 rad for better readability 
                end
            end
        end
    end
end

begin # Insert Ma data
    ErrorMa = load(string("./images/",name,"/ResultMa.jld"))["Error"]
    n_sampleMa = load(string("./images/",name,"/ResultMa.jld"))["N_sample"]
    unitMa = load(string("./images/",name,"/ResultMa.jld"))["unit"]
    num_anglesMa = size(ErrorMa,1); num_sampleMa = size(ErrorMa,2);
    DIM_ERROR_MA = size(ErrorMa,3); n_timesMa = size(ErrorMa,4);
    
    @assert unitMa == unit
    @assert num_anglesMa == num_angles;
    @assert num_sample == num_sampleMa
    @assert DIM_ERROR == DIM_ERROR_MA
    @assert n_times == n_timesMa

    # valid = fill(false,num_angles,num_sample,n_times)
    # @inbounds for i ∈ eachindex(n_sample)
    #     @inbounds for j ∈ eachindex(num_angles)
    #         @inbounds for s = 2
	# 			valid[j,i,:] = ((.!any(ErrorMa[j,i,[1,3,5],:,s].>0.1,dims=1)) .&& (.!any(ErrorMa[j,i,[2,4,6],:,s].>0.1,dims=1)))
	# 		end
	# 	end
	# end
    local MaFailCount = 0;
    @inbounds for i = 1:num_sample
        @inbounds for j = 1:num_angles
            @inbounds for s = 1
                @inbounds for m =1:DIM_ERROR
                    valid = (.!any(abs.(ErrorMa[j,i,[1,3,5],:,s]).>0.2,dims=1)).&& (.!any(abs.(ErrorMa[j,i,[2,4,6],:,s]).>0.2,dims=1))
                    MaFailCount = MaFailCount + n_timesMa - sum(valid)
                    Err = ErrorMa[j,i,m,valid[:],s]
                    # Err = ErrorMa[j,i,m,:,1]
                    E[j,i,m,n_methods+1] = Statistics.mean(Err)
                    E[j,i,m,n_methods+1] = E[j,i,m,n_methods+1]*1000            # Convert to millimeter & 1e-3 rad for better readability 
                end
            end
        end
    end
    @printf "%d out of %d dataset on Ma fails\n" MaFailCount/3 num_angles*n_times*size(n_sampleMa,1)
end

# begin # Insert Ma data
#     ErrorMaOne = load(string("./images/",name,"/ResultMaOne.jld"))["Error"]
#     n_sampleMa2 = load(string("./images/",name,"/ResultMaOne.jld"))["N_sample"]
#     unitMa2 = load(string("./images/",name,"/ResultMaOne.jld"))["unit"]
#     num_anglesMa2 = size(ErrorMaOne,1); num_sampleMa2 = size(ErrorMaOne,2);
#     DIM_ERROR_MA2 = size(ErrorMaOne,3); n_timesMa2 = size(ErrorMaOne,4);
    
#     @assert unitMa2 == unit
#     @assert num_anglesMa2 == num_angles;
#     @assert num_sample == num_sampleMa2
#     @assert DIM_ERROR == DIM_ERROR_MA2
#     @assert n_times == n_timesMa2

#     # valid = fill(false,num_angles,num_sample,n_times)
#     # @inbounds for i ∈ eachindex(n_sample)
#     #     @inbounds for j ∈ eachindex(num_angles)
#     #         @inbounds for s = 2
# 	# 			valid[j,i,:] = ((.!any(ErrorMa[j,i,[1,3,5],:,s].>0.1,dims=1)) .&& (.!any(ErrorMa[j,i,[2,4,6],:,s].>0.1,dims=1)))
# 	# 		end
# 	# 	end
# 	# end
#     local MaFailCount2 = 0;
#     @inbounds for i = 1:num_sample
#         @inbounds for j = 1:num_angles
#             @inbounds for s = 1
#                 @inbounds for m =1:DIM_ERROR
#                     valid2 = (.!any(ErrorMaOne[j,i,[1,3,5],:,s].>0.2,dims=1)).&& (.!any(ErrorMaOne[j,i,[2,4,6],:,s].>0.2,dims=1))
#                     MaFailCount2 = MaFailCount2 + n_timesMa2 - sum(valid2)
#                     Err = ErrorMaOne[j,i,m,valid2[:],s]
#                     # Err = ErrorMa[j,i,m,:,1]
#                     E[j,i,m,n_methods+2] = Statistics.mean(Err)
#                     E[j,i,m,n_methods+2] = E[j,i,m,n_methods+2]*1000            # Convert to millimeter & 1e-3 rad for better readability 
#                 end
#             end
#         end
#     end
#     @printf "%d out of %d dataset on Ma fails" MaFailCount2/3 num_angles*n_times*size(n_sampleMa2,1)
# end


begin
    Title = [string(level," Noise ", s,mode," Err") for level in ["low","mid","high"] for s in ["X","Y","Z",""] for mode in [" Angle"," Dis"]]
    # Title = [string(level," Noise ",mode," Err") for level in ["low","mid","high"] for mode in [" Angle"," Dis"]]
    ylabels = [string(s,mode[1]," Err / ",mode[2]) for level in ["low","mid","high"] for s in ["X","Y","Z","AXB=YCZ"] for mode in [(" Rot"," 10^{-3} Rad"),(" Dis","mm")]]
    # ylabels = latexstring.(ylabels);

    # ylabels = [string(dir,mode) for dir in ["X","Y","Z"] for mode in [" Angle","Dis"]]
    s = 1;
    order = [1,5,2,6,3,7,4,8];
    for i = 1:3     # 3 noise levels
        # l = @layout [a b; c d; e f]
        l = @layout [a b c; d e f];
        Ps = Vector{Plots.Plot{Plots.GRBackend}}(undef,6)
        for j = 1:6     # 6 errors angle-x, transl-x, angle-y, transl-y, angle-z, transl-z
            if n_methods == 5    # Without GA Method
                if j == 2 #j == 2 || j == 4 || j == 6
                    # Wang
                    Ps[j] = plot(n_sample[s:2:end],E[i,s:2:end,j,2],linestyle=:dash,linecolor=RGB(109/255, 174/255, 209/255),
                    label = labels[2],linewidth=2,markershape=:utriangle,markercolor=RGB(109/255, 174/255, 209/255),markersize=5,#yscale=:log10,
                        ylabel = ylabels[(i-1)*DIM_ERROR+j],guidefont = (14,"times"), tickfont = (12), legendfontsize = 12, grid=true, legend=:right)
                    # Liao Wu
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,1],linestyle=:dashdotdot,linecolor=RGB(17/255, 70/255, 127/255),label = labels[1],
                        linewidth=2,markershape=:star5,markercolor=RGB(17/255, 70/255, 127/255),markersize=5)
                    # Fu+Bayro
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,4],linestyle=:dash,linecolor=RGB(245/255, 201/255, 91/255),label = labels[4],
                        linewidth=2,markershape=:hexagon,markercolor=RGB(245/255, 201/255, 91/255),markersize=5)
                    # Fu+Bayro 1/3
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,5],linestyle=:solid, linecolor=RGB(245/255, 201/255, 91/255),label = labels[5],
                        linewidth=2,markershape=:rect,markercolor=RGB(245/255, 201/255, 91/255),markersize=5)
                    # Ma
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,6],linestyle=:solid,linecolor=RGB(248/255, 178/255, 147/255),label = labels[6],
                        linewidth=2,markershape=:dtriangle,markercolor=RGB(248/255, 178/255, 147/255),markersize=5)
                    # G3
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,3],linestyle=:solid,linecolor=RGB(182/255, 35/255, 48/255),label = labels[3],
                        linewidth=2,markershape=:circle,markercolor=RGB(182/255, 35/255, 48/255),markersize=5)
                else
                    # Wang
                    Ps[j] = plot(n_sample[s:2:end],E[i,s:2:end,j,2],linestyle=:dash,linecolor=RGB(109/255, 174/255, 209/255),label = false,
                        linewidth=2,markershape=:utriangle,markercolor=RGB(109/255, 174/255, 209/255),markersize=5,#yscale=:log10,
                        ylabel = ylabels[(i-1)*DIM_ERROR+j],guidefont = (14,"times"), tickfont = (12), grid=true)
                    # Liao Wu
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,1],linestyle=:dashdotdot,linecolor=RGB(17/255, 70/255, 127/255),label = false,
                        linewidth=2,markershape=:star5,markercolor=RGB(17/255, 70/255, 127/255),markersize=5)
                    # Fu+Bayro
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,4],linestyle=:dash,linecolor=RGB(245/255, 201/255, 91/255),label = false,
                        linewidth=2,markershape=:hexagon,markercolor=RGB(245/255, 201/255, 91/255),markersize=5)
                    # Fu+Bayro1/3
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,5],linestyle=:solid, linecolor=RGB(245/255, 201/255, 91/255),label = false,
                        linewidth=2,markershape=:rect,markercolor=RGB(245/255, 201/255, 91/255),markersize=5)
                    # Ma
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,6],linestyle=:solid,linecolor=RGB(248/255, 178/255, 147/255),label = false,
                        linewidth=2,markershape=:dtriangle,markercolor=RGB(248/255, 178/255, 147/255),markersize=5)
                    # G3
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,3],linestyle=:solid,linecolor=RGB(182/255, 35/255, 48/255),label = false,
                        linewidth=2,markershape=:circle,markercolor=RGB(182/255, 35/255, 48/255),markersize=5)
                    # plot!(n_sample[s:2:end],E[i,s:2:end,j,7],linestyle=:dashdotdot,linecolor=RGB(182/255, 35/255, 48/255),label = false,
                    #     linewidth=2,markershape=:circle,markercolor=RGB(182/255, 35/255, 48/255),markersize=5)
                end
            else            # With GA Method
                if j == 7 #j == 2 || j == 4 || j == 6
                    # Wang
                    Ps[j] = plot(n_sample[s:2:end],E[i,s:2:end,j,2],linestyle=:dash,linecolor=RGB(109/255, 174/255, 209/255),
                    label = labels[2],linewidth=2,markershape=:utriangle,markercolor=RGB(109/255, 174/255, 209/255),markersize=5,#yscale=:log10,
                        ylabel = ylabels[(i-1)*DIM_ERROR+j],guidefont = (14,"times"), tickfont = (12), legendfontsize = 12, grid=true, legend=:right)
                    # Liao Wu
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,1],linestyle=:dashdotdot,linecolor=RGB(17/255, 70/255, 127/255),label = labels[1],
                        linewidth=2,markershape=:star5,markercolor=RGB(17/255, 70/255, 127/255),markersize=5)
                    # Fu+Bayro
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,4],linestyle=:dash,linecolor=RGB(245/255, 201/255, 91/255),label = labels[4],
                        linewidth=2,markershape=:hexagon,markercolor=RGB(245/255, 201/255, 91/255),markersize=5)
                    # Fu+Bayro 1/3
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,5],linestyle=:solid, linecolor=RGB(245/255, 201/255, 91/255),label = labels[5],
                        linewidth=2,markershape=:rect,markercolor=RGB(245/255, 201/255, 91/255),markersize=5)
                    # Ma
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,7],linestyle=:solid,linecolor=RGB(248/255, 178/255, 147/255),label = labels[6],
                        linewidth=2,markershape=:dtriangle,markercolor=RGB(248/255, 178/255, 147/255),markersize=5)
                    # GA
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,6],linestyle=:dash,linecolor=RGB(182/255, 35/255, 1/255),label = "GA",
                        linewidth=2,markershape=:rect,markercolor=RGB(182/255, 35/255, 1/255),markersize=5)
                    # G3
                    plot!(n_sample[s:2:end],E[i,s:2:end,j,3],linestyle=:solid,linecolor=RGB(182/255, 35/255, 48/255),label = labels[3],
                        linewidth=2,markershape=:circle,markercolor=RGB(182/255, 35/255, 48/255),markersize=5)
                else
                    # Wang
                    idxs = [1,3,5,2,4,6];
                    idx = idxs[j]; 
                    Ps[j] = plot(n_sample[s:2:end],E[i,s:2:end,idx,2],linestyle=:dash,linecolor=RGB(109/255, 174/255, 209/255),label = false,
                        linewidth=2,markershape=:utriangle,markercolor=RGB(109/255, 174/255, 209/255),markersize=5,#yscale=:log10,
                        ylabel = ylabels[(i-1)*DIM_ERROR+idx],guidefont = (14,"times"), tickfont = (12), grid=true)
                    # Liao Wu
                    plot!(n_sample[s:2:end],E[i,s:2:end,idx,1],linestyle=:dashdotdot,linecolor=RGB(17/255, 70/255, 127/255),label = false,
                        linewidth=2,markershape=:star5,markercolor=RGB(17/255, 70/255, 127/255),markersize=5)
                    # Fu+Bayro
                    plot!(n_sample[s:2:end],E[i,s:2:end,idx,4],linestyle=:dash,linecolor=RGB(245/255, 201/255, 91/255),label = false,
                        linewidth=2,markershape=:hexagon,markercolor=RGB(245/255, 201/255, 91/255),markersize=5)
                    # Fu+Bayro1/3
                    plot!(n_sample[s:2:end],E[i,s:2:end,idx,5],linestyle=:solid, linecolor=RGB(245/255, 201/255, 91/255),label = false,
                        linewidth=2,markershape=:rect,markercolor=RGB(245/255, 201/255, 91/255),markersize=5)
                    # Ma
                    # plot!(n_sample[s:2:end],E[i,s:2:end,idx,7],linestyle=:solid,linecolor=RGB(248/255, 178/255, 147/255),label = false,
                    #     linewidth=2,markershape=:dtriangle,markercolor=RGB(248/255, 178/255, 147/255),markersize=5)
                    # GA
                    plot!(n_sample[s:2:end],E[i,s:2:end,idx,6],linestyle=:dash,linecolor=RGB(182/255, 35/255, 1/255),label = false,
                        linewidth=2,markershape=:rect,markercolor=RGB(182/255, 35/255, 1/255),markersize=5)
                    # G3
                    plot!(n_sample[s:2:end],E[i,s:2:end,idx,3],linestyle=:solid,linecolor=RGB(182/255, 35/255, 48/255),label = false,
                        linewidth=2,markershape=:circle,markercolor=RGB(182/255, 35/255, 48/255),markersize=5)
                end
            end
        end
        plot(Ps...,layout=l,size=(800,500))
        if n_methods == 5
            savefig(string("./images/",name,"/",levels[i],".svg"))
        else
            savefig(string("./images/",name,"/",levels[i],"_GA.svg"))
        end
    end
end

## Evaluate Result from AXB(YCZ)^-1
## This error ceriteria is not used in the paper, because it does not truthfully reflect the Result
## The error from A,B,C caused the fluctuation of the error.
begin
    Title2 = [string(level," Noise ", s,mode," Err") for level in ["low","mid","high"] for s in ["X","Y","Z",""] for mode in [" Angle"," Dis"]]
    # Title = [string(level," Noise ",mode," Err") for level in ["low","mid","high"] for mode in [" Angle"," Dis"]]
    ylabels2 = [string(s,mode[1]," Err / ",mode[2]) for level in ["low","mid","high"] for s in ["AXB=YCZ"] for mode in [(" Rot"," 10^{-3} Rad"),(" Dis","mm")]]
    # ylabels = latexstring.(ylabels);

    # ylabels = [string(dir,mode) for dir in ["X","Y","Z"] for mode in [" Angle","Dis"]]
    i = 2
    Ps_ = Vector{Plots.Plot{Plots.GRBackend}}(undef,2)
    for j = 7:8    # 2 errors under AXBYCZ cross check
            Ps_[j-6] = plot(n_sample[s:2:end],E[i,s:2:end,j,2],linestyle=:dash,linecolor=RGB(109/255, 174/255, 209/255),label = false,
                linewidth=2,markershape=:utriangle,markercolor=RGB(109/255, 174/255, 209/255),markersize=5,
                ylabel = ylabels2[j-6],guidefont = (14,"times"), tickfont = (12), grid=true)
            plot!(n_sample[s:2:end],E[i,s:2:end,j,1],linestyle=:dashdotdot,linecolor=RGB(17/255, 70/255, 127/255),label = false,
                linewidth=2,markershape=:star5,markercolor=RGB(17/255, 70/255, 127/255),markersize=5)
            plot!(n_sample[s:2:end],E[i,s:2:end,j,4],linestyle=:dash,linecolor=RGB(245/255, 201/255, 91/255),label = false,
                linewidth=2,markershape=:hexagon,markercolor=RGB(245/255, 201/255, 91/255),markersize=5)
            plot!(n_sample[s:2:end],E[i,s:2:end,j,5],linestyle=:solid, linecolor=RGB(245/255, 201/255, 91/255),label = false,
                linewidth=2,markershape=:rect,markercolor=RGB(245/255, 201/255, 91/255),markersize=5)
            plot!(n_sample[s:2:end],E[i,s:2:end,j,3],linestyle=:solid,linecolor=RGB(182/255, 35/255, 48/255),label = false,
                linewidth=2,markershape=:circle,markercolor=RGB(182/255, 35/255, 48/255),markersize=5)
            # plot!(n_sample[s:2:end],E[i,s:2:end,j,6],linestyle=:solid,linecolor=RGB(248/255, 178/255, 147/255),label = false,
            #     linewidth=2,markershape=:dtriangle,markercolor=RGB(248/255, 178/255, 147/255),markersize=5)
    end
    plot(Ps_...,layout=(@layout [a b]),size=(800,300))
end


savefig(string("./images/",name,"/AXBYCZ_",levels[i],".svg"))

end

main(ARGS)