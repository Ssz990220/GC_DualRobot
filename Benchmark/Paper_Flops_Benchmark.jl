name = "Benchmark_Flops"
labels = ["Wu" "Wang" "Proposed"]
discribe = "Benchmark: Computation Openations"


begin
	using Test, StaticArrays, BenchmarkTools, LIKWID, CurveFit, Printf
	include("../project_utils/DualRobotCalibrator.jl")
	using .DualRobotCalibrator
	using .DualRobotCalibrator: iter_solve, get_i_grad_reg_AXBYCZG3, get_full_grad_regG3, ErrG3ₙ, E_normG3, toVG3
end
unit = "mm"
begin
    X = SMatrix{4,4,Float64,16}([1. 0 0 0;0 1 0 0;0 0 1 270;0 0 0 1]*tfrotx(π/3)*tfroty(π/6));
	Z = tfxyz([0.,0.,100])*tfrotx(π)*tfroty(0.1)
	Y = tfxyz([1700.,0,200])*tfrotz(π)*tfroty(-π/6)*tfrotx(π/12);
	if unit == "mm"
	else
		X = tometer(X); Y = tometer(Y); Z = tometer(Z)
	end
	Xₜ = tfm2motorS(X); Yₜ = tfm2motorS(Y); Zₜ = tfm2motorS(Z);
	Ans = [X,Y,Z]; Ansₘ = [Xₜ,Yₜ,Zₜ];
end

function get_VPₘ(n)
	Aₘₛ,Bₘₛ,Cₘₛ = get_ABCSₘ(n,Xₜ,Yₜ,Zₜ,unit);
	Aₙₘₛ = add_noiseSₘ.(Aₘₛ,angle_noise,dis_noise);
	Bₙₘₛ = add_noiseSₘ.(Bₘₛ,angle_noise,dis_noise);
	Cₙₘₛ = add_noiseSₘ.(Cₘₛ,angle_noise,dis_noise);
	Xᵢₙᵢₜₘ = add_noiseSₘ(Xₜ,xyz_noise_angle,xyz_dis_noise[1]);
	Yᵢₙᵢₜₘ = add_noiseSₘ(Yₜ,xyz_noise_angle,xyz_dis_noise[2]);
	Zᵢₙᵢₜₘ = add_noiseSₘ(Zₜ,xyz_noise_angle,xyz_dis_noise[3]);
	Initₘ = [Xᵢₙᵢₜₘ,Yᵢₙᵢₜₘ,Zᵢₙᵢₜₘ];
	V = scale_motorS.(match_signS.(Ansₘ,Initₘ),scale);
	Pₘ = [scale_motor.(Aₙₘₛ,scale),scale_motorS.(Bₙₘₛ,scale),scale_motor.(Cₙₘₛ,scale)];
	V = toVᵥ(motor2vectorS.(V))
	return V,Pₘ	
end

function get_VPG3(n)
	Aₘₛ,Bₘₛ,Cₘₛ = get_ABCSₘ(n,Xₜ,Yₜ,Zₜ,unit);
	Aₙₘₛ = add_noiseSₘ.(Aₘₛ,angle_noise,dis_noise);
	Bₙₘₛ = add_noiseSₘ.(Bₘₛ,angle_noise,dis_noise);
	Cₙₘₛ = add_noiseSₘ.(Cₘₛ,angle_noise,dis_noise);
	Xᵢₙᵢₜₘ = add_noiseSₘ(Xₜ,xyz_noise_angle,xyz_dis_noise[1]);
	Yᵢₙᵢₜₘ = add_noiseSₘ(Yₜ,xyz_noise_angle,xyz_dis_noise[2]);
	Zᵢₙᵢₜₘ = add_noiseSₘ(Zₜ,xyz_noise_angle,xyz_dis_noise[3]);
	Initₘ = [Xᵢₙᵢₜₘ,Yᵢₙᵢₜₘ,Zᵢₙᵢₜₘ];
	V = scale_motorS.(match_signS.(Ansₘ,Initₘ),scale);
	Pₘ = [Aₙₘₛ,Bₙₘₛ,Cₙₘₛ];
	V = toVG3(motorS2rotorG3S.(V))
	Pₘ = [motorS2rotorG3S.(Pₘ[i]) for i∈eachindex(Pₘ)]
	return V,Pₘ	
end

function get_VP(n)
	Aₘₛ,Bₘₛ,Cₘₛ = get_ABCSₘ(n,Xₜ,Yₜ,Zₜ,unit);
	Aₙₘₛ = add_noiseSₘ.(Aₘₛ,angle_noise,dis_noise);
	Bₙₘₛ = add_noiseSₘ.(Bₘₛ,angle_noise,dis_noise);
	Cₙₘₛ = add_noiseSₘ.(Cₘₛ,angle_noise,dis_noise);
	Xᵢₙᵢₜₘ = add_noiseSₘ(Xₜ,xyz_noise_angle,xyz_dis_noise[1]);
	Yᵢₙᵢₜₘ = add_noiseSₘ(Yₜ,xyz_noise_angle,xyz_dis_noise[2]);
	Zᵢₙᵢₜₘ = add_noiseSₘ(Zₜ,xyz_noise_angle,xyz_dis_noise[3]);
	Aₙ₂ = motorS2tfm.(Aₙₘₛ);Bₙ₂ = motorS2tfm.(Bₙₘₛ);Cₙ₂ = motorS2tfm.(Cₙₘₛ);
	Xᵢₙᵢₜ₂ = motorS2tfm(Xᵢₙᵢₜₘ); Yᵢₙᵢₜ₂ = motorS2tfm(Yᵢₙᵢₜₘ); Zᵢₙᵢₜ₂ = motorS2tfm(Zᵢₙᵢₜₘ);
	Init = [Xᵢₙᵢₜ₂,Yᵢₙᵢₜ₂,Zᵢₙᵢₜ₂]; P₂ = [Aₙ₂,Bₙ₂,Cₙ₂];
	return Init,P₂
end

function get_VPₗ(n)
	V,P = get_VP(n);
	RA = [P[1][i][SOneTo(3),SOneTo(3)] for i = 1:n];
	RB = [P[2][i][SOneTo(3),SOneTo(3)] for i = 1:n];
	RC = [P[3][i][SOneTo(3),SOneTo(3)] for i = 1:n];
	ta = [P[1][i][SOneTo(3),4] for i = 1:n];
	tb = [P[2][i][SOneTo(3),4] for i = 1:n];
	tc = [P[3][i][SOneTo(3),4] for i = 1:n];
	RXᵢₙᵢₜ = V[1][SOneTo(3),SOneTo(3)];
	RYᵢₙᵢₜ = V[2][SOneTo(3),SOneTo(3)];
	RZᵢₙᵢₜ = V[3][SOneTo(3),SOneTo(3)];
	return RA,RB,RC,RXᵢₙᵢₜ,RYᵢₙᵢₜ,RZᵢₙᵢₜ,ta,tb,tc
end	


begin
    scale = 1
    N_iter = 100;
    angle_noise = 0.3/180*π; dis_noise = 0.3; #Parameters
    xyz_noise_angle = 0.5/180*π; xyz_dis_noise = ones(1,3)*3.0;
end
# begin
# 	conf_gc = ((max_iter = N_iter, m = 3, η=0.1,μ=0.1,scale=1e-2, err = "NAN", stop = false, reg=true,τ=0.3))
#     println("GC_Katyusha computation Time:")
# 	@btime katyushaX(Vₛ,Pₛ,$get_i_grad_reg_AXBYCZᵥ₂, $get_full_grad_regᵥ₂, $Errₙ,$normfuncᵥ₂, $conf_gc) setup=(begin; Vₛ,Pₛ=get_VPₘ(256); end;)
# end;

begin
	nData = 30:30:300
	Flops = zeros(4,size(nData,1))
	conf_wang = ((max_iter = N_iter, μ = [1,1e-5,1,1,1],m=3, η=0.1,stop = false,err="NAN",reg = true,τ=0.3))
	conf_g3ₖ = ((max_iter = N_iter, m = 3, η=0.3,μ=0.1,scale = 1e-2, stop = false, err="NAN",τ=0.3,reg=true))
	confₗ = ((max_iter = N_iter,))

    ## Warm Up
    n_ = 30;
    Vₛ,Pₛ=get_VP(n_)
    @perfmon "FLOPS_DP" iter_solve(Vₛ,Pₛ,conf_wang);
    VG3ₛ,PG3ₛ=get_VPG3(n_)
    @perfmon "FLOPS_DP" katyushaX(VG3ₛ,PG3ₛ,get_i_grad_reg_AXBYCZG3,get_full_grad_regG3,ErrG3ₙ,E_normG3,conf_g3ₖ);
    RA,RB,RC,RXᵢₙᵢₜ,RYᵢₙᵢₜ,RZᵢₙᵢₜ=get_VPₗ(n_)
    @perfmon "FLOPS_DP" liao_iter(RA,RB,RC,RXᵢₙᵢₜ,RYᵢₙᵢₜ,RZᵢₙᵢₜ,confₗ);
    RA,RB,RC,RX_sln, RY_sln, RZ_sln,ta,tb,tc=get_VPₗ(n_)
    @perfmon "FLOPS_DP" T_SVD(RX_sln,RY_sln,RZ_sln,RA,RB,RC,ta,tb,tc);
	for i ∈ eachindex(nData);
        n = nData[i];

        Vₛt,Pₛt=get_VP(n)
        metrics, events = @perfmon "FLOPS_DP" iter_solve(Vₛt,Pₛt,conf_wang);
        flops_per_second = first(metrics["FLOPS_DP"])["DP [MFLOP/s]"] * 1e6
        runtime = first(metrics["FLOPS_DP"])["Runtime (RDTSC) [s]"]
        NFLOPs_actual = round(Int, flops_per_second * runtime)
        Flops[1,i] = NFLOPs_actual;


        VG3ₛt,PG3ₛt=get_VPG3(n)
        metrics, events = @perfmon "FLOPS_DP" katyushaX(VG3ₛt,PG3ₛt,get_i_grad_reg_AXBYCZG3,get_full_grad_regG3,ErrG3ₙ,E_normG3,conf_g3ₖ);
        flops_per_second = first(metrics["FLOPS_DP"])["DP [MFLOP/s]"] * 1e6
        runtime = first(metrics["FLOPS_DP"])["Runtime (RDTSC) [s]"]
        NFLOPs_actual = round(Int, flops_per_second * runtime)
        Flops[2,i] = NFLOPs_actual;

        RAt,RBt,RCt,RXᵢₙᵢₜt,RYᵢₙᵢₜt,RZᵢₙᵢₜt=get_VPₗ(n)
        metrics, events = @perfmon "FLOPS_DP" liao_iter(RAt,RBt,RCt,RXᵢₙᵢₜt,RYᵢₙᵢₜt,RZᵢₙᵢₜt,confₗ);
        flops_per_second = first(metrics["FLOPS_DP"])["DP [MFLOP/s]"] * 1e6
        runtime = first(metrics["FLOPS_DP"])["Runtime (RDTSC) [s]"]
        NFLOPs_actual = round(Int, flops_per_second * runtime)
        Flops[3,i] = NFLOPs_actual;

        RAt,RBt,RCt,RX_slnt, RY_slnt, RZ_slnt,tat,tbt,tct=get_VPₗ(n)
        metrics, events = @perfmon "FLOPS_DP" T_SVD(RX_slnt,RY_slnt,RZ_slnt,RAt,RBt,RCt,tat,tbt,tct);
        flops_per_second = first(metrics["FLOPS_DP"])["DP [MFLOP/s]"] * 1e6
        runtime = first(metrics["FLOPS_DP"])["Runtime (RDTSC) [s]"]
        NFLOPs_actual = round(Int, flops_per_second * runtime)
        Flops[4,i] = NFLOPs_actual;
	end
end;


begin
	using JLD, HDF5
    mkpath(string("./images/",name))
    save(string("./images/",name,"/Result.jld"),"Flops",Flops)
end
begin
	Flops[2,:].=Flops[2,:].+Flops[4,:];
	Flops[3,:].=Flops[3,:].+Flops[4,:];
	using Plots
	plot(nData,Flops[1,:].*1e-6,linestyle=:dash,linecolor=RGB(109/255, 174/255, 209/255),
		label = false,linewidth=2,markershape=:utriangle,markercolor=RGB(109/255, 174/255, 209/255),markersize=5,
		ylabel = "FLOPs / 10^{6}",xlabel = "Data Set Size", guidefont = (14,"times"), tickfont = (12), legendfontsize = 14, grid=true, legend=:right)
	plot!(nData,Flops[3,:].*1e-6,linestyle=:dashdotdot,linecolor=RGB(17/255, 70/255, 127/255),label = false,
		linewidth=2,markershape=:star5,markercolor=RGB(17/255, 70/255, 127/255),markersize=5)
	plot!(nData,Flops[2,:].*1e-6,linestyle=:solid,linecolor=RGB(182/255, 35/255, 48/255),label = false,
		linewidth=2,markershape=:circle,markercolor=RGB(182/255, 35/255, 48/255),markersize=5)
	savefig(string("./images/",name,"/","Flops.svg"))
end

begin
	a,b = linear_fit(nData,Flops[2,:])
	@printf "Proposed method: NFlops ≈ %.2f e5 * nData + %.2f e6\n" b/1e5 a/1e6
	a,b = linear_fit(nData,Flops[1,:])
	@printf "Wang's method: NFlops ≈ %.2f e5 * nData + %.2f e6\n" b/1e5 a/1e6
	a,b = linear_fit(nData,Flops[3,:])
	@printf "Wu's method: NFlops ≈ %.2f e5 * nData + %.2f e6\n" b/1e5 a/1e6
end