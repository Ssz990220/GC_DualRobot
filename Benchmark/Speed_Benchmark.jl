### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ c7f66c54-6e8d-4dcf-8e9a-9996f230309d
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.activate("../.")
	using PlutoUI, Reexport, StaticArrays, Plots, Printf, Statistics, Distributions
	using LinearAlgebra, Random, Test, BenchmarkTools, StatsPlots, StatsPlots.PlotMeasures
end

# ╔═╡ dc0e522e-f3dc-4bf6-96f5-d402f09fd4e6
# ╠═╡ show_logs = false
begin
	include("../project_utils/DualRobotCalibrator.jl")
	@reexport  using .DualRobotCalibrator
	import .DualRobotCalibrator: iter_solve as iter_wang,AXB_YCZ_Close, post_process_Wang, err,iter_solve_Katyusha as iter_wang_k
	import .DualRobotCalibrator: toVG3, toXYZ as toXYZG3, ErrG3, iter_solve_Katyusha_G3
	md"Load Packages"
end

# ╔═╡ 209b3ee9-69a0-4f6b-aefe-4af20d2e583d
md"
# Benchmark: Wang, Wu VS Proposed in convergence rate and range

**This notebook compares Wang's & Wu's iterative method against G3 head to head.** Code optimization has been done as much as I can.

It will take several minutes to install required package due to pluto package management protocal. Loading time can be signaficantly reduced if you preinstalled listed packages in default julia enviornment.

But around one minute is still needed to precompile packages and included libraries everytime you open the notebook.
"

# ╔═╡ 1c7f849b-4f1c-421f-bf34-ce2dcd46e955
md"## Parameters"

# ╔═╡ 94499cf7-66b3-4920-af7b-0bfeaa3f658b
unit = "m"

# ╔═╡ c338bab3-7d44-400c-8aea-e5e7820d9ba6
begin
	X = SMatrix{4,4,Float64,16}([1. 0 0 0;0 1 0 0;0 0 1 270;0 0 0 1]*tfrotx(π/3)*tfroty(π/6));
	Z = tfxyz([0.,0.,100])*tfrotx(π)#*tfroty(0.1)
	Y = tfxyz([1700.,0,200])*tfrotz(π)#*tfroty(-π/6)*tfrotx(π/12);
	if unit == "mm"
	else
		X = tometer(X); Y = tometer(Y); Z = tometer(Z)
	end
	Xₜ = tfm2motorS(X); Yₜ = tfm2motorS(Y); Zₜ = tfm2motorS(Z);
	Ans = [X,Y,Z]; Ansₘ = [Xₜ,Yₜ,Zₜ];
	md"True Values Here"
end

# ╔═╡ bf930527-f561-4537-904b-c2cba3474b54
md"## Descent Speed Test

This part shows how efficient each algorithm is for each iteration it runs."

# ╔═╡ 77683183-1f55-4436-ab26-934e452625c3
begin
	N_slider = @bind N_iter Slider(10:5:200,30, true)
	md"Number of Iteration: $(N_slider) "
end

# ╔═╡ 61d732b0-19f3-4d34-8fd3-ab224300f828
begin
	israndom = @bind random CheckBox()
	md"Random Initial Value? $(israndom)"
end

# ╔═╡ f492eecf-b9a1-47a5-873c-1b6f23088664
begin
	isshared = @bind shared CheckBox(default=true)
	md"Shared Initial Value? $(isshared)"
end

# ╔═╡ 95692758-52e4-4d84-818d-53abc3209fe0
begin
	r_view = @bind r Slider(2:20,11,true)
	md"Number of steps in zoomed view : $(r_view), \
	available when Random Initial Value is true"
end

# ╔═╡ ca5f77b1-9b44-4496-9f53-8dd04467f934
md"Wang's iterative method along does not refine its rotational component. This fig is not shown in the paper. Translational component for Wu and proposed method is solved in each step by SVD."

# ╔═╡ 56079b1e-a80b-4818-8fbd-b2c0138185c6
begin
	n = 32;
	samples = sample(1:n,32,replace=false);
	# samples = 1:n
	angle_noise = 1.2/180*π; dis_noise = 1.2; # mm
	dis_noise = (unit=="mm") ? dis_noise : dis_noise/1000;
	xyz_noise_angle = 0.3/180*π; xyz_dis_noise = [5.0,5.0,5.0]; # mm
	xyz_dis_noise = (unit=="mm") ? xyz_dis_noise : xyz_dis_noise/1000;
	md"Data Gen Config, 
	$n measurements are simulated, angle noise level is $angle_noise rad, translational noise level is $dis_noise $unit"
end

# ╔═╡ 3f2be549-1fe8-45c8-9b23-c33aadb3ed93
begin
	Aₘₛ,Bₘₛ,Cₘₛ = get_ABCSₘ(n,Xₜ,Yₜ,Zₜ,unit);
	Aₙₘₛ = add_noiseSₘ.(Aₘₛ,angle_noise,dis_noise,true);
	Bₙₘₛ = add_noiseSₘ.(Bₘₛ,angle_noise,dis_noise,true);
	Cₙₘₛ = add_noiseSₘ.(Cₘₛ,angle_noise,dis_noise,true);
	Pₘ = [Aₙₘₛ,Bₙₘₛ,Cₙₘₛ];

	## For Matrix
	Aₙ₂ = motorS2tfm.(Aₙₘₛ);Bₙ₂ = motorS2tfm.(Bₙₘₛ);Cₙ₂ = motorS2tfm.(Cₙₘₛ);
	P₂ = [Aₙ₂,Bₙ₂,Cₙ₂];
	md"Data Gen"
end

# ╔═╡ 71792be6-7a13-4f8d-aead-c327f1b4fc36
begin
	RA = [Aₙ₂[i][SOneTo(3),SOneTo(3)] for i = 1:n];
	RB = [Bₙ₂[i][SOneTo(3),SOneTo(3)] for i = 1:n];
	RC = [Cₙ₂[i][SOneTo(3),SOneTo(3)] for i = 1:n];
	ta = [Aₙ₂[i][SOneTo(3),4] for i = 1:n];
	tb = [Bₙ₂[i][SOneTo(3),4] for i = 1:n];
	tc = [Cₙ₂[i][SOneTo(3),4] for i = 1:n];
end;

# ╔═╡ b1ca5ec0-d291-49ed-b898-43da84854477
md"## Descent Test with 0 init translation value"

# ╔═╡ 591c9131-f99e-421a-be3c-92e0b379a2f9
@bind xyz_noise_angle₂ Slider(0.0:0.5:10,default=3)

# ╔═╡ 060f90e0-9d9d-4557-9dfb-5ce6aecea51a
md"Angle Noise Added to init Values: $xyz_noise_angle₂ °"

# ╔═╡ b7e4e47b-2c43-4e5c-b74a-c10c6062225c
md"Random Initial Value in All SO3 space: $(@bind randSO3 CheckBox())"

# ╔═╡ ca204e58-57a9-4a96-acd5-c08497e13317
begin
	N_slider₂ = @bind N_iter₂ Slider(10:10:200,20, true)
	md"Number of Iteration: $(N_slider₂) "
end

# ╔═╡ 43ec9469-619e-4467-876e-e49f3f8955e5
md"## Functions"

# ╔═╡ 6a8b84bf-1de4-4f7b-99e8-4ca87021202f
function rand_SE3()
	R = SMatrix(rand(AngleAxis));
	T = @SMatrix rand(3,1);
	unit == "mm" ? T = T*100 : T = T
	return rt2T(R,T)
end

# ╔═╡ b7392bcf-d3cc-4838-879f-657d385c9037
begin
	if randSO3
		Xᵢₙᵢₜ₃ = rand_SE3(); Yᵢₙᵢₜ₃ = rand_SE3(); Zᵢₙᵢₜ₃ = rand_SE3();
		RXᵢₙᵢₜ₀ = Xᵢₙᵢₜ₃[SOneTo(3),SOneTo(3)];
		RYᵢₙᵢₜ₀ = Yᵢₙᵢₜ₃[SOneTo(3),SOneTo(3)];
		RZᵢₙᵢₜ₀ = Zᵢₙᵢₜ₃[SOneTo(3),SOneTo(3)];
		Xᵢₙᵢₜₘ₂ = tfm2motorS(Xᵢₙᵢₜ₃); 
		rand() > 0.5 ? Xᵢₙᵢₜₘ₂ = Xᵢₙᵢₜₘ₂ : Xᵢₙᵢₜₘ₂ = -Xᵢₙᵢₜₘ₂
		Yᵢₙᵢₜₘ₂ = tfm2motorS(Yᵢₙᵢₜ₃);
		rand() > 0.5 ? Yᵢₙᵢₜₘ₂ = Yᵢₙᵢₜₘ₂ : Yᵢₙᵢₜₘ₂ = -Yᵢₙᵢₜₘ₂
		Zᵢₙᵢₜₘ₂ = tfm2motorS(Zᵢₙᵢₜ₃);
		rand() > 0.5 ? Zᵢₙᵢₜₘ₂ = Zᵢₙᵢₜₘ₂ : Zᵢₙᵢₜₘ₂ = -Zᵢₙᵢₜₘ₂
		Init₂ = [Xᵢₙᵢₜ₃,Yᵢₙᵢₜ₃,Zᵢₙᵢₜ₃]
		Initₘ₂ = motorS2rotorG3S.([Xᵢₙᵢₜₘ₂,Yᵢₙᵢₜₘ₂,Zᵢₙᵢₜₘ₂])
	else
		Xᵢₙᵢₜₘ₂ = add_noiseSₘ(Xₜ,xyz_noise_angle₂/180*π,xyz_dis_noise[1],true);
		Yᵢₙᵢₜₘ₂ = add_noiseSₘ(Yₜ,xyz_noise_angle₂/180*π,xyz_dis_noise[2],true);
		Zᵢₙᵢₜₘ₂ = add_noiseSₘ(Zₜ,xyz_noise_angle₂/180*π,xyz_dis_noise[3],true);
		Xᵢₙᵢₜₘ₂ = motorS2rotorG3S(Xᵢₙᵢₜₘ₂);Yᵢₙᵢₜₘ₂ = motorS2rotorG3S(Yᵢₙᵢₜₘ₂);Zᵢₙᵢₜₘ₂ = motorS2rotorG3S(Zᵢₙᵢₜₘ₂);
		Initₘ₂ = [Xᵢₙᵢₜₘ₂,Yᵢₙᵢₜₘ₂,Zᵢₙᵢₜₘ₂];
	
		## For Matrix
		Xᵢₙᵢₜ₃ = r2t(rotorG32rm(Xᵢₙᵢₜₘ₂)); Yᵢₙᵢₜ₃ = r2t(rotorG32rm(Yᵢₙᵢₜₘ₂)); Zᵢₙᵢₜ₃ = r2t(rotorG32rm(Zᵢₙᵢₜₘ₂));
		Init₂ = [Xᵢₙᵢₜ₃,Yᵢₙᵢₜ₃,Zᵢₙᵢₜ₃];
		RXᵢₙᵢₜ₀ = Xᵢₙᵢₜ₃[SOneTo(3),SOneTo(3)];
		RYᵢₙᵢₜ₀ = Yᵢₙᵢₜ₃[SOneTo(3),SOneTo(3)];
		RZᵢₙᵢₜ₀ = Zᵢₙᵢₜ₃[SOneTo(3),SOneTo(3)];
	end
	md"Gen Data"
end

# ╔═╡ ba20da2e-c3c5-416e-ace2-3d2555aebbcb
errᵢₙᵢₜ₀ = get_error_XYZ(Init₂,Ansₘ)

# ╔═╡ 20fccdd0-1f3d-421a-8fe2-a13176df32c2
function get_VP(n)
	Aₘₛ,Bₘₛ,Cₘₛ = get_ABCₘ(n,Xₜ,Yₜ,Zₜ,unit);
	Aₙₘₛ = add_noiseₘ.(Aₘₛ,angle_noise,dis_noise);
	Bₙₘₛ = add_noiseₘ.(Bₘₛ,angle_noise,dis_noise);
	Cₙₘₛ = add_noiseₘ.(Cₘₛ,angle_noise,dis_noise);
	Xᵢₙᵢₜₘ = add_noiseₘ(Xₜ,xyz_noise_angle,xyz_dis_noise[1]);
	Yᵢₙᵢₜₘ = add_noiseₘ(Yₜ,xyz_noise_angle,xyz_dis_noise[2]);
	Zᵢₙᵢₜₘ = add_noiseₘ(Zₜ,xyz_noise_angle,xyz_dis_noise[3]);
	Aₙ₂ = motor2tfm.(Aₙₘₛ);Bₙ₂ = motor2tfm.(Bₙₘₛ);Cₙ₂ = motor2tfm.(Cₙₘₛ);
	Xᵢₙᵢₜ₂ = motor2tfm(Xᵢₙᵢₜₘ); Yᵢₙᵢₜ₂ = motor2tfm(Yᵢₙᵢₜₘ); Zᵢₙᵢₜ₂ = motor2tfm(Zᵢₙᵢₜₘ);
	Init = [Xᵢₙᵢₜ₂,Yᵢₙᵢₜ₂,Zᵢₙᵢₜ₂]; P₂ = [Aₙ₂,Bₙ₂,Cₙ₂];
	return Init,P₂
end

# ╔═╡ d7f4c0f3-296e-49cd-b453-bf0778a1608f
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

# ╔═╡ ba87cb0a-697a-462c-941a-ab9605d35a16
function get_VPₛ(n)
	RA,RB,RC,RXᵢₙᵢₜ,RYᵢₙᵢₜ,RZᵢₙᵢₜ,ta,tb,tc = get_VPₗ(n);
	RX_sln,RY_sln,RZ_sln = liao_iter(RA,RB,RC,RXᵢₙᵢₜ,RYᵢₙᵢₜ,RZᵢₙᵢₜ,((max_iter=200,)))
	return RX_sln, RY_sln, RZ_sln, RA,RB,RC,ta,tb,tc
end

# ╔═╡ ad34b9bf-7162-4682-aa71-7ec2aadefd27
function get_Pₘ(n)
	Aₘₛ,Bₘₛ,Cₘₛ = get_ABCₘ(n,Xₜ,Yₜ,Zₜ,unit);
	return [motor2rotor.(Aₘₛ),motor2rotor.(Bₘₛ),motor2rotor.(Cₘₛ)]
end

# ╔═╡ 3269f5e3-c216-4a27-a6ea-4e9cc297dced
function T_Wang(X,Y,Z)
	return T_SVD(X[SOneTo(3),SOneTo(3)],Y[SOneTo(3),SOneTo(3)],Z[SOneTo(3),SOneTo(3)],RA,RB,RC,ta,tb,tc)
end

# ╔═╡ ebc58e1b-fae3-4e1e-9451-b150417b76ca
begin
	function ErrG3_Full(V,P,∇ᵥgṼₜ,ηₜ,conf)
	    Result = toXYZG3(V);
		Xᵣ = rotorG32rm(Result[1]);
		Yᵣ = rotorG32rm(Result[2]);
		Zᵣ = rotorG32rm(Result[3]);
		Xᵣ,Yᵣ,Zᵣ = T_Wang(Xᵣ,Yᵣ,Zᵣ);
	    Errx = get_error(Xᵣ,Xₜ);
	    Erry = get_error(Yᵣ,Yₜ);
	    Errz = get_error(Zᵣ,Zₜ);
	    norme = norm(∇ᵥgṼₜ)
	    err = ErrG3(V,P,∇ᵥgṼₜ,ηₜ,conf)
	    return [Errx[1] Errx[2] Erry[1] Erry[2] Errz[1] Errz[2] ηₜ norme err];
	end
	function err_Wang_Full(V,P,∇ᵥgṼₜ,ηₜ,conf)
		Xᵣ,Yᵣ,Zᵣ = post_process_Wang(V)
		A = P[1]; B = P[2]; C = P[3];
		Xᵣ = tometer(Xᵣ); Yᵣ = tometer(Yᵣ); Zᵣ = tometer(Zᵣ);
		# Xᵣ,Yᵣ,Zᵣ = T_Wang(Xᵣ,Yᵣ,Zᵣ)
		Errx = get_error(Xᵣ,X);
		Erry = get_error(Yᵣ,Y);
		Errz = get_error(Zᵣ,Z);
		norme = norm(∇ᵥgṼₜ,2)
		error = sum(err.(A,B,C,[Xᵣ],[Yᵣ],[Zᵣ]))/size(P[1],1)
		return [Errx[1] Errx[2] Erry[1] Erry[2] Errz[1] Errz[2] ηₜ norme error]
	end
	function Err_Liao(V)
		Xᵣ,Yᵣ,Zᵣ = T_Wang(V[1],V[2],V[3]);
		Xᵣ = SMatrix{4,4,Float64,16}([ProjectToSO3(Xᵣ[SOneTo(3),SOneTo(3)]) Xᵣ[SOneTo(3),4];zeros(1,4)])
		Yᵣ = SMatrix{4,4,Float64,16}([ProjectToSO3(Yᵣ[SOneTo(3),SOneTo(3)]) Yᵣ[SOneTo(3),4];zeros(1,4)])
		Zᵣ = SMatrix{4,4,Float64,16}([ProjectToSO3(Zᵣ[SOneTo(3),SOneTo(3)]) Zᵣ[SOneTo(3),4];zeros(1,4)])
		# Xᵣ = [Xᵣ[SOneTo(3),SOneTo(3)] Xᵣ[SOneTo(3),4];zeros(1,4)]
		# Yᵣ = [Yᵣ[SOneTo(3),SOneTo(3)] Yᵣ[SOneTo(3),4];zeros(1,4)]
		# Zᵣ = [Zᵣ[SOneTo(3),SOneTo(3)] Zᵣ[SOneTo(3),4];zeros(1,4)]
		Errx = get_error(Xᵣ,X);
		Erry = get_error(Yᵣ,Y);
		Errz = get_error(Zᵣ,Z);
		return [Errx[1] Errx[2] Erry[1] Erry[2] Errz[1] Errz[2] 0.0 0.0 0.0]
	end
	md"Error Functions"
end

# ╔═╡ 26c72e43-07bb-4381-b96e-bec3519a0469
begin
	scale = 1e-2;
	conf = ((max_iter = N_iter, μ = [1,1e-5,1,1,1],m=3, η=0.3,stop = false,err="Ext",errfunc = err_Wang_Full,errdim=9,τ=0.3,reg = true,svd=false))
	confₘ₃ = ((max_iter = N_iter, m = 3, η=0.3,μ=0.1,scale = scale, stop = false, err="Ext",errfunc = ErrG3_Full,errdim=9,τ=0.3,reg = true))
	conf₅ = ((max_iter = N_iter,err="Ext",errfunc = Err_Liao,errdim=9))
	md"Solver Config"
end

# ╔═╡ 0a661256-be96-4359-9515-eff2cc0ebbd5
begin
	if random
		Xᵢₙᵢₜₘ = add_noiseSₘ(Xₜ,xyz_noise_angle,xyz_dis_noise[1],true);
		Yᵢₙᵢₜₘ = add_noiseSₘ(Yₜ,xyz_noise_angle,xyz_dis_noise[2],true);
		Zᵢₙᵢₜₘ = add_noiseSₘ(Zₜ,xyz_noise_angle,xyz_dis_noise[3],true);
		Initₘ = [Xᵢₙᵢₜₘ,Yᵢₙᵢₜₘ,Zᵢₙᵢₜₘ];
		Xᵢₙᵢₜ₂ = motorS2tfm(Xᵢₙᵢₜₘ); Yᵢₙᵢₜ₂ = motorS2tfm(Yᵢₙᵢₜₘ); Zᵢₙᵢₜ₂ = motorS2tfm(Zᵢₙᵢₜₘ);
		Initₘ = motorS2rotorG3S.(Initₘ)
		Init = [Xᵢₙᵢₜ₂,Yᵢₙᵢₜ₂,Zᵢₙᵢₜ₂];
		RXᵢₙᵢₜ = Xᵢₙᵢₜ₂[SOneTo(3),SOneTo(3)]
		RYᵢₙᵢₜ = Yᵢₙᵢₜ₂[SOneTo(3),SOneTo(3)]
		RZᵢₙᵢₜ = Zᵢₙᵢₜ₂[SOneTo(3),SOneTo(3)]
	else
		# Xᵢₙᵢₜ,Yᵢₙᵢₜ,Zᵢₙᵢₜ = AXB_YCZ_Close(Aₙ₂,Bₙ₂,Cₙ₂,conf);
		Xᵢₙᵢₜₘ,Yᵢₙᵢₜₘ,Zᵢₙᵢₜₘ = get_XYZ_rotor_Close(Aₙₘₛ,Bₙₘₛ,Cₙₘₛ)
		Initₘ = [Xᵢₙᵢₜₘ,Yᵢₙᵢₜₘ,Zᵢₙᵢₜₘ];
		if shared
			Init = r2t.(rotorG32rm.(Initₘ))
			Init = T_Wang(Init[1],Init[2],Init[3])
			RXᵢₙᵢₜ,_ = T2rt(Init[1]);
			RYᵢₙᵢₜ,_ = T2rt(Init[2]);
			RZᵢₙᵢₜ,_ = T2rt(Init[3]);
		else
			Xᵢₙᵢₜ,Yᵢₙᵢₜ,Zᵢₙᵢₜ = AXB_YCZ_Close(Aₙ₂,Bₙ₂,Cₙ₂,conf);
			Init = [Xᵢₙᵢₜ,Yᵢₙᵢₜ,Zᵢₙᵢₜ];
		 	Init .= ProjectToSE3.(Init)
			Init = T_Wang(Init[1],Init[2],Init[3])
			RXᵢₙᵢₜ,RYᵢₙᵢₜ,RZᵢₙᵢₜ = liao_close(motorS2rotorG3S.(Aₙₘₛ[samples]),motorS2rotorG3S.(Bₙₘₛ[samples]),motorS2rotorG3S.(Cₙₘₛ[samples]))
		end
		# RXᵢₙᵢₜ,RYᵢₙᵢₜ,RZᵢₙᵢₜ = liao_close(Aₙₘₛ,Bₙₘₛ,Cₙₘₛ)
		# RXᵢₙᵢₜ = rotor2rm(motor2rotor(Xᵢₙᵢₜₘ));
		# RYᵢₙᵢₜ = rotor2rm(motor2rotor(Yᵢₙᵢₜₘ));
		# RZᵢₙᵢₜ = rotor2rm(motor2rotor(Zᵢₙᵢₜₘ));
	end
	md"Solve Init Value"
end

# ╔═╡ 94b298dc-8c34-44fa-ba3c-f0104ff14d22
begin
	Vₘ₃,errₘ₃ₜ= iter_solve_Katyusha_G3(match_signS.(motorS2rotorG3S.(Ansₘ),Initₘ),Pₘ,confₘ₃);
	X₅,Y₅,Z₅,err₅ₜ = liao_iter(RA,RB,RC,RXᵢₙᵢₜ,RYᵢₙᵢₜ,RZᵢₙᵢₜ,conf₅)
	Pm = [tommeter.(P₂[i]) for i ∈ eachindex(P₂)]
	Result₂,Err₂ₜ = iter_wang(tommeter.(Init),Pm,conf)
	Result₂ = post_process_Wang(Result₂);
	# errₘ₃ = errₘ₃[1:5:size(errₘ₃,1),:];
	# Err₂ = Err₂[1:5:size(Err₂,1),:];
	# err₅ = err₅[1:5:size(err₅,1),:];
	md"Solve Here"
end

# ╔═╡ 2b6c7140-acb0-4ab4-a2b9-ad2bfc46182b
begin
	# Err G3 Init
	Xᵣ,Yᵣ,Zᵣ = T_Wang(rotorG32rm(Initₘ[1]), rotorG32rm(Initₘ[2]), rotorG32rm(Initₘ[3]));
	errₘ₃ = vcat(reshape(append!(get_error_XYZ([Xᵣ,Yᵣ,Zᵣ],[X,Y,Z]),[0,0,0]),1,9),errₘ₃ₜ)
	# errₘ₄ = vcat(reshape(append!(get_error_XYZ([Xᵣ,Yᵣ,Zᵣ],[X,Y,Z]),[0,0,0]),1,9),errₘ₄ₜ)
	# Err Wang Init
	Err₂ = vcat(reshape(append!(get_error_XYZ(Init,[X,Y,Z]),[0,0,0]),1,9),Err₂ₜ)
	# Err Liao Init
	Xᵣ,Yᵣ,Zᵣ = T_Wang(RXᵢₙᵢₜ,RYᵢₙᵢₜ,RZᵢₙᵢₜ);
	err₅ = vcat(reshape(append!(get_error_XYZ([Xᵣ,Yᵣ,Zᵣ],[X,Y,Z]),[0,0,0]),1,9),err₅ₜ)
	md"Append Init Value Error"
end

# ╔═╡ 8b927af6-1001-4006-b59c-0495f1ba2e33
function one_plot(i,r,title,labels,scale=1e3;ylabel=false,zoom=true)
	max = maximum([maximum(Err₂[:,i]),maximum(errₘ₃[:,i]),maximum(err₅[:,i])])
	min = minimum([minimum(Err₂[:,i]),minimum(errₘ₃[:,i]),minimum(err₅[:,i])])
	step = 1;
	P = plot(1:step:size(Err₂,1),Err₂[1:step:end,i]*scale,
		label=labels[1],
		linecolor=RGB(109/255, 174/255, 209/255),
		linestyle=:dash,
		lw=2,
		# markershape=:utriangle,
		# markercolor=RGB(109/255, 174/255, 209/255),
		# markersize=5,
		guidefont = (14,"times"), tickfont = (12), legendfontsize = 12, grid=true)
	plot!(1:step:size(err₅,1),err₅[1:step:end,i]*scale,
		label=labels[3],
		linecolor=RGB(17/255, 70/255, 127/255),
		linestyle=:dashdotdot,
		lw=2,
		# markershape=:star5,
		# markercolor=RGB(17/255, 70/255, 127/255),
		# markersize=5
	)
	plot!(1:step:size(errₘ₃,1),errₘ₃[1:step:end,i]*scale,
		label=labels[2],
		linecolor=RGB(182/255, 35/255, 48/255),
		linestyle=:solid,
		linewidth=2,
		# markershape=:circle,
		# markercolor=RGB(182/255, 35/255, 48/255),
		# markersize=3
	)
	if isa(ylabel,String)
		plot!(ylabel=ylabel)
	end
	# plot!(1:size(errₘ₄,1),errₘ₄[:,i]*scale,label=labels[2],linecolor=:red,linestyle=:dash,lw=2)
	# plot!(ylim=[min-(max-min)*0.05,min+(max-min)*1.1]*scale)
	# plot!(0:size(errₘ₄,1),vcat(errᵢₙᵢₜ[i],errₘ₄[:,i]),label=labels[4])
	# plot!(title=title)
	# plot!(ylim=[0,max*scale*1.01])
	if !random & (N_iter>10) & zoom
		max = maximum([maximum(Err₂[end-r:end,i]),maximum(errₘ₃[end-r:end,i]),maximum(err₅[end-r:end,i])])
		min = minimum([minimum(Err₂[end-r:end,i]),minimum(errₘ₃[end-r:end,i]),minimum(err₅[end-r:end,i])])
		
		plot!(size(Err₂,1)-r:size(Err₂,1),Err₂[end-r:end,i]*scale,
			label=false,
			inset = (1, bbox(0.05, 0.4, 0.6, 0.25, :right)),
			subplot=2,
			linecolor=RGB(109/255, 174/255, 209/255),
			linestyle=:dash,
			# markershape=:utriangle,
			# markercolor=RGB(109/255, 174/255, 209/255),
			# markersize=5
		)
		
		plot!(size(err₅,1)-r:size(err₅,1),err₅[end-r:end,i]*scale,
			label=false,
			subplot=2,
			linecolor=RGB(17/255, 70/255, 127/255),
			linestyle=:dashdotdot,
			# markershape=:star5,
			# markercolor=RGB(17/255, 70/255, 127/255),
			# markersize=5
		)
		
		plot!(size(errₘ₃,1)-r:size(errₘ₃,1),errₘ₃[end-r:end,i]*scale,
			label=false,
			subplot=2,
			linecolor=RGB(182/255, 35/255, 48/255),
			linestyle=:solid,
			# markershape=:circle,
			# markercolor=RGB(182/255, 35/255, 48/255),
			# markersize=3
		)
		
		# plot!(size(errₘ₄,1)-r:size(errₘ₄,1),errₘ₄[end-r:end,i]*scale,label=false,subplot=2)
		
		xticklabels = 
			[mod(x,5)==0 ? @sprintf("%d",x) : " " for x in N_iter-r+1:N_iter ]
		
		plot!(
			xticks=(collect(N_iter-r+1:N_iter),xticklabels),
			xtickfont=font(8),subplot=2)
		
		plot!(ylim=[min-(max-min)*0.15,min+(max-min)*1.15]*scale,subplot=2)
		
		plot!(yticks = false,subplot=2)
		# if scale != 1
		# 	annotate!(N_iter-r-1,(((max-min)*1.1+min)*scale),("×10⁻³",:top,:left,7),subplot=2)
		# end
	end
	return P
end

# ╔═╡ b03f82f4-8901-4ae2-897c-8a28b70e3990
begin
	l = @layout [a b c]
	P1 = one_plot(1,r,"Angle Error X",[false,false,false,false];ylabel="X Rot Err / 10^{-3} rad",zoom=false)
	# plot!(xlabel="#Iter",ylabel="X Angle Error/rad")
	P2 = one_plot(3,r,"Angle Error Y",[false,false,false,false];ylabel="Y Rot Err / 10^{-3} rad",zoom=false)
	# plot!(xlabel="#Iter",ylabel="Y Angle Error/rad")
	P3 = one_plot(5,r,"Angle Error Z",[false,false,false,false];ylabel="Z Rot Err / 10^{-3} rad",zoom=false)
	# plot!(xlabel="#Iter",ylabel="Z Angle Error/rad")
	plot!(legend=:topright)
	plot(P1,P2,P3,layout = l,size=(800,250))
end

# ╔═╡ 83ee5a08-751b-4afc-be0f-320bd33e8a26
begin
	P4 = one_plot(2,r,"Dis Error X",[false,false,false,false],1e3;ylabel="X Dis Err / mm",zoom=false)
	P5 = one_plot(4,r,"Dis Error Y",[false,false,false,false],1e3;ylabel="Y Dis Err / mm",zoom=false)
	P6 = one_plot(6,r,"Dis Error Z",[false,false,false,false],1e3;ylabel="Z Dis Err / mm",zoom=false)
	plot!(legend=:bottomright)
	plot(P4,P5,P6,layout = l,size=(800,250))
end

# ╔═╡ 07d07297-7426-4f7e-a3ec-0551cb669da6
function get_VPG3(n)
	Aₘₛ,Bₘₛ,Cₘₛ = get_ABCₘ(n,Xₜ,Yₜ,Zₜ,unit);
	Aₙₘₛ = add_noiseₘ.(Aₘₛ,angle_noise,dis_noise);
	Bₙₘₛ = add_noiseₘ.(Bₘₛ,angle_noise,dis_noise);
	Cₙₘₛ = add_noiseₘ.(Cₘₛ,angle_noise,dis_noise);
	Xᵢₙᵢₜₘ = add_noiseₘ(Xₜ,xyz_noise_angle,xyz_dis_noise[1]);
	Yᵢₙᵢₜₘ = add_noiseₘ(Yₜ,xyz_noise_angle,xyz_dis_noise[2]);
	Zᵢₙᵢₜₘ = add_noiseₘ(Zₜ,xyz_noise_angle,xyz_dis_noise[3]);
	Initₘ = [Xᵢₙᵢₜₘ,Yᵢₙᵢₜₘ,Zᵢₙᵢₜₘ];
	V = scale_motor.(match_sign.(Ansₘ,Initₘ),scale);
	Pₘ = [Aₙₘₛ,Bₙₘₛ,Cₙₘₛ];
	V = toVG3(motor2rvectorS.(V))
	Pₘ = [motor2rvectorS.(Pₘ[i]) for i∈eachindex(Pₘ)]
	return V,Pₘ	
end

# ╔═╡ 3fcada44-8c14-417d-b028-1fb969fbe1f5
begin
	conf₂ = ((max_iter = N_iter₂, μ = [1,1e-5,1,1,1],m=3, η=0.1,stop = false,err="Ext",errfunc = err_Wang_Full,errdim=9,τ=0.3,reg = true,svd=false))
	confₘ₃₂ = ((max_iter = N_iter₂, m = 3, η=0.3,μ=0.1,scale = scale, stop = false, err="Ext",errfunc = ErrG3_Full,errdim=9,τ=0.3,reg = true))
	conf₅₂ = ((max_iter = N_iter₂,err="Ext",errfunc = Err_Liao,errdim=9))
	md"Config Here"
end

# ╔═╡ 5ebf840c-8af1-4c36-b492-108440a59a1b
begin
	Vₘ₃₀,errₘ₃₀= iter_solve_Katyusha_G3(match_signS.(motorS2rotorG3S.(Ansₘ),Initₘ₂),Pₘ,confₘ₃₂);
	X₅₀,Y₅₀,Z₅₀,err₅₀ = liao_iter(RA,RB,RC,RXᵢₙᵢₜ₀,RYᵢₙᵢₜ₀,RZᵢₙᵢₜ₀,conf₅₂)
	Result₂₀,Err₂₀ = iter_wang(Init₂,P₂,conf₂)
	Result₂₀ = post_process_Wang(Result₂₀);
	md"Solve Here"
end

# ╔═╡ f14e175b-eb1f-4efa-a678-429af57431de
function one_plot₂(i,r,title,labels,scale=1e3)
	P = plot(0:size(Err₂₀,1),vcat(errᵢₙᵢₜ₀[i],Err₂₀[:,i]),
		label=labels[1],
		linecolor=RGB(109/255, 174/255, 209/255),
		linestyle=:dash,
		lw=2,
		guidefont = (14,"times"), tickfont = (12), legendfontsize = 12, grid=true)
	# plot!(0:size(errₘ₄,1),vcat(errᵢₙᵢₜ[i],errₘ₄[:,i]),label=labels[3])
	plot!(0:size(err₅₀,1),vcat(errᵢₙᵢₜ₀[i],err₅₀[:,i]),
		label=labels[3],
		lw = 2,
		linestyle=:dashdotdot,
		linecolor=RGB(17/255, 70/255, 127/255))
	plot!(0:size(errₘ₃₀,1),vcat(errᵢₙᵢₜ₀[i],errₘ₃₀[:,i]),
		label=labels[2],
		lw = 2,
		linecolor=RGB(182/255, 35/255, 48/255))
	plot!(title=title)
	if (N_iter>100)
		max = maximum([maximum(Err₂₀[end-r:end,i]),maximum(errₘ₃₀[end-r:end,i]),maximum(err₅₀[end-r:end,i])])
		min = minimum([minimum(Err₂₀[end-r:end,i]),minimum(errₘ₃₀[end-r:end,i]),minimum(err₅₀[end-r:end,i])])
		plot!(size(Err₂₀,1)-r:size(Err₂₀,1),Err₂[end-r:end,i]*scale,label=false,inset = (1, bbox(0.05, 0.05, 0.7, 0.3, :top, :right)),subplot=2,linecolor=:red)
		plot!(size(errₘ₃₀,1)-r:size(errₘ₃₀,1),errₘ₃₀[end-r:end,i]*scale,label=false,subplot=2,linecolor=:green)
		# plot!(size(errₘ₄,1)-r:size(errₘ₄,1),errₘ₄[end-r:end,i]*scale,label=false,subplot=2)
		plot!(size(err₅₀,1)-r:size(err₅₀,1),err₅₀[end-r:end,i]*scale,label=false,subplot=2,linecolor=:blue)
		ticklabels = [mod(x,5)==0 ? @sprintf("%d",x) : " " for x in N_iter-r:N_iter ]
		plot!(xticks=(collect(N_iter-9:N_iter),ticklabels),xtickfont=font(8),subplot=2)
		annotate!(N_iter-r-1,max,(@sprintf("%.3f",max),:top,:left,7),subplot=2)
	end
	return P
end

# ╔═╡ 502fb1ee-7384-4d1e-8cb7-cf24d737d31b
begin
	P10 = one_plot₂(1,r,"Angle Error X",[false,false,false],1)
	P20 = one_plot₂(3,r,"Angle Error Y",[false,false,false],1)
	P30 = one_plot₂(5,r,"Angle Error Z",["Wang","G3","Wu"],1)
	plot!(legend=:right)
	plot(P10,P20,P30,layout = l)
end

# ╔═╡ 19afc054-ec13-40ea-8683-15abc7c9b5c6
begin
	P40_ = one_plot₂(2,r,"Transl Error X",[false,false,false],1)
	P50_ = one_plot₂(4,r,"Transl Error Y",[false,false,false],1)
	P60_ = one_plot₂(6,r,"Transl Error Z",["Wang","G3","Wu"],1)
	plot!(legend=:right)
	plot(P40_,P50_,P60_,layout = l)
end

# ╔═╡ Cell order:
# ╟─209b3ee9-69a0-4f6b-aefe-4af20d2e583d
# ╟─c7f66c54-6e8d-4dcf-8e9a-9996f230309d
# ╟─dc0e522e-f3dc-4bf6-96f5-d402f09fd4e6
# ╟─1c7f849b-4f1c-421f-bf34-ce2dcd46e955
# ╟─94499cf7-66b3-4920-af7b-0bfeaa3f658b
# ╟─c338bab3-7d44-400c-8aea-e5e7820d9ba6
# ╟─bf930527-f561-4537-904b-c2cba3474b54
# ╟─77683183-1f55-4436-ab26-934e452625c3
# ╟─61d732b0-19f3-4d34-8fd3-ab224300f828
# ╟─f492eecf-b9a1-47a5-873c-1b6f23088664
# ╟─95692758-52e4-4d84-818d-53abc3209fe0
# ╟─b03f82f4-8901-4ae2-897c-8a28b70e3990
# ╟─ca5f77b1-9b44-4496-9f53-8dd04467f934
# ╟─83ee5a08-751b-4afc-be0f-320bd33e8a26
# ╟─56079b1e-a80b-4818-8fbd-b2c0138185c6
# ╟─3f2be549-1fe8-45c8-9b23-c33aadb3ed93
# ╟─71792be6-7a13-4f8d-aead-c327f1b4fc36
# ╟─0a661256-be96-4359-9515-eff2cc0ebbd5
# ╟─94b298dc-8c34-44fa-ba3c-f0104ff14d22
# ╟─2b6c7140-acb0-4ab4-a2b9-ad2bfc46182b
# ╟─26c72e43-07bb-4381-b96e-bec3519a0469
# ╟─b1ca5ec0-d291-49ed-b898-43da84854477
# ╟─060f90e0-9d9d-4557-9dfb-5ce6aecea51a
# ╟─591c9131-f99e-421a-be3c-92e0b379a2f9
# ╟─b7e4e47b-2c43-4e5c-b74a-c10c6062225c
# ╟─ca204e58-57a9-4a96-acd5-c08497e13317
# ╟─502fb1ee-7384-4d1e-8cb7-cf24d737d31b
# ╟─19afc054-ec13-40ea-8683-15abc7c9b5c6
# ╟─b7392bcf-d3cc-4838-879f-657d385c9037
# ╟─3fcada44-8c14-417d-b028-1fb969fbe1f5
# ╟─5ebf840c-8af1-4c36-b492-108440a59a1b
# ╠═ba20da2e-c3c5-416e-ace2-3d2555aebbcb
# ╟─43ec9469-619e-4467-876e-e49f3f8955e5
# ╟─8b927af6-1001-4006-b59c-0495f1ba2e33
# ╟─f14e175b-eb1f-4efa-a678-429af57431de
# ╟─6a8b84bf-1de4-4f7b-99e8-4ca87021202f
# ╟─07d07297-7426-4f7e-a3ec-0551cb669da6
# ╟─20fccdd0-1f3d-421a-8fe2-a13176df32c2
# ╟─d7f4c0f3-296e-49cd-b453-bf0778a1608f
# ╟─ba87cb0a-697a-462c-941a-ab9605d35a16
# ╟─ad34b9bf-7162-4682-aa71-7ec2aadefd27
# ╟─ebc58e1b-fae3-4e1e-9451-b150417b76ca
# ╟─3269f5e3-c216-4a27-a6ea-4e9cc297dced
