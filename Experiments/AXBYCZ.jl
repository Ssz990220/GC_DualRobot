### A Pluto.jl notebook ###
# v0.19.27

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

# ╔═╡ bfc7e3b2-6509-455e-9f0e-a6b6e01ffe27
begin
	using Pkg
	Pkg.activate("../.")
end

# ╔═╡ 43daea20-166d-11ed-1b81-e3a8fde6e0fe
begin
	using MAT
	using StaticArrays
	using LinearAlgebra
	using Distributions
	using Random
	using Statistics
	using Plots
	using Rotations
	using Printf
	using StatsPlots
	using JLD
	using HypothesisTests
	using Reexport
end

# ╔═╡ d8a31e25-427c-43b5-83fb-7a4b4fbcf2b1
using PlutoUI

# ╔═╡ 0c8f94cd-331f-4f31-b641-702c31d55800
begin
	include("../project_utils/DualRobotCalibrator.jl")
	@reexport using .DualRobotCalibrator
end

# ╔═╡ 18ecf252-8608-43ae-8ea7-517303a9ad8d
begin
	include("../utils/ABB_4600_60_205.jl")
	include("../utils/JAKA_Zu12.jl")
	md"Import robots"
end

# ╔═╡ 3330f278-e9b7-4ecd-9a27-7376e8d170a4
md"# Load Aₛ Bₛ & Cₛ"

# ╔═╡ e9025d63-49f0-46be-9ade-796ddae8b9c6
md"**SET UNIT FIRST**, bad choice of UNIT may lead to divergence due to instability of some algorithms!!"

# ╔═╡ 566de0a8-2aac-4fef-ba60-d42d82e3f392
UNIT = "m";

# ╔═╡ e3774fab-c07d-41ab-8361-bea25c0f92ce
dataNames = ["EXPN_240"]

# ╔═╡ 664e3a06-38d7-4ec3-9e83-dfd2f9d02d13
MaDataPath3 = ["Ma/Ma_Fixedeye_40_6","Ma/Ma_Handeye_40_6","Ma/Ma_Fixedeye_40_7","Ma/Ma_Handeye_40_7","Ma/Ma_Fixedeye_40_8","Ma/Ma_Handeye_40_8"]

# ╔═╡ e87f8f36-594a-4772-8fd7-eb174a949a07
begin
	_,_,Qs_ABBs,Qs_JAKAs,Ass,Bss,Css = load_data(dataNames,UNIT);
	DataAXXB = load_data(["AXXB240R"],UNIT); 	# Load Handeye Calibration data for Fu
	NN = 3; nn = 80;
	DataMa = load_data(MaDataPath3,UNIT);
end

# ╔═╡ d2acbbbf-73ab-41f6-b279-46a2d08a4c77
# Match Sign for Quaternion/GA/G3 Based Algorithms
Aₛ,Bₛ,Cₛ,x,y,z = match_sign_ABC(Ass,Bss,Css,Qs_JAKAs,Qs_ABBs,JAKA,ABB4600,"m");

# ╔═╡ cca2762b-4edb-474e-bcbb-39efd7cfd680
md"# Config"

# ╔═╡ f86679b5-86a9-49c7-8842-bf6eeea18825
N_iter = 200;

# ╔═╡ 066168f4-a9ca-4ab6-8f6f-5bf98b7b7a37
md"Solver Configure defined in the next cell"

# ╔═╡ 4ed5d4ad-0815-452d-8bdf-9c1c567be310
begin
	conf = ((max_iter = N_iter, μ = [1.0,1e-5,1.0,1.0,1.0],m=3, η=0.1, err = "NAN", stop = false, reg=true,svd=true,τ=0.3)); 	# Wang
	confₗ = ((max_iter = N_iter,))				# Wu
	confₘ = ((max_iter = N_iter, m = 3, η=0.1,μ=0.1,scale=1e-4, err = "NAN", stop = false, reg=true))			# GA-SVRG
	confₘ₂ = ((max_iter = N_iter, m = 3, η=0.1,μ=0.1,scale=1e-5, err = "NAN", stop = false, reg=true,τ=0.2))		# GA-Katyusha
	confₘ₃ = ((max_iter = N_iter, m = 3, η=0.3,μ=0.1,scale=1e-2, err = "NAN", stop = false, reg=true,τ=0.2))		# G3
end

# ╔═╡ c08040e8-df14-4cb2-af0c-207ae85b1c52
md"# Four Fold Cross Validation
According to the rules proposed by Ma and in this article, two datasets of size 240 are generated. Like in Ma's paper, Ma's dataset consists of 3 pairs of sub-datasets, each pair containing two datasets of size 40. In each pair of datasets, half of the data fixes the JAKA robot and moves the ABB robot in a small range, while the other half fixes the ABB robot and moves the JAKA robot in a small range (following Ma's paper). In addition, 60 additional data points are collected for Fu's method to correct the AX=XB problem.

Four-fold cross-validation is performed using the above data. However, because Ma requires the dataset to be \"as compact as possible,\" the validation set data is too close to the training set data and is not suitable for cross-validation. Therefore, in each cross-validation, the part of the dataset collected using the method proposed in this article is used for \"cross-validation\" for testing purposes.
"

# ╔═╡ b05de384-cd1f-4b43-89c9-c7d501582f50
md"
```julia
for i = 1:200
	splitData
	runTest
	CrossCheck
	collectResult
end
```
"

# ╔═╡ f1015ab1-23b7-43ce-bc0f-60ce068bc615
N_CrossCheck = 200

# ╔═╡ f3878dfc-2b49-42d7-91ed-14f77319033e
md"The following cell conducting 4-fold cross check takes ~4min on my workstation"

# ╔═╡ ae16cffc-e369-473f-b2b5-9bf1d4ff4d81
md"Toggle the following checkbox to display result of PGA"

# ╔═╡ 55e8a4dc-5ecf-47fb-9751-5a7e3882d16c
@bind displayGA CheckBox()

# ╔═╡ 78938e52-22c1-4910-91dc-8616c9a51240
thresholdGA = 5e-3 ## Filter out failed samples

# ╔═╡ c1bd5255-02ee-449e-a1e4-e235533a7994
md"## Solve Once
All data collected are used to solve the problem"

# ╔═╡ 42ddfa21-d01f-4102-8b7e-c041207378d4
NDATA = 240;

# ╔═╡ f4761271-2719-474a-b150-4eaa645e537b
md"## Save Result"

# ╔═╡ c6324192-1a16-4c26-8476-b00ff52acf2e
function saveResult(WANG,G3,WU,MAI,MAH,FU,PGA,Name)

	X,Y,Z = WANG
	Xg3,Yg3,Zg3 = G3
	Xₗ,Yₗ,Zₗ = WU
	XMa,YMa,ZMa = MAI
	XMa2,YMa2,ZMa2 = MAH
	XFu,YFu,ZFu = FU
	XGA,YGA,ZGA = PGA
	
if UNIT == "mm"
	matwrite("./result/XYZ/"*Name*".mat", Dict(
		"Xm" => Matrix(X),
		"Ym" => Matrix(Y),
		"Zm" => Matrix(Z),
		"Xg3" => Matrix(Xg3),
		"Yg3" => Matrix(Yg3),
		"Zg3" => Matrix(Zg3),
		"Xl" => Matrix(Xₗ),
		"Yl" => Matrix(Yₗ),
		"Zl" => Matrix(Zₗ),
		"XFu" => Matrix(XFu),
		"YFu" => Matrix(YFu),
		"ZFu" => Matrix(ZFu),
		"XMa" => Matrix(XMa),
		"YMa" => Matrix(YMa),
		"ZMa" => Matrix(ZMa),
		"XGA" => Matrix(XGA),
		"YGA" => Matrix(YGA),
		"ZGA" => Matrix(ZGA)
	); compress = true)
else
	matwrite("./result/XYZ/"*Name*".mat", Dict(
		"Xm" => Matrix(tommeter(X)),
		"Ym" => Matrix(tommeter(Y)),
		"Zm" => Matrix(tommeter(Z)),
		"Xg3" => Matrix(tommeter(Xg3)),
		"Yg3" => Matrix(tommeter(Yg3)),
		"Zg3" => Matrix(tommeter(Zg3)),
		"Xl" => Matrix(tommeter(Xₗ)),
		"Yl" => Matrix(tommeter(Yₗ)),
		"Zl" => Matrix(tommeter(Zₗ)),
		"XFu" => Matrix(tommeter(XFu)),
		"YFu" => Matrix(tommeter(YFu)),
		"ZFu" => Matrix(tommeter(ZFu)),
		"XMa" => Matrix(tommeter(XMa)),
		"YMa" => Matrix(tommeter(YMa)),
		"ZMa" => Matrix(tommeter(ZMa)),
		"XGA" => Matrix(tommeter(XGA)),
		"YGA" => Matrix(tommeter(YGA)),
		"ZGA" => Matrix(tommeter(ZGA))
	); compress = true)
	md"Results are saved at ./result/XYZ/$Name.mat"
end
end

# ╔═╡ 31cda4bd-4a74-40b9-8199-81d0bde9f614
md"### You can run \"validation\_set\_processing.m\" with MATLAB to view the final result for measurement experiment"

# ╔═╡ ef6a15fc-7a25-4084-8360-fe3c35f0f87b
md"**The followings are temporary functions in this script**

## Cross Validation Functions"

# ╔═╡ 4ad28737-4d10-4240-9810-e405cfc2aade
function setBuffer(N_CrossCheck,NFold,sizeTestSet)
	dimBuffer = N_CrossCheck*NFold*sizeTestSet
	eWs = zeros(dimBuffer);eWts = zeros(dimBuffer);
	eWus = zeros(dimBuffer);eWuts = zeros(dimBuffer);
	eFus = zeros(dimBuffer);eFuts = zeros(dimBuffer);
	eG3s = zeros(dimBuffer);eG3ts = zeros(dimBuffer);
	eMas = zeros(dimBuffer);eMats = zeros(dimBuffer);
	eGAs = zeros(dimBuffer);eGAts = zeros(dimBuffer);
	return eWs, eWts, eWus, eWuts, eFus, eFuts, eG3s, eG3ts, eMas, eMats, eGAs, eGAts
end

# ╔═╡ 0162dd87-f953-4f2f-90c7-d9d389543a73
function extractDataByIndex(Data,idx,NFold,i)
	NData = size(Data[1],1); BatchSize = NData÷NFold;
	@assert mod(NData,NFold) == 0;
	A = Data[1]; B = Data[2]; C = Data[3]; Aₘ = Data[4]; Bₘ = Data[5]; Cₘ = Data[6];
	As = A[idx[(i-1)*BatchSize+1:i*BatchSize]]
	Bs = B[idx[(i-1)*BatchSize+1:i*BatchSize]]
	Cs = C[idx[(i-1)*BatchSize+1:i*BatchSize]]
	Asₘ = Aₘ[idx[(i-1)*BatchSize+1:i*BatchSize]]
	Bsₘ = Bₘ[idx[(i-1)*BatchSize+1:i*BatchSize]]
	Csₘ = Cₘ[idx[(i-1)*BatchSize+1:i*BatchSize]]
	return (As,Bs,Cs,Asₘ,Bsₘ,Csₘ)
end

# ╔═╡ 5116e753-cd4c-4064-ac6f-676e0327ff4e
function extractDataByIndexAXXB(Data,idx,NFold,i)
	NData = size(Data[1],1); BatchSize = NData÷NFold;
	@assert mod(NData,NFold) == 0;
	A = Data[1]; B = Data[2];
	As = A[idx[(i-1)*BatchSize+1:i*BatchSize]]
	Bs = B[idx[(i-1)*BatchSize+1:i*BatchSize]]
	return (As,Bs)
end

# ╔═╡ e4de16cd-8345-4191-b383-7fb92dd18247
function extractDataByIndexMa(Data,idx,NMa,NFold,j)
	Ntotal = size(DataMa[1],1); nn  = Ntotal ÷ NMa ÷ 2;
	BatchSize = nn ÷ NFold;
	@assert mod(nn,NFold) == 0
	A = [Data[1][(i-1)*nn+1:i*nn] for i = 1:NMa*2]
	B = [Data[2][(i-1)*nn+1:i*nn] for i = 1:NMa*2]
	C = [Data[3][(i-1)*nn+1:i*nn] for i = 1:NMa*2]
	ASelected = [A[i][idx[i][(j-1)*BatchSize+1:j*BatchSize]] for i = 1:NMa*2]
	BSelected = [B[i][idx[i][(j-1)*BatchSize+1:j*BatchSize]] for i = 1:NMa*2]
	CSelected = [C[i][idx[i][(j-1)*BatchSize+1:j*BatchSize]] for i = 1:NMa*2]
	return (ASelected,BSelected,CSelected)
end

# ╔═╡ 23e2848a-e6ba-44c7-afa3-5cdccd4a2522
function splitData(Data::T,DataMa,DataAXXB,NMa,NFold) where {T<:Vector{Vector}}
	NData = size(Data[1],1);
	s = sample(1:NData,NData,replace=false)
	sMa = [sample(1:NData÷NMa÷2,NData÷NMa÷2,replace=false) for i = 1:NMa*2]
	Datasets = [extractDataByIndex(Data,s,NFold,i) for i = 1:NFold]
	DatasetsAXXB = [extractDataByIndexAXXB(DataAXXB,s,NFold,i) for i = 1:NFold]
	DatasetsMa = [extractDataByIndexMa(DataMa,sMa,NMa,NFold,j) for j = 1:NFold]
	return Datasets, DatasetsMa, DatasetsAXXB
end

# ╔═╡ fa417aca-2aa1-4396-9292-efd7adb77c09
function MergeData(dataset, idx, NFold)
	Nidx = 1:NFold
	Nidx = Nidx[Nidx.!=idx]
	DataA = Vector{SMatrix{4,4,Float64,16}}(undef,0)
	DataB = Vector{SMatrix{4,4,Float64,16}}(undef,0)
	DataC = Vector{SMatrix{4,4,Float64,16}}(undef,0)
	DataAₘ = Vector{SVector{8,Float64}}(undef,0)
	DataBₘ = Vector{SVector{8,Float64}}(undef,0)
	DataCₘ = Vector{SVector{8,Float64}}(undef,0)
	for id ∈ Nidx
		DataA = vcat(DataA,dataset[id][1])
		DataB = vcat(DataB,dataset[id][2])
		DataC = vcat(DataC,dataset[id][3])
		DataAₘ = vcat(DataAₘ,dataset[id][4])
		DataBₘ = vcat(DataBₘ,dataset[id][5])
		DataCₘ = vcat(DataCₘ,dataset[id][6])
	end
	DataTrain = [DataA,DataB,DataC,DataAₘ,DataBₘ,DataCₘ]
	DataTest = dataset[idx]
	return DataTrain, DataTest
end

# ╔═╡ f28a1168-5c38-4b9a-b4be-8adf06ee8fcb
function MergeDataAXXB(dataset, idx, NFold)
	Nidx = 1:NFold
	Nidx = Nidx[Nidx.!=idx]
	DataA = Vector{SMatrix{4,4,Float64,16}}(undef,0)
	DataB = Vector{SMatrix{4,4,Float64,16}}(undef,0)
	for id ∈ Nidx
		DataA = vcat(DataA,dataset[id][1])
		DataB = vcat(DataB,dataset[id][2])
	end
	DataTrain = [DataA,DataB]
	return DataTrain
end

# ╔═╡ 79ae6ae1-e0d8-4e6e-9dbe-d45da8e291e6
function MergeDataMa(dataset,idx,NFold,NMa)
	Nidx = 1:NFold; Nidx = Nidx[Nidx.!=idx]
	DataA = Vector{SMatrix{4,4,Float64,16}}(undef,0)
	DataB = Vector{SMatrix{4,4,Float64,16}}(undef,0)
	DataC = Vector{SMatrix{4,4,Float64,16}}(undef,0)
	for poseId = 1:NMa*2
		for id ∈ Nidx
			DataA = vcat(DataA,dataset[id][1][poseId])
			DataB = vcat(DataB,dataset[id][2][poseId])
			DataC = vcat(DataC,dataset[id][3][poseId])
		end
	end
	DataTrainMa = [DataA,DataB,DataC]
	return DataTrainMa
end

# ╔═╡ 45eb7b65-16db-4940-a0e0-98520c84e4c1
md"## Functions"

# ╔═╡ 8e3a0481-d79e-4f32-9c80-bdc223c7edcd
function selectData(Data,n)
	return Data[1:n]
end

# ╔═╡ 216eca4f-0181-40cf-b3d7-78cc28a2d92f
function selectDataMa(Data,nsn,ns,NN)
	counter = 0
	dataS = Vector{SMatrix{4,4,Float64,16}}(undef,0)
	for i = 1:(NN*2)
		dataLocal = Data[counter .+ (1:nsn>>1)]
		dataS = vcat(dataS,dataLocal)
		counter = counter + ns>>1
	end
	return dataS
end

# ╔═╡ 2fbdb4e0-a962-47cd-a8f1-1d85f258d462
function solveAll(As,Bs,Cs,Aₛ,Bₛ,Cₛ,NN,ns,As_Ma,Bs_Ma,Cs_Ma,AXXBAs,AXXBBs)

	## GA
	
	Xga, Yga, Zga = GAClose(Aₛ,Bₛ,Cₛ)
	GA,_ = iter_solve_GC_katyushaᵥ([Xga, Yga, Zga],[Aₛ,Bₛ,Cₛ],confₘ₂)
	GA = motorS2tfm.(GA);
	
	## G3
	Xg3ᵣ,Yg3ᵣ,Zg3ᵣ = get_XYZ_rotor_Close(Aₛ,Bₛ,Cₛ)
	Result₃,err₃ = AXBYCZ_G3([Xg3ᵣ,Yg3ᵣ,Zg3ᵣ],[Aₛ,Bₛ,Cₛ],confₘ₃)
	Xg3 = Result₃[1]; Yg3 = Result₃[2]; Zg3 = Result₃[3];

	G3 = [Xg3,Yg3,Zg3];

	## Fu + Bayro
	XFu = AXXBs(tfm2motorS.(AXXBAs),tfm2motorS.(AXXBBs))
	YFu,ZFu = AXB_YCZ_FUS(Aₛ,Bₛ,Cₛ,XFu)
	XFu = motorS2tfm(XFu); YFu = motorS2tfm(YFu); ZFu = motorS2tfm(ZFu)
	FU = [XFu,YFu,ZFu];

	## Ma Initial Value
	XMai,YMai,ZMai = axbyczProb1(As_Ma[1:ns>>1],Bs_Ma[1:ns>>1],Cs_Ma[1:ns>>1],
							As_Ma[ns>>1+1:ns],Bs_Ma[ns>>1+1:ns],Cs_Ma[ns>>1+1:ns])
	# @show XMai
	XMAI = [XMai,YMai,ZMai]
	#### Ma Iterative
	XMa,YMa,ZMa = axbyczProb3(NN,ns,As_Ma,Bs_Ma,Cs_Ma,XMai,YMai,ZMai)
	# XMAI = [XMa,YMa,ZMa]

	#### Ma Hybrid
	AMa_mean = [meanCov(As_Ma[(i-1)*(ns>>1)+1:i*(ns>>1)])[1] for i = 1:2*NN]
	BMa_mean = [meanCov(Bs_Ma[(i-1)*(ns>>1)+1:i*(ns>>1)])[1] for i = 1:2*NN]
	CMa_mean = [meanCov(Cs_Ma[(i-1)*(ns>>1)+1:i*(ns>>1)])[1] for i = 1:2*NN]
	XWM,YWM,ZWM = AXBYCZ_Wang2(AMa_mean,BMa_mean,CMa_mean,XMai,YMai,ZMai)
	XMAH = [XWM,YWM,ZWM];

	Asm = tommeter.(As); Bsm = tommeter.(Bs); Csm = tommeter.(Cs);
	
	XW,YW,ZW,E = AXB_YCZ_Wang(Asm,Bsm,Csm,conf); 				# Wang
	WANG = [XW,YW,ZW]; WANG = tometer.(WANG);
	Xₗ,Yₗ,Zₗ = AXBYCZ_Liao(As,Bs,Cs,Aₛ,Bₛ,Cₛ,confₗ); 		# Wu
	WU = [Xₗ,Yₗ,Zₗ];

	return G3,FU,XMAI,XMAH,WANG,WU,GA
end

# ╔═╡ 884f3c97-a20f-4110-b3d9-0f41f50a59d9
begin
	ns = NDATA÷NN;

	# Solve
	G3,FU,MAI,MAH,WANG,WU,GA = solveAll(Ass,Bss,Css,Aₛ,Bₛ,Cₛ,NN,ns,DataMa[5],DataMa[6],DataMa[7],DataAXXB[5],DataAXXB[6]);
	Name = "XYZ_ALL"
	saveResult(WANG,G3,WU,MAI,MAH,FU,GA,Name)
end

# ╔═╡ fbf0769a-ab4c-4ed2-864c-deb0755d0f11
function get_cross_error(DataTest,X,Y,Z)
	return get_cross_error(DataTest[1],DataTest[2],DataTest[3],X,Y,Z)
end

# ╔═╡ ee33537c-a2b5-4f07-a73d-de21d08298aa
function get_cross_error(As,Bs,Cs,X,Y,Z)
	n = size(As,1)
	errᵣ = zeros(n); errₜ = zeros(n);
	for i = 1:n
		# e = As[i]*X*Bs[i]/(Y*Cs[i]*Z);
		# eᵣ = get_rot_error(e[SOneTo(3),SOneTo(3)],eye(3));
		# eₜ = norm(e[1:3,4]);
		eᵣ,eₜ = get_error(As[i]*X*Bs[i],Y*Cs[i]*Z)
		errᵣ[i] = eᵣ;
		errₜ[i] = eₜ;
	end
	return sum(errᵣ)/n, sum(errₜ)/n, errᵣ,errₜ
end

# ╔═╡ 61b42980-f0d1-4687-9a76-2f64ea0fd78c
function FourFoldCrossCheck(Data,DataMa,DataAXXBAll,NMa,N_CrossCheck)
	NFold = 4; 
	NData = size(Data[1],1); nn = size(DataMa[1],1)÷NMa;
	ns = nn÷NFold*(NFold-1); sizeTestSet = NData÷NFold;

	eWs, eWts, eWus, eWuts, eFus, eFuts, eG3s, eG3ts, eMas, eMats, eGAs, eGAts = setBuffer(N_CrossCheck,NFold,sizeTestSet)
	for i = 1:N_CrossCheck
		DataSets, DataSetsMa,DataSetsAXXB = splitData(Data,DataMa,DataAXXBAll,NMa,NFold)
		for j = 1:NFold
			DataTrain, DataTest = MergeData(DataSets,j,NFold)
			DataTrainAXXB = MergeDataAXXB(DataSetsAXXB,j,NFold)
			DataTrainMa = MergeDataMa(DataSetsMa,j,NFold,NMa)
			G3,FU,MAI,_,WANG,WU,GA = solveAll(DataTrain[1],DataTrain[2],DataTrain[3],DataTrain[4],DataTrain[5],DataTrain[6],NMa,ns,DataTrainMa[1],DataTrainMa[2],DataTrainMa[3],DataTrainAXXB[1],DataTrainAXXB[2]);
			
			X,Y,Z = WANG
			Xg3,Yg3,Zg3 = G3
			Xₗ,Yₗ,Zₗ = WU
			XMa,YMa,ZMa = MAI
			XFu,YFu,ZFu = FU
			XGA,YGA,ZGA = GA

			s = (i-1)*NFold*sizeTestSet + (j-1)*sizeTestSet + 1;
			e = (i-1)*NFold*sizeTestSet + j * sizeTestSet;
			
			_,_,eWs[s:e],eWts[s:e] = get_cross_error(DataTest,X,Y,Z)
			_,_,eG3s[s:e],eG3ts[s:e] = get_cross_error(DataTest,Xg3,Yg3,Zg3)

			wuFailedCount = 0
			while sum(isnan.(Xₗ)) > 0 			
## Sometimes Wu may find NaN result, a recompute will typically solve the problem
				Xₗ,Yₗ,Zₗ = AXBYCZ_Liao(DataTrain[1],DataTrain[2],DataTrain[3],DataTrain[4],DataTrain[5],DataTrain[6],confₗ)
				wuFailedCount = wuFailedCount + 1
				if wuFailedCount > 5 
					save(string("../images/","PaperV2","/","WuFailedData",".jld"),"Data",DataTrain)
					error("Wu Failed to Converge")
				end
			end
			_,_,eWus[s:e],eWuts[s:e] = get_cross_error(DataTest,Xₗ,Yₗ,Zₗ)
			_,_,eFus[s:e],eFuts[s:e] = get_cross_error(DataTest,XFu,YFu,ZFu)
			_,_,eMas[s:e],eMats[s:e] = get_cross_error(DataTest,XMa,YMa,ZMa)
			_,_,eGAs[s:e],eGAts[s:e] = get_cross_error(DataTest,XGA,YGA,ZGA)
		end
	end
	return eWs, eWts, eWus, eWuts, eFus, eFuts, eG3s, eG3ts, eMas, eMats, eGAs, eGAts
end

# ╔═╡ 3fa29a27-3abb-4097-a539-296e5b0f3127
md"# Import Packages"

# ╔═╡ Cell order:
# ╟─3330f278-e9b7-4ecd-9a27-7376e8d170a4
# ╟─e9025d63-49f0-46be-9ade-796ddae8b9c6
# ╠═566de0a8-2aac-4fef-ba60-d42d82e3f392
# ╟─e3774fab-c07d-41ab-8361-bea25c0f92ce
# ╠═664e3a06-38d7-4ec3-9e83-dfd2f9d02d13
# ╠═e87f8f36-594a-4772-8fd7-eb174a949a07
# ╠═d2acbbbf-73ab-41f6-b279-46a2d08a4c77
# ╟─cca2762b-4edb-474e-bcbb-39efd7cfd680
# ╠═f86679b5-86a9-49c7-8842-bf6eeea18825
# ╟─066168f4-a9ca-4ab6-8f6f-5bf98b7b7a37
# ╠═4ed5d4ad-0815-452d-8bdf-9c1c567be310
# ╟─c08040e8-df14-4cb2-af0c-207ae85b1c52
# ╟─b05de384-cd1f-4b43-89c9-c7d501582f50
# ╠═f1015ab1-23b7-43ce-bc0f-60ce068bc615
# ╟─61b42980-f0d1-4687-9a76-2f64ea0fd78c
# ╟─f3878dfc-2b49-42d7-91ed-14f77319033e
# ╠═46c79cc2-9a9d-4665-a136-2113720e3a58
# ╠═fc258ba9-e2ad-45de-9bb1-ccc03907137e
# ╟─ae16cffc-e369-473f-b2b5-9bf1d4ff4d81
# ╟─55e8a4dc-5ecf-47fb-9751-5a7e3882d16c
# ╠═db42db8f-b52b-4853-a4f4-1dd478074397
# ╟─6068c8cd-9d80-488b-b339-7beb0788bd03
# ╟─078410c5-9e1e-4980-a363-dc591474708c
# ╠═78938e52-22c1-4910-91dc-8616c9a51240
# ╠═f23fdff3-0eb0-4361-a303-894a6ff3b0d1
# ╟─c1bd5255-02ee-449e-a1e4-e235533a7994
# ╠═42ddfa21-d01f-4102-8b7e-c041207378d4
# ╠═884f3c97-a20f-4110-b3d9-0f41f50a59d9
# ╟─f4761271-2719-474a-b150-4eaa645e537b
# ╠═c6324192-1a16-4c26-8476-b00ff52acf2e
# ╟─31cda4bd-4a74-40b9-8199-81d0bde9f614
# ╟─ef6a15fc-7a25-4084-8360-fe3c35f0f87b
# ╠═4ad28737-4d10-4240-9810-e405cfc2aade
# ╟─23e2848a-e6ba-44c7-afa3-5cdccd4a2522
# ╟─0162dd87-f953-4f2f-90c7-d9d389543a73
# ╟─5116e753-cd4c-4064-ac6f-676e0327ff4e
# ╟─e4de16cd-8345-4191-b383-7fb92dd18247
# ╟─fa417aca-2aa1-4396-9292-efd7adb77c09
# ╟─f28a1168-5c38-4b9a-b4be-8adf06ee8fcb
# ╟─79ae6ae1-e0d8-4e6e-9dbe-d45da8e291e6
# ╟─45eb7b65-16db-4940-a0e0-98520c84e4c1
# ╟─8e3a0481-d79e-4f32-9c80-bdc223c7edcd
# ╟─216eca4f-0181-40cf-b3d7-78cc28a2d92f
# ╟─2fbdb4e0-a962-47cd-a8f1-1d85f258d462
# ╟─fbf0769a-ab4c-4ed2-864c-deb0755d0f11
# ╟─ee33537c-a2b5-4f07-a73d-de21d08298aa
# ╟─3fa29a27-3abb-4097-a539-296e5b0f3127
# ╟─bfc7e3b2-6509-455e-9f0e-a6b6e01ffe27
# ╠═43daea20-166d-11ed-1b81-e3a8fde6e0fe
# ╠═d8a31e25-427c-43b5-83fb-7a4b4fbcf2b1
# ╟─18ecf252-8608-43ae-8ea7-517303a9ad8d
# ╟─0c8f94cd-331f-4f31-b641-702c31d55800
