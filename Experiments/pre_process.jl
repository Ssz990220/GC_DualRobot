using DelimitedFiles
# using Rotations
using MAT
# include("../utils/Robotics.jl")
# include("../project_utils/utils.jl")
# include("../project_utils/AXBYCZ_utils.jl")

"""
	read_Data(dir)
Load `A,B,C`s from directory `dir`
"""
function read_Data(dir)
    Data_ABB = readdlm(dir*"/ABBTraj.txt",',',Float64,'\n')
    n = size(Data_ABB,1);
    Qs_ABB = [([Data_ABB[i,j] for j = 1:6]./180*π) for i = 1:n]
    XYZ_ABB = [Data_ABB[i,7:9] for i = 1:n]
    Quat_ABB = [Data_ABB[i,10:end] for i = 1:n]

    R_ABB = Quat2rm.(Quat_ABB)
    T_ABB = rt2T.(R_ABB,XYZ_ABB)

    Data_JAKA = readdlm(dir*"/JAKATraj.txt",',',Float64,'\n')
    XYZ_JAKA = [Data_JAKA[i,1:3] for i = 1:n]
    RXYZ_JAKA = [(Data_JAKA[i,4:6]./180*π) for i = 1:n]
    Qs_JAKA = [(Data_JAKA[i,7:12]./180*π) for i = 1:n]

    R_JAKA = [RotZYX(RXYZ_JAKA[i][3],RXYZ_JAKA[i][2],RXYZ_JAKA[i][1]) for i = 1:n]
    T_JAKA = rt2T.(R_JAKA,XYZ_JAKA)
    return T_ABB,T_JAKA,Qs_ABB,Qs_JAKA
end

"""
	clean_up(As₁,Bs₁,Cs₁)
Use pre calibrated data (X,Y,Z) to filter bad data points.
"""
function clean_up(As₁,Bs₁,Cs₁,prefix = ".")
	file = matopen(prefix * "/result/XYZ/XYZ.mat");
	X = read(file,"Xga"); Y = read(file,"Yga"); Z = read(file,"Zga");
	X = SMatrix{4,4,Float64,16}(X);
	Y = SMatrix{4,4,Float64,16}(Y);
	Z = SMatrix{4,4,Float64,16}(Z);
	close(file);
	n = size(As₁,1);
	eᵣ = zeros(n); eₜ = zeros(n);
	good = Vector{Bool}(undef, n)
	for i = 1:n
		e = get_error(As₁[i]*X*Bs₁[i],Y*Cs₁[i]*Z)
		if (e[1] > 0.1) | (e[2] > 10)
			good[i] = false;
		else
			good[i] = true;
		end
	end
	
	if sum(good) != n
		println("Found ",n-sum(good)," bad data points")
	end
	
	# for i = 1:n
	# 	good[i] = true;
	# end
	good
end

"""
	load_dataS(names, n, unit)
Load and Select n data from the datasets listed in `Names`
"""
function load_dataS(names,n,unit="mm")
	T_ABB, T_JAKA, Qs_ABB,Qs_JAKA, As,Bs,Cs = load_data(names,unit);
	# @printf "%d data points were found, %d were requested" size(T_ABB,1) n
	if n > size(T_ABB,1)
		println("n Given is larger than the size of dataset...")
		return T_ABB, T_JAKA, Qs_ABB,Qs_JAKA, As,Bs,Cs
	else
		s = sample(1:size(T_ABB,1),n,replace=false)
		Asn, Bsn, Csn = unselect_dataset(As,Bs,Cs,s)
		return T_ABB[s], T_JAKA[s], Qs_ABB[s], Qs_JAKA[s], As[s], Bs[s], Cs[s], Asn, Bsn, Csn
	end
end

"""
	unselect_dataset(As,Bs,Cs,s)
Collect data points that is not selected by `s`
"""
function unselect_dataset(As,Bs,Cs,s)
	num = size(As,1);
	selector = 1:num;
	in_s = [~(selector[i] in s) for i = 1:num]
	if sum(in_s) == num
		selector = [1]
	else
		selector = selector[in_s]
	end
	n = minimum([size(selector,1),50])
	k = sample(1:size(selector,1),n,replace=false)
	return As[selector[k]],Bs[selector[k]],Cs[selector[k]]
end

"""
	load_data(names, unit)
Load data from the datasets listed in `Names`
"""
function load_data(names,unit="mm",prefix = ".")
	T_ABB = Vector{SMatrix{4,4,Float64,16}}(undef,0);
	T_JAKA =  Vector{SMatrix{4,4,Float64,16}}(undef,0);
	Qs_ABB = Vector{SVector{6,Float64}}(undef,0);
	Qs_JAKA = Vector{SVector{6,Float64}}(undef,0); 
	As = Vector{SMatrix{4,4,Float64,16}}(undef,0);
	Bs = Vector{SMatrix{4,4,Float64,16}}(undef,0);
	Cs = Vector{SMatrix{4,4,Float64,16}}(undef,0);
	for name in names
		data_dir = prefix * "/traj/"*name; Bs_dir = prefix * "/result/"*name;
		T_ABB₁,T_JAKA₁,Qs_ABB₁,Qs_JAKA₁ = read_Data(data_dir);
		# T_ABB₁ = abb.fkine.(Qs_ABB₁); T_JAKA₁ = jaka.fkine.(Qs_JAKA₁)
		As₁ = T_JAKA₁; Cs₁ = T_ABB₁;
		file = matopen(Bs_dir*"/Ts.mat");
		Bs₁ = read(file, "Ts");
		Bs₁ = cvtSMatrix(Bs₁)
		close(file);
		good = clean_up(T_JAKA₁,Bs₁,T_ABB₁,prefix);
		append!(T_ABB,T_ABB₁[good]); append!(T_JAKA,T_JAKA₁[good]);
		append!(Qs_ABB,Qs_ABB₁[good]); append!(Qs_JAKA,Qs_JAKA₁[good]);
		append!(As,As₁[good]); append!(Cs,Cs₁[good]); append!(Bs,Bs₁[good]);
	end
	if unit == "m"
		T_ABB = tometer.(T_ABB); T_JAKA = tometer.(T_JAKA);
		As = tometer.(As); Bs = tometer.(Bs); Cs = tometer.(Cs)
	end
	return T_ABB, T_JAKA, Qs_ABB,Qs_JAKA, As, Bs, Cs
end