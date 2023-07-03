"""
	AXBYCZ_G3(V,P,conf)
A iterative solver for AXB=YCZ problem with Geometric Algebra and Geometric Calculus.
G3 is used in this algorithm to represent rotations.

...
# Arguments
- `V::Union{Vector{Matrix},Vector{SMatrix{1,4,Float64,4}}}`: Vector of initial `X,Y,Z` 
`V[1:3]` is `X,Y,Z` repectively.
- `P::Union{Vector{Vector{Matrix}},Vector{Vector{SMatrix{1,4,Float64,4}}}`: Vector of data sets, 
`P[1]` are `X`s, `P[2]` are `Y`s, `P[3]` are set of `C`s

**Input by SMatrix is highly recommended for best Performance**

...
"""
function AXBYCZ_G3(V::Vector{SVector{4,Float64}},P,conf)
	V,err = iter_solve_Katyusha_G3(V,P,conf);
	A = motorS2tfm.(P[1]); B = motorS2tfm.(P[2]); C = motorS2tfm.(P[3]);
	n = size(A,1)
	RA = [A[i][SOneTo(3),SOneTo(3)] for i = 1:n];
	RB = [B[i][SOneTo(3),SOneTo(3)] for i = 1:n];
	RC = [C[i][SOneTo(3),SOneTo(3)] for i = 1:n];
	tA = [A[i][SOneTo(3),4] for i = 1:n];
	tB = [B[i][SOneTo(3),4] for i = 1:n];
	tC = [C[i][SOneTo(3),4] for i = 1:n];	
	RX = rotorG32rm(V[1]); RY = rotorG32rm(V[2]); RZ = rotorG32rm(V[3]);
	X,Y,Z = T_SVD(RX,RY,RZ,RA,RB,RC,tA,tB,tC);
	# X,Y,Z = T_EKF(RX,RY,RZ,RA,RB,RC,tA,tB,tC)
	return [X,Y,Z],err
end

function iter_solve_Katyusha_G3(V,P,conf)
	V = toVG3(V)
	P = [motorS2rotorG3S.(P[i]) for i∈eachindex(P)]
	if conf.err == "NAN"
		errfunc = ErrG3ₙ
	elseif conf.err == "LS"
		errfunc = ErrG3
	elseif conf.err == "Ext"
		errfunc = conf.errfunc;
	end
	V,err = katyushaX(V,P,get_i_grad_reg_AXBYCZG3, get_full_grad_regG3, errfunc,E_normG3, conf)
	V = rotorG3Norm.(toXYZ(V))
	return V,err
end

function iter_solve_GD_G3(V,P,conf)
	V = toVG3(V)
	P = [motorS2rotorG3S.(P[i]) for i∈eachindex(P)]
	if conf.err == "NAN"
		errfunc = ErrG3ₙ
	elseif conf.err == "LS"
		errfunc = ErrG3
	elseif conf.err == "Ext"
		errfunc = conf.errfunc;
	end
	V,err = GradientDescent(V,P,get_i_grad_reg_AXBYCZG3, get_full_grad_regG3, errfunc,E_normG3, conf)
	V = rotorG3Norm.(toXYZ(V))
	return V,err

end

## Get Gradient
function get_i_grad_AXBYCZG3(V,P,i,conf)
    A = P[1][i]; B = P[2][i]; C = P[3][i];
	return get_one_gradG3(A,B,C,V)
end

function get_one_gradG3(A,B,C,V)
	E = E_ABCXYZG3(A,B,C,V)
	Y = V[StaticArrays.SUnitRange(5,8)];
	Z = V[StaticArrays.SUnitRange(9,12)];
	gradx = ecmulG3(ereversionG3(B),E,ereversionG3(A));
	grady = -ecmulG3(ereversionG3(ecmulG3(Z,C)),E);
	gradz = -ecmulG3(E,ereversionG3(ecmulG3(C,Y)));
	return vcat(gradx,grady,gradz)
end

function get_one_grad_regG3(A,B,C,V,μ)
	return get_one_gradG3(A,B,C,V).+get_regulatorG3(V,μ)
end

function get_i_grad_reg_AXBYCZG3(V,P,i,conf)
	return get_i_grad_AXBYCZG3(V,P,i,conf).+get_regulatorG3(V,conf.μ)
end

function get_regulatorG3(V,μ)
	X = V[StaticArrays.SUnitRange(1,4)];
	Y = V[StaticArrays.SUnitRange(5,8)];
	Z = V[StaticArrays.SUnitRange(9,12)];
	E₁ = E_normG3(X); E₂ = E_normG3(Y); E₃ = E_normG3(Z)
	gradx = μ.*ecmulG3((ereversionG3(E₁).+E₁),X);
	grady = μ.*ecmulG3((ereversionG3(E₂).+E₂),Y);
	gradz = μ.*ecmulG3((ereversionG3(E₃).+E₃),Z);
	return vcat(gradx,grady,gradz)
end

function get_full_gradG3(V,P,conf)
	grad = sum(get_one_gradG3.(P[1],P[2],P[3],[V]))/size(P[1],1)
	return grad
end

function get_full_grad_regG3(V,P,conf)
	μ = conf.μ
	grad = sum(get_one_grad_regG3.(P[1],P[2],P[3],[V],μ))/size(P[1],1)
	return grad
end

## Variable Step size
function get_ηₜG3(Ṽₜ::SVector,Ṽₜ₁::SVector,∇ᵥgṼₜ::SVector,∇ᵥgṼₜ₁::SVector,conf)
	m = conf.m
	a1 = Ṽₜ.-Ṽₜ₁; b1 = ∇ᵥgṼₜ.-∇ᵥgṼₜ₁;
	n1 = norm(a1);
	n2 = sum(a1.*b1)
	ηₜ = 1/m*n1^2/n2;
	return ηₜ
end

## Err Functions
E_ABCXYZG3(A,B,C,V)=ecmulG3(B,V[StaticArrays.SUnitRange(1,4)],A).-ecmulG3(V[StaticArrays.SUnitRange(9,12)],C,V[StaticArrays.SUnitRange(5,8)])

function ErrG3(V,P,∇ᵥgṼₜ,ηₜ,conf)
	return norm(E_ABCXYZG3.(P[1],P[2],P[3],[V]))
end

"""
	ErrG3ₙ(V,P,∇ᵥgṼₜ,ηₜ,conf)
Null Error Function that returns 0.
"""
function ErrG3ₙ(V,P,∇ᵥgṼₜ,ηₜ,conf)
	return 0.0
end

function E_normG3(X) 
	v1 = @SVector [1.0,0.0,0.0,0.0]
	return (ecmulG3(X,ereversionG3(X)).-v1);
end

function ErrG3_Full(V,P,∇ᵥgṼₜ,ηₜ,conf)
    Result = toXYZ(V);
    Errx = get_error(Result[1],Xₜ);
    Erry = get_error(Result[2],Yₜ);
    Errz = get_error(Result[3],Zₜ);
    norme = norm(∇ᵥgṼₜ)
    err = ErrG3(V,P,∇ᵥgṼₜ,ηₜ,conf)
    return [Errx[1] 0.0 Erry[1] 0.0 Errz[1] 0.0 ηₜ norme err];
end

# Other Functions
toVG3(Init) = vcat(Init[1], Init[2], Init[3])
function toXYZ(V)
    Result = [V[StaticArrays.SUnitRange(4*(i-1)+1,4*i)] for i = 1:3]
    return Result
end