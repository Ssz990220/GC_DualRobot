# GA Operations Definitions
# Most of the functions are defined with StaticArrays (SVectors, SMatrices) for performance reasons.

@inline rotorG3S(θ,n::SVector) = @SVector [cos(θ/2),sin(θ/2)*n[3],-sin(θ/2)*n[2],sin(θ/2)*n[1]]
@inline rotorS(θ,n::SVector) = @SVector [cos(θ/2),sin(θ/2)*n[3],-sin(θ/2)*n[2],0.0,sin(θ/2)*n[1],0.0,0.0,0.0]
@inline translatorS(tₛ::Union{Matrix,Vector,SVector,SMatrix}) = @SVector [1.0,0.0,0.0,-tₛ[1]/2,0.0,-tₛ[2]/2,-tₛ[3]/2,0.0]
@inline elineS(L::SVector) = @SVector [L[5],-L[3],L[2]];
@inline emomentS(L::SVector) = @SVector [-L[4],-L[6],-L[7]];

"""
	motorS2tₛ(M::SVector)

Extract the translation vector from a motor.
"""
function motorS2tₛ(M::SVector)
	return rotorS2vec(2*ecmul(ereversion(M),motordS2rotorS(M)))[StaticArrays.SUnitRange(2,4)]
end

"""
	Qᵣ(p)
Matrix form of the right quaternion multiplication.

q∘p = Qᵣ(p)vec(q)
"""
Qᵣ(p)=
@SMatrix [p[1] -p[2] -p[3]  -p[4];
 p[2] p[1] -p[4] p[3];
p[3] p[4]  p[1] -p[2];
 p[4] -p[3]  p[2]  p[1]]

"""
	Qₗ(p)
Matrix form of the left quaternion multiplication.

p∘q = Qᵣ(p)vec(q)
"""
 Qₗ(p)=
@SMatrix [p[1] -p[2] -p[3]  -p[4];
  p[2] p[1] p[4] -p[3];
 p[3] -p[4]  p[1] p[2];
  p[4] p[3]  -p[2]  p[1]]

#########################################################################
####                      FUNCTIONS FOR AXB=YCZ                      ####
#########################################################################

function approx(x::Float64,tol=1e-10)
	a = abs(x)>tol ? x : 0;
end

function approx(x::Union{Matrix,Vector,SVector{8,Float64},SVector{4,Float64}},tol=1e-10)
	a = approx.(x)
end
# Need Fixing Here: Dual Quaternion Form
# @inline quat2rotor(q::Union{SVector,Vector,Matrix}) = q[1]*v+q[2]*v₂₃-q[3]*v₁₃+q[4]*v₁₂
@inline rotorS2quat(x::SVector) = @SVector [x[1],x[5],-x[3],x[2]];
@inline rotorG32quat(x::SVector) = @SVector [x[1],x[4],-x[3],x[2]];
@inline rotorS2vec(x::SVector) = @SVector [x[1],x[5],-x[3],x[2]];
@inline drotorS2quat(x::SVector) =@SVector [x[8],-x[4],x[6],-x[7]]
@inline motorS2rotorS(x::SVector) = @SVector [x[1],x[2],x[3],0.0,x[5],0.0,0.0,0.0];
@inline motordS2rotorS(x::SVector) = @SVector [x[8],-x[7],x[6],0.0,-x[4],0.0,0.0,0.0];

"""
	Mₗ(p)
Matrix form of the left motor multiplication.

pq = Mₗ(p)vec(q)

See also [`Mᵣ`](@ref).
"""
@inline function Mₗ(A::SVector)
	M = @SMatrix [A[1] -A[2] -A[3] 0.0 -A[5] 0.0 0.0 0.0;
	A[2] A[1] A[5] 0.0 -A[3] 0.0 0.0 0.0;
	A[3] -A[5] A[1] 0.0 A[2] 0.0 0.0 0.0;
	A[4] -A[6] -A[7] A[1] -A[8] A[2] A[3] -A[5];
	A[5] A[3] -A[2] 0.0 A[1] 0.0 0.0 0.0;
	A[6] A[4] A[8] -A[2] -A[7] A[1] A[5] A[3];
	A[7] -A[8] A[4] -A[3] A[6] -A[5] A[1] -A[2];
	A[8] A[7] -A[6] A[5] A[4] -A[3] A[2] A[1]]
	return M
end

"""
	Mᵣ(p)
Matrix form of the right motor multiplication.

qp = Mᵣ(p)vec(q)

See also [`Mₗ`](@ref).
"""
@inline function Mᵣ(B::SVector)
	M = @SMatrix [B[1] -B[2] -B[3] 0.0 -B[5] 0.0 0.0 0.0;
	B[2] B[1] -B[5] 0.0 B[3] 0.0 0.0 0.0;
	B[3] B[5] B[1] 0.0 -B[2] 0.0 0.0 0.0;
	B[4] B[6] B[7] B[1] -B[8] -B[2] -B[3] -B[5];
	B[5] -B[3] B[2] 0.0 B[1] 0.0 0.0 0.0;
	B[6] -B[4] B[8] B[2] B[7] B[1] -B[5] B[3];
	B[7] -B[8] -B[4] B[3] -B[6] B[5] B[1] -B[2];
	B[8] B[7] -B[6] B[5] B[4] -B[3] B[2] B[1]];
	return M
end

"""
	rand_rotorS(angle_noise,fix=false,dis="normal")
Uniform Random Rotor Generator. 

...
# Arguments
- `angle_noise` : Noise level
- `fix` : True for fixed length noise but in random direction.
		  False for random length noise in random direction
- `dis` : noise distribution, "normal" for normal distribution,
		  "Uniform" for uniform distribution
...
"""
function rand_rotorS(angle_noise,fix=false,dis="uniform")
	if dis=="uniform"
		if fix
			θ = angle_noise + (rand()-0.5)*0.2*angle_noise
		else
			angle_noise == 0 ? θ = 0 : θ = rand(Uniform(0,angle_noise))
		end
	else
		θ = rand(Normal(0,angle_noise/3))
	end
	vec = (@SVector rand(3)) .-0.5; vec = vec/norm(vec);
	return rotorS(θ,vec)
end

"""
	rand_translatorS(dis_noise,fix = false, dist="uniform")
Uniform Random Translator Generator. 

...
# Arguments
- `dis_noise` : Noise level
- `fix` : True for fixed length noise but in random direction.
		  False for random length noise in random direction
- `dis` : noise distribution, "normal" for normal distribution,
			"uniform" for uniform distribution
...
"""
function rand_translatorS(dis_noise, fix=false, dist="uniform")
	vec = (@SVector rand(3)) .-0.5; vec = vec/norm(vec);
	if dist == "uniform"
		if fix
			dis = dis_noise + (rand()-0.5)*0.2*dis_noise;
		else
			dis_noise == 0 ? dis = 0 : dis = rand(Uniform(0,dis_noise))
		end
	else
		dis = rand(Normal(0,dis_noise/3))
	end
	return translatorS(dis*vec)
end
#########################################################################
####                      STATIC ARRAY FUNCTIONS                     ####
#########################################################################

## PGA ##
"""
	ecmul(A,B)

Geometric Product of two Motors.

# Argument
- `A` : Motor, preferred to be coded with `SVector{8,Float64}`
- `B` : Motor, preferred to be coded with `SVector{8,Float64}`
"""
@inline function ecmul(A::SVector,B::SVector)
	X₀ = A[1]*B[1]-A[2]*B[2]-A[3]*B[3]-A[5]*B[5];
	X₁ = A[2]*B[1]+A[1]*B[2]+A[5]*B[3]-A[3]*B[5];
	X₂ = A[3]*B[1]-A[5]*B[2]+A[1]*B[3]+A[2]*B[5];
	X₃ = A[4]*B[1]-A[6]*B[2]-A[7]*B[3]+A[1]*B[4]-A[8]*B[5]+A[2]*B[6]+A[3]*B[7]-A[5]*B[8];
	X₄ = A[5]*B[1]+A[3]*B[2]-A[2]*B[3]+A[1]*B[5];
	X₅ = A[6]*B[1]+A[4]*B[2]+A[8]*B[3]-A[2]*B[4]-A[7]*B[5]+A[1]*B[6]+A[5]*B[7]+A[3]*B[8];
	X₆ = A[7]*B[1]-A[8]*B[2]+A[4]*B[3]-A[3]*B[4]+A[6]*B[5]-A[5]*B[6]+A[1]*B[7]-A[2]*B[8];
	X₇ = A[8]*B[1]+A[7]*B[2]-A[6]*B[3]+A[5]*B[4]+A[4]*B[5]-A[3]*B[6]+A[2]*B[7]+A[1]*B[8];
	return @SVector [X₀,X₁,X₂,X₃,X₄,X₅,X₆,X₇]
end

function ecmul(A::SVector,B::SVector,C::SVector)
	ecmul(ecmul(A,B),C)
end

"""
	ereversion(A)

Reversion of a Motor, ̃A.

# Argument
- `A` : Motor, preferred to be coded with `SVector{8,Float64}`
"""
@inline function ereversion(X::SVector)
    return @SVector [X[1],-X[2],-X[3],-X[4],-X[5],-X[6],-X[7], X[8]]
end

@inline function eMₗ(A)
	M = @SMatrix [A[1] -A[2] -A[3] 0.0 -A[5] 0.0 0.0 0.0;
		A[2] A[1] A[5] 0.0 -A[3] 0.0 0.0 0.0;
		A[3] -A[5] A[1] 0.0 A[2] 0.0 0.0 0.0;
		A[4] -A[6] -A[7] A[1] -A[8] A[2] A[3] -A[5];
		A[5] A[3] -A[2] 0.0 A[1] 0.0 0.0 0.0;
		A[6] A[4] A[8] -A[2] -A[7] A[1] A[5] A[3];
		A[7] -A[8] A[4] -A[3] A[6] -A[5] A[1] -A[2];
		A[8] A[7] -A[6] A[5] A[4] -A[3] A[2] A[1]];
	return M
end


## G3 ##

"""
	motorS2rotorG3S(X)

Extract the rotor part from a Motor in SVector{8,Float64}.
"""
@inline function motorS2rotorG3S(X::SVector)
	return @SVector [X[1],X[2],X[3],X[5]];
end

"""
	rotorG3S2motorS(X)

Envelope a rotor in G3 SVector{4,Float64} to a Motor in PGA SVector{8,Float64}.
"""
@inline function rotorG3S2motorS(X::SVector)
	return @SVector [X[1],X[2],X[3],0.0,X[4],0.0,0.0,0.0];
end

"""
	ecmulG3(A,B)

Geometric Product of rotors in G3.
"""
function ecmulG3(A::SVector,B::SVector)
	X₀ = A[1]*B[1]-A[2]*B[2]-A[3]*B[3]-A[4]*B[4];
	X₁ = A[2]*B[1]+A[1]*B[2]+A[4]*B[3]-A[3]*B[4];
	X₂ = A[3]*B[1]-A[4]*B[2]+A[1]*B[3]+A[2]*B[4];
	X₃ = A[4]*B[1]+A[3]*B[2]-A[2]*B[3]+A[1]*B[4];
	return @SVector [X₀,X₁,X₂,X₃]
	# return Rₗ(A)*B
end

function ecmulG3(A::SVector,B::SVector,C::SVector)
	return ecmulG3(ecmulG3(A,B),C)
end

"""
	ereversionG3(X)

Reversion of a rotor in G3.
"""
@inline function ereversionG3(X)
    return @SVector [X[1],-X[2],-X[3],-X[4]]
end

"""
	Rᵣ(m)

Right multiplication matrix of a rotor in G3.
r₁r₂ = Rᵣ(r₂)*vec(r₁)
"""
@inline function Rᵣ(m)
	return @SMatrix [m[1] -m[2] -m[3] -m[4];
	m[2] m[1] -m[4] m[3];
	m[3] m[4] m[1] -m[2];
	m[4] -m[3] m[2] m[1]]
end

"""
	Rₗ(m)

Left multiplication matrix of a rotor in G3.
r₁r₂ = Rₗ(r₁)*vec(r₂)
"""
@inline function Rₗ(m)
	return @SMatrix [m[1] -m[2] -m[3] -m[4];
          m[2] m[1] m[4] -m[3];
          m[3] -m[4] m[1] m[2];
          m[4] m[3] -m[2] m[1]]
end