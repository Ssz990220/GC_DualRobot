"""
	scale_motorS(M::Union{SMatrix{1,8,Float64,8},SVector{8,Float64}}, scale)
Scale the translational part of a motor by `scale`. The function can be used to convert unit of length for a motor.
"""
@inline function scale_motorS(M::Union{SMatrix{1,8,Float64,8},SVector{8,Float64}}, scale)
	return @SVector [M[1],M[2],M[3],M[4]*scale,M[5],M[6]*scale,M[7]*scale,M[8]*scale];
end

"""
	tfm2motorS(T)
Convert a transformation matrix to a motor. 

...
# Argument
- `T::SMatrix{4,4,Float64,16}`: Transformation matrix. The type of the matrix must be *SMatrix{4,4,Float64,16}*.

# return
- `M::SVector{8,Float64}`: Motor. The type of the motor is *SVector{8,Float64}*.

See also [`motorS2tfm`](@ref).
"""
function tfm2motorS(T)
	return ecmul(rm2rotorS(T[SOneTo(3),SOneTo(3)]),translatorS(T[SOneTo(3),4]))
end
"""
	motorS2tfm(M::SVector)
Convert a motor to a transformation matrix.

...
# Argument
- `M::SVector{8,Float64}`: Motor. The type of the motor is *SVector{8,Float64}*.

# return
- `T::SMatrix{4,4,Float64,16}`: Transformation matrix. The type of the matrix is *SMatrix{4,4,Float64,16}*.

See also [`tfm2motorS`](@ref).
"""
function motorS2tfm(M::SVector)
	tₛ = motorS2tₛ(M);
	rₛ = rotorS2rm(M);
	return vcat(hcat(rₛ,(@SMatrix [tₛ[1];tₛ[2];tₛ[3]])),@SMatrix [0.0 0.0 0.0 1.0]);
end

"""
	rotorS2rm(r::SVector)
Convert a rotor (in PGA) to a rotation matrix.

...
# Argument
- `r::SVector{8,Float64}`: Rotor in PGA. e.g. `@SVector [1,0,0,0,0,0,0,0]`

# return
- `R::SMatrix{3,3,Float64,9}`: Rotation matrix.

See also [`rotorG32rm`](@ref), [`rm2rotorS`](@ref), [`rm2rotorG3S`](@ref).
"""
function rotorS2rm(r::SVector)
	Quat2rm(rotorS2quat(r))
end

"""
	rotorG32rm(r::SVector)

Convert a rotor (in G3) to a rotation matrix.

...
# Argument
- `r::SVector{4,Float64}`: Rotor in G3. e.g. `@SVector [1,0,0,0]`

# return
- `R::SMatrix{3,3,Float64,9}`: Rotation matrix.

See also [`rotorS2rm`](@ref), [`rm2rotorS`](@ref), [`rm2rotorG3S`](@ref).
"""
function rotorG32rm(r::SVector)
	Quat2rm(rotorG32quat(r))
end

"""
	rm2rotorS(R::Union{Matrix,Vector,SVector,SMatrix,Rotation})

Convert a rotation matrix to a rotor (in PGA).

...
# Argument
- `R::Union{Matrix,Vector,SVector,SMatrix,Rotation}`: Rotation matrix.

# return
- `r::SVector{8,Float64}`: Rotor in PGA.

See also [`rm2rotorG3S`](@ref), [`rotorS2rm`](@ref), [`rotorG32rm`](@ref).
"""
function rm2rotorS(R::Union{Matrix,Vector,SVector,SMatrix,Rotation})
	R_ = R[SOneTo(3),SOneTo(3)];
	quat = rm2Quat(R_);
	return @SVector [quat[1],quat[4],-quat[3],0.0,quat[2],0.0,0.0,0.0] 
end

"""
	rm2rotorG3S(R::Union{Matrix,Vector,SVector,SMatrix,Rotation})

Convert a rotation matrix to a rotor (in G3).

...
# Argument
- `R::Union{Matrix,Vector,SVector,SMatrix,Rotation}`: Rotation matrix.

# return
- `r::SVector{4,Float64}`: Rotor in G3.

See also [`rm2rotorS`](@ref), [`rotorS2rm`](@ref), [`rotorG32rm`](@ref).
"""
function rm2rotorG3S(R::Union{Matrix,Vector,SVector,SMatrix,Rotation})
	R_ = R[SOneTo(3),SOneTo(3)];
	quat = rm2Quat(R_);
	return @SVector [quat[1],quat[4],-quat[3],quat[2]] 
end

function MᵢⱼS(θ,i,DH)	
	c₁ = cos(1/2*θ - 1/2*DH[i,4]); s₁ = sin(1/2*θ - 1/2*DH[i,4]);
	c₂ = cos(1/2*θ + 1/2*DH[i,4]); s₂ = sin(1/2*θ + 1/2*DH[i,4]);
	return @SVector [0.5*c₁ + 0.5*c₂,
	0.5*s₂ + 0.5*s₁,
	-0.5*c₁ + 0.5*c₂,
	-0.25*DH[i,2]*c₁ - 0.25*DH[i,2]*c₂ + 0.25*DH[i,3]*c₁ - 0.25*DH[i,3]*c₂,
	(0.5*s₂ - 0.5*s₁),
	-0.25*DH[i,2]*s₂ - 0.25*DH[i,2]*s₁ - 0.25*DH[i,3]*s₂ + 0.25*DH[i,3]*s₁,
	0.25*DH[i,2]*c₁ - 0.25*DH[i,2]*c₂ - 0.25*DH[i,3]*c₁ - 0.25*DH[i,3]*c₂,
	-0.25*DH[i,2]*s₂ + 0.25*DH[i,2]*s₁ - 0.25*DH[i,3]*s₂ - 0.25*DH[i,3]*s₁]
end
"""
	FkineG3S(θₛ,DH)
Forward kinematics of robot with PGA algebra.
Result is represented in *Motor* form.
"""
function FkineGAS(θₛ,DH)
	M = @SVector [1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
	@inbounds for i = 1:size(DH,1)
		M = ecmul(MᵢⱼS(θₛ[i]+DH[i,1],i,DH),M)
	end
	return M
end

@inline rotorG3X(θ) = rotorG3S(θ,@SVector [1.0,0.0,0.0])
@inline rotorG3Z(θ) = rotorG3S(θ,@SVector [0.0,0.0,1.0])
"""
	FkineG3S(θₛ,DH)
Forward kinematics of robot with G3 algebra.
Rotation part and translation part are solved seperately.
Rotation is represented in *rotor* and trasnlation is represented in *vector*.
"""
function FkineG3S(θₛ,DH)
	local t = @SVector [0.0,0.0,0.0]		# Init Translation
	local R = @SVector [1.0,0.0,0.0,0.0]	# Init Rotation
	local x = @SVector [1.0,0.0,0.0]
	local z = @SVector [0.0,0.0,1.0]
	@inbounds for n = 1:size(DH,1)
		i = size(DH,1) - n + 1
		RZ = rotorG3Z(θₛ[i]+DH[i,1]); RX = rotorG3X(DH[i,4]);
		R = ecmulG3(R,RX,RZ)
		t = rotPoint(RZ,rotPoint(RX,t + x*DH[i,2]))+DH[i,3]*z
	end
	return R,t
end

"""
	rotorG3Norm(X)

Normalize a rotor `X` in G3.
"""
@inline function rotorG3Norm(X::SVector{4,Float64})
	return X./norm(X)
end

"""
	motorNormS(X)

Normalize a motor `X` in PGA.
"""
function motorNormS(X::SVector{8,Float64})
	core = ecmul(X,ereversion(X))
	return ecmul(X/sqrt(core[1]),(@SVector [1.0,0.0,0.0,0.0,0.0,0.0,0.0,-core[8]/(2*core[1])]));
end


"""
	rotPoint(r::SVector{4,Float64},p::SVector{3,Float64})
Rotate a point `p=xv₁+yv₂+zv₃` by a rotor `r=r₀+r₁v₁₂+r₂v₁₃+r₃v₂₃` with sandwich operation:
	`p' = r̃pr`

...
# Argument
- `r::SVector{4,Float64}` : rotor
- `p::SVector{3,Float64}` : point

Formula is simplified by symbolic algebra.
"""
@inline function rotPoint(r::SVector{4,Float64},p::SVector{3,Float64})
	@SVector[p[1]*r[1]^2 - p[1]*r[2]^2 - p[1]*r[3]^2 + p[1]*r[4]^2 - 2*p[2]*r[1]*r[2] - 2*p[2]*r[3]*r[4] - 2*p[3]*r[1]*r[3] + 2*p[3]*r[2]*r[4],
	2*p[1]*r[1]*r[2] - 2*p[1]*r[3]*r[4] + p[2]*r[1]^2 - p[2]*r[2]^2 + p[2]*r[3]^2 - p[2]*r[4]^2 - 2*p[3]*r[1]*r[4] - 2*p[3]*r[2]*r[3],
	2*p[1]*r[1]*r[3] + 2*p[1]*r[2]*r[4] + 2*p[2]*r[1]*r[4] - 2*p[2]*r[2]*r[3] + p[3]*r[1]^2 + p[3]*r[2]^2 - p[3]*r[3]^2 - p[3]*r[4]^2]
end