# This function implements the AXB=YCZ algorithm based on the paper "A Dual Quaternion-Based Approach for 
# Coordinate Calibration of Dual Robots in Collaborative Motion" by Fu et al., (2020).
# Technically speaking, this is a two step method that requires the user to first solve AX=XB.
#########################################################################
####                     FUNCTIONS FOR AXB=YCZ FU                    ####
#########################################################################
"""
	AXB_YCZ_FU(A,B,C,X)

Solve AXB=YCZ problem with Fu's [A Dual Quaternion-Based Approach for Coordinate Calibration of Dual Robots in Collaborative Motion](https://ieeexplore.ieee.org/document/9072583).

...

`A,B,C<Vector{SVector{8,Float64}}`: are motors that collected in the experiment.

`X::SVector{8,Float64}` is motor that solved by AX=XB.

...

Output: `Y,Z` are transformations that represented in motors.
"""
function AXB_YCZ_FUS(A::T,B::T,C::T,X::SVector) where {T<:Vector{SVector{8,Float64}}}
	n = size(A,1);
	Mₛ = ecmul.(ecmul.(B,[X]),A);
	D = zeros(8*n,16);
	# qₘ₁ = rotorS2quat.(Mₛ); qₘ₂ = drotorS2quat.(Mₛ);
	# qₙ₁ = rotorS2quat.(C); qₙ₂ = drotorS2quat.(C);
	@inbounds for i = 1:n
		D[8*(i-1)+1:8*i,:] = [Mₗ(Mₛ[i]) Mᵣ(C[i])]
	end
	U, Σ, V = svd(D);
	V = SMatrix{16,2}(V[:,15:16])
	V₁ = vcat(rotorS2quat(V[SOneTo(8),1]),drotorS2quat(V[SOneTo(8),1]))
	V₂ = vcat(rotorS2quat(V[SOneTo(8),2]),drotorS2quat(V[SOneTo(8),2]))
	α,β = solve_αβ(V₁,V₂)
	Y = ereversion(postprocS(α,β,V₁,V₂))
	
	V₁ = vcat(rotorS2quat(V[StaticArrays.SUnitRange(9,16),1]),drotorS2quat(V[StaticArrays.SUnitRange(9,16),1]))
	V₂ = vcat(rotorS2quat(V[StaticArrays.SUnitRange(9,16),2]),drotorS2quat(V[StaticArrays.SUnitRange(9,16),2]))
	α,β = solve_αβ(V₁,V₂)
	Z = postprocS(α,β,V₁,V₂)
	return Y,Z
end

function solve_αβ(V₁,V₂)
	ω₁₁=V₁[SOneTo(4)]; 	#U1
	ω₁₂=V₁[StaticArrays.SUnitRange(5,8)]; 	#V1
	ω₂₁=V₂[SOneTo(4)]; 	#U2
	ω₂₂=V₂[StaticArrays.SUnitRange(5,8)]; 	#V2
	den(μ) = μ^2*ω₁₁'*ω₁₁+2*μ*(ω₁₁'*ω₂₁)+ω₂₁'*ω₂₁
	m = ω₁₁'*ω₁₂  #U1'*V1
	n = ω₁₁'*ω₂₂+ω₂₁'*ω₁₂ #U1'*V2+U2'*V1
	q = ω₂₁'*ω₂₂ 	#U2'*V2
	b²4ac = n^2-4*m*q
	μ = [(-n+sqrt(b²4ac))/(2*m),(-n-sqrt(b²4ac))/(2*m)]
	idx = argmax(den.(μ))
	β = sqrt(1/den(μ[idx]))
	α = β*μ[idx]
	return α,β
end

@inline function postprocS(α,β,V₁,V₂)
	R = α.*V₁.+β.*V₂
	return @SVector [R[1],R[4],-R[3],-R[6],R[2],R[7],-R[8],R[5]]
end