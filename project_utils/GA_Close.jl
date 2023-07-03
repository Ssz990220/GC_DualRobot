##########################################################
####                Solve Close Form                  ####
##########################################################

"""
    GAClose(A,B,C,Xᵢ,Yᵢ,Zᵢ,Q₁,Q₂)
	
Close form AXB=YCZ in GA.
With Given initial value
"""
function GAClose(A,B,C,Xᵢ,Yᵢ,Zᵢ)
	# Aₛₘ,Bₛₘ,Cₛₘ = match_signs_ABC(A,B,C,Q₁,Q₂,Xᵢ,Yᵢ,Zᵢ,unit);
	Xₑ,Yₑ,Zₑ = get_GAClose_XYZ(A,B,C)
	Xₑ = match_signS(Xᵢ,Xₑ); Yₑ = match_signS(Yᵢ,Yₑ); Zₑ = match_signS(Zᵢ,Zₑ);
	return Xₑ,Yₑ,Zₑ
end
"""
    GAClose(A,B,C,Q₁,Q₂)
	
Close form AXB=YCZ in GA.
With Unknown Initial sign for XYZ
Sign of B should be correct in this case.
"""
function GAClose(A,B,C)
	Xₑ,Yₑ,Zₑ = get_GAClose_XYZ(A,B,C)
	return Xₑ,Yₑ,Zₑ
end

"""
    get_GAClose_XYZ(A,B,C)

Find X,Y,Z in close form based on motor algebra.
A,B,C's sign should be matched.
"""
function get_GAClose_XYZ(A::T,B::T,C::T) where {T<:Vector{SVector{8,Float64}}}
	X = get_one_Motor(A,B,C);
	Y = get_one_Motor(ereversion.(A),C,B);
	Z = get_one_Motor(C,ereversion.(B),A)
	return X,Y,Z
end

function get_one_Motor(A,B,C)
	n = size(A,1)
	WAB = zeros(8*n,8); WC = zeros(8*n,64);
	for i = 1:n
		WAB[8*i-7:8*i,:] = Mₗ(B[i])*Mᵣ(A[i]);
		WC[8*i-7:8*i,:] = get_Wc(C[i]);
	end
	WABC = [WAB -WC];
	F = eigen(WABC'*WABC)
	V1 = F.vectors;
	return motorNormS(SVector{8,Float64}(V1[1:8,1]))
end

function get_Wc(c::Union{Vector,SVector{8,Float64}})
    # DO NOT TRY TO READ LINE BELOW
    # IT IS GENERATED BY MAPLE IN ADVANCE AND LONGER THAN YOU THOUGHT
    m = [c[1] -c[2] -c[3] 0.0 -c[5] 0.0 0.0 0.0 -c[2] -c[1] c[5] 0.0 -c[3] 0.0 0.0 0.0 -c[3] -c[5] -c[1] 0.0 c[2] 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -c[5] c[3] -c[2] 0.0 -c[1] 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
	c[2] c[1] -c[5] 0.0 c[3] 0.0 0.0 0.0 c[1] -c[2] -c[3] 0.0 -c[5] 0.0 0.0 0.0 c[5] -c[3] c[2] 0.0 c[1] 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -c[3] -c[5] -c[1] 0.0 c[2] 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
	c[3] c[5] c[1] 0.0 -c[2] 0.0 0.0 0.0 -c[5] c[3] -c[2] 0.0 -c[1] 0.0 0.0 0.0 c[1] -c[2] -c[3] 0.0 -c[5] 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 c[2] c[1] -c[5] 0.0 c[3] 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
	c[4] c[6] c[7] c[1] -c[8] -c[2] -c[3] -c[5] -c[6] c[4] -c[8] -c[2] -c[7] -c[1] c[5] -c[3] -c[7] c[8] c[4] -c[3] c[6] -c[5] -c[1] c[2] c[1] -c[2] -c[3] 0.0 -c[5] 0.0 0.0 0.0 -c[8] -c[7] c[6] -c[5] -c[4] c[3] -c[2] -c[1] c[2] c[1] -c[5] 0.0 c[3] 0.0 0.0 0.0 c[3] c[5] c[1] 0.0 -c[2] 0.0 0.0 0.0 -c[5] c[3] -c[2] 0.0 -c[1] 0.0 0.0 0.0;
	c[5] -c[3] c[2] 0.0 c[1] 0.0 0.0 0.0 c[3] c[5] c[1] 0.0 -c[2] 0.0 0.0 0.0 -c[2] -c[1] c[5] 0.0 -c[3] 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 c[1] -c[2] -c[3] 0.0 -c[5] 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
	c[6] -c[4] c[8] c[2] c[7] c[1] -c[5] c[3] c[4] c[6] c[7] c[1] -c[8] -c[2] -c[3] -c[5] c[8] c[7] -c[6] c[5] c[4] -c[3] c[2] c[1] -c[2] -c[1] c[5] 0.0 -c[3] 0.0 0.0 0.0 -c[7] c[8] c[4] -c[3] c[6] -c[5] -c[1] c[2] c[1] -c[2] -c[3] 0.0 -c[5] 0.0 0.0 0.0 c[5] -c[3] c[2] 0.0 c[1] 0.0 0.0 0.0 c[3] c[5] c[1] 0.0 -c[2] 0.0 0.0 0.0;
	c[7] -c[8] -c[4] c[3] -c[6] c[5] c[1] -c[2] -c[8] -c[7] c[6] -c[5] -c[4] c[3] -c[2] -c[1] c[4] c[6] c[7] c[1] -c[8] -c[2] -c[3] -c[5] -c[3] -c[5] -c[1] 0.0 c[2] 0.0 0.0 0.0 c[6] -c[4] c[8] c[2] c[7] c[1] -c[5] c[3] -c[5] c[3] -c[2] 0.0 -c[1] 0.0 0.0 0.0 c[1] -c[2] -c[3] 0.0 -c[5] 0.0 0.0 0.0 -c[2] -c[1] c[5] 0.0 -c[3] 0.0 0.0 0.0;
	c[8] c[7] -c[6] c[5] c[4] -c[3] c[2] c[1] c[7] -c[8] -c[4] c[3] -c[6] c[5] c[1] -c[2] -c[6] c[4] -c[8] -c[2] -c[7] -c[1] c[5] -c[3] c[5] -c[3] c[2] 0.0 c[1] 0.0 0.0 0.0 c[4] c[6] c[7] c[1] -c[8] -c[2] -c[3] -c[5] -c[3] -c[5] -c[1] 0.0 c[2] 0.0 0.0 0.0 c[2] c[1] -c[5] 0.0 c[3] 0.0 0.0 0.0 c[1] -c[2] -c[3] 0.0 -c[5] 0.0 0.0 0.0];
    return m
end

function toVec(y::Vector,z::Vector)
	return [y[i]*z[j] for i ∈ eachindex(y) for j ∈ eachindex(z)]
end

##########################################################
####                  Sign Fixing                     ####
##########################################################

function match_sign_norm(A,B)
	Aₘₛ = norm(A-B) < norm(A+B) ? A : -A
	return Aₘₛ
end

