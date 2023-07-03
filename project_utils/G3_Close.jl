##
function match_sign_ABC(A,B,C,Q₁,Q₂,robot1,robot2,unit="mm")
	Aₛ = match_signS.(robot1.fkineGAS.(Q₁,unit),tfm2motorS.(A));
	Cₛ = match_signS.(robot2.fkineGAS.(Q₂,unit),tfm2motorS.(C));
	Bₛ = tfm2motorS.(B);
	N = size(Aₛ,1);
	## Find a rough guess on X Y Z
	n = min(10,N);
	local X,Y,Z
	while true
		samples = sample(1:N,n,replace=false)
		Aₛₛ = Aₛ[samples]; Bₛₛ = Bₛ[samples]; Cₛₛ = Cₛ[samples];
		idx = find_best_5_Bcs(Aₛₛ,Bₛₛ,Cₛₛ)
		err = zeros(5)
		for i = 1:5
			Bₖ = correct_B_Sign(Bₛₛ,idx[i])
			X,Y,Z = get_XYZ_rotor_Close(Aₛₛ,Bₖ,Cₛₛ)
			err[i] = get_error_XYZ_Close(Aₛₛ,Bₖ,Cₛₛ,X,Y,Z)
		end
		if minimum(err) < 1
			# println(err)
			Bₖ = correct_B_Sign(Bₛₛ,idx[argmin(err)])
			X,Y,Z = get_XYZ_rotor_Close(Aₛₛ,Bₖ,Cₛₛ)
			break
		end
	end
	## Use X Y Z to match B's Sign
	Bₛ = match_B_signS.(Aₛ,Bₛ,Cₛ,[X],[Y],[Z]);
	return Aₛ,Bₛ,Cₛ,X,Y,Z
end

function get_error_XYZ_Close(A,B,C,X,Y,Z)
	A = motorS2rotorG3S.(A); B = motorS2rotorG3S.(B); C = motorS2rotorG3S.(C);
	err = sum(norm.(ecmulG3.(B,[X],A).-ecmulG3.([Z],C,[Y])));
	return err;
end

function get_one_rotor(P)
	a = motorS2rotorG3S.(P[1]);
	b = motorS2rotorG3S.(P[2]);
	c = motorS2rotorG3S.(P[3]);
	WAB,WC = get_WAB_WC(a,b,c);
	WABₘ = cat_AB(WAB,0);
	WCₛ = cat_C(WC);
	WABCₘ = hcat(WABₘ,WCₛ);
	F = eigen(WABCₘ'*WABCₘ)
	V1 = F.vectors;
	# return motorNorm(vector2rotor(V1[1:4,1]))
	return SVector{4,Float64}(V1[1:4,1]/norm(V1[1:4,1]))
end

function find_best_5_Bcs(Aₛ,Bₛ,Cₛ)
	## Aₛ,Cₛ's sign should be calibrated.
	## This function returns 5 best B's Signs' combination coded in bit form
	a = motorS2rotorG3S.(Aₛ);
	b = motorS2rotorG3S.(Bₛ);
	c = motorS2rotorG3S.(Cₛ);
	n = size(a,1)
	WAB,WC = get_WAB_WC(a,b,c);
	WCₛ = cat_C(WC);
	min_Singular = zeros(2^n)
	@inbounds @simd for j = 0:(2^n-1)
		WABₛ = cat_AB(WAB,j);
		WABC = hcat(WABₛ,-WCₛ);
		u,s,k = svd(WABC)
		min_Singular[j+1] = s[end];
	end
	idx = sortperm(min_Singular)
	return idx[1:5].-1      # bit index starts from 0, so -1
end

function correct_B_Sign(Bₛ,bitidx)
	n = size(Bₛ,1);
	return [Bₛ[i]*(2*bitget(bitidx,i)-1) for i = 1:n]
end

function get_Bₛ(Aₛ,Bₛ,Cₛ)
	idx = find_best_5_Bcs(Aₛ,Bₛ,Cₛ)
	Bₛ = correct_B_Sign(Bₛ,idx[1])
	return Bₛ
end

function get_XYZ_rotor_Close(A::T,B::T,C::T) where {T<:Vector{SVector{8,Float64}}}
	X = get_one_rotor([A,B,C])
	Y = get_one_rotor([ereversion.(A),C,B]);
	Z = get_one_rotor([C,ereversion.(B),A]);
	return X,Y,Z
end

## Method Utils
function bitget(n,i)
	D = digits(Int(n), base = 2)
	if i > length(D)
		return 0
	else
		return D[i]
	end
end

function cat_AB(WAB,k)
	n = size(WAB,1);
	WABₛ = zeros(n*4,4);
	for i = 1:n
		if bitget(k,i)==0
			WABₛ[4*(i-1)+1:4*i,:] = WAB[i]
		else
			WABₛ[4*(i-1)+1:4*i,:] = -WAB[i]
		end
	end
	return WABₛ
end

function cat_C(WC)
	n = size(WC,1);
	WCₛ = zeros(n*4,16);
	for i = 1:n
		WCₛ[4*(i-1)+1:4*i,:] = WC[i]
	end
	return WCₛ
end

function get_WAB_WC(A,B,C)
	n = size(A,1);
	WAB = Rₗ.(B).*Rᵣ.(A);
	WC = get_Wc.(C)
	return WAB, WC
end

function get_Wc(c)
	return [c[1] -c[2] -c[3] -c[4] -c[2] -c[1] c[4] -c[3] -c[3] -c[4] -c[1] c[2] -c[4] c[3] -c[2] -c[1];
c[2] c[1] -c[4] c[3] c[1] -c[2] -c[3] -c[4] c[4] -c[3] c[2] c[1] -c[3] -c[4] -c[1] c[2];
c[3] c[4] c[1] -c[2] -c[4] c[3] -c[2] -c[1] c[1] -c[2] -c[3] -c[4] c[2] c[1] -c[4] c[3];
c[4] -c[3] c[2] c[1] c[3] c[4] c[1] -c[2] -c[2] -c[1] c[4] -c[3] c[1] -c[2] -c[3] -c[4]];
end

"""
    match_B_signS(A::T,B::T,C::T,X::K,Y::K,Z::K) where {T<:SVector{8,Float64},K<:SVector{4,Float64}}

A,C's signs are correct and X,Y,Z's signs are given and fixed.
B's sign will be found by Comparing ||BXA-ZCY|| and ||BXA+ZCY||
"""
function match_B_signS(A::T,B::T,C::T,X::K,Y::K,Z::K) where {T<:SVector{8,Float64},K<:SVector{4,Float64}}
	Ar = motorS2rotorG3S(A);
	Br = motorS2rotorG3S(B);
	Cr = motorS2rotorG3S(C);
	Bₛ = norm(ecmulG3(Br,X,Ar) .- ecmulG3(Z,Cr,Y)) < norm(ecmulG3(Br,X,Ar) .+ ecmulG3(Z,Cr,Y)) ? B : -B
	return Bₛ
end