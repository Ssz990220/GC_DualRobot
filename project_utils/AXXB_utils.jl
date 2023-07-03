# This function implements the AX=XB algorithm
# based on the paper "Hand-Eye Calibration Using Dual Quaternions" by Daniilidis,
# and the paper "Motor algebra for 3D kinematics: the case of the hand-eye calibration" by Bayro-Corrochano.
# The algorithm uses dual quaternion/motor to solve the hand-eye calibration in closed form.

function AXXBs(As::Vector{SVector{8,Float64}},Bs::Vector{SVector{8,Float64}})
	n = size(As,1);
	a = Vector{SVector{3,Float64}}(undef,n-1);
	a_ = Vector{SVector{3,Float64}}(undef,n-1);
	b = Vector{SVector{3,Float64}}(undef,n-1);
	b_ = Vector{SVector{3,Float64}}(undef,n-1);
	@inbounds for i = 1:n-1
		A = ecmul(As[i],ereversion(As[i+1]));
		a[i] = elineS(A); a_[i] = emomentS(A);
		B = ecmul(ereversion(Bs[i]),Bs[i+1])
		b[i] = elineS(B); b_[i] = emomentS(B);
	end
	
	D = zeros(6*(n-1),8)
	@inbounds for i in 1:(n-1)
		D[6*(i-1)+1:6*(i-1)+6,:] = [(a[i] - b[i]) mcross(a[i] + b[i]) zeros(3,4);
									(a_[i] - b_[i]) mcross(a_[i] + b_[i]) (a[i] - b[i]) mcross(a[i] + b[i])]
	end
	U, Σ, V = svd(D);
	V = SMatrix{8,2,Float64,16}(V[:,7:8]);
	α,β = solve_αβ(V[SOneTo(8),1],V[SOneTo(8),2])
	return postproc_axxbS(α,β,V[SOneTo(8),1],V[SOneTo(8),2])
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

function solv_quadS(V)
	α,β = solve_αβ(V[SOneTo(8),7],V[SOneTo(8),8])
	return postproc_axxbS(α,β,V)
end

@inline function postproc_axxbS(α,β,V₁,V₂)
	R = α.*V₁.+β.*V₂
	return @SVector [R[1],R[4],-R[3],-R[6],R[2],-R[7],-R[8],R[5]]
end


function axxb_movements_noiseₛ(n,angle_noise, dis_noise,unit="mm",fix = false)
	Aₘ,Bₘ = axxb_movementsₛ(n,Xₘ,Yₘ,Zₘ,unit);
	Aₘₙ = add_noiseSₘ.(Tₛ,angle_noise,dis_noise,fix);
	Bₘₙ = add_noiseSₘ.(Tₒ,angle_noise,dis_noise,fix);
	return Aₘₙ,Bₘₙ
end

function axxb_movementsₛ(n,Xₘ,Yₘ,Zₘ,unit="mm")
	# Hold C in still and move A randomly
	Q₁ = [SVector{6,Float64}(rand(Uniform(-π,π),1,6)) for i = 1:n]
	Qc = SVector{6,Float64}(rand(Uniform(-π,π),1,6));
	ZCYₘ = ecmul(Zₘ,UR10.fkineGAS(Qc,unit),Yₘ);
	Aₘ = UR10.fkineGAS.(Q₁,unit);
	Bₘ = ecmul.([ZCYₘ],ereversion.(ecmul.([Xₘ],Aₘ)))
	return Aₘ,Bₘ
end