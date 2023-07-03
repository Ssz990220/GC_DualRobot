# Kinematics for Robotics

function FkineDH(θ,DH)
	T = eye(4);
	n = size(DH,1)
	@inbounds for i in 1:n
		T = T*Tᵢⱼ(θ[i],i,DH)
	end
	return T
end

function Tᵢⱼ(θᵢ,i,DH)
	θ = θᵢ + DH[i,1];
	Tᵢⱼ₁=@SMatrix [cos(θ) 	-sin(θ) 	0 	0;
		 sin(θ) 	cos(θ) 	0 	0;
		 0 			0 			1 	DH[i,3];
		 0 			0 			0 	1];
	Tᵢⱼ₂=@SMatrix [1 0 				0 				DH[i,2];
		  0 cos(DH[i,4]) 	-sin(DH[i,4])  	0;
		  0 sin(DH[i,4]) 	cos(DH[i,4]) 	0;
		  0 0 				0 				1];
	return Tᵢⱼ₁*Tᵢⱼ₂
end

function ikine_iter(Te,DH,q0=rand(1,6),TOL::Float64=1e-9;)
	i=0;
	jacobe(q) = Jacobe(q,DH)
	e = zeros(6,1);
	T = FkineDH(q0,DH)
	q = q0;
	e = Inf*ones(6,1);
	while maximum(abs.(e))>TOL
		err = MatrixLog6(T\Te);
		e[1:3,1]=err[1:3,4];
		e[4:6,1] = dmcross(err[1:3,1:3])
		JB=jacobe(q);
		dq = JB\e;
		q = q + dq';
		i=i+1;
		T = FkineDH(q,DH)
		if i > 200
			println(maximum(abs.(e)))
			error("Max Iteration Reached")
		end
	end
	return wraptoπ.(q)
end

function Jacobe(qs,DH)
	J = zeros(6,6)
	Tb = FkineDH(qs,DH)
	pe = Tb[1:3,4]
	Ts = zeros(4,4,6)
	for i = 1:6
		if i == 1	
			Ts[:,:,1] = eye(4)
		else
			Ts[:,:,i] = Ts[:,:,i-1]*Tᵢⱼ(qs[i-1],i-1,DH)
		end
		z = Ts[1:3,3,i]
		p = Ts[1:3,4,i]
		Jp = mcross(z)*(pe-p)
		Jo=z;
		J[:,i] = [Jp;Jo];
	end
	Rb = Tb[1:3,1:3];
	JB = [Rb' zeros(3,3);zeros(3,3) Rb']*J;
	return JB
end

@inline tfrotx(θ) = @SMatrix [1 0 0 0;0 cos(θ) -sin(θ) 0;0 sin(θ) cos(θ) 0;0 0 0 1]
@inline tfroty(θ) = @SMatrix [cos(θ) 0 sin(θ) 0;0 1 0 0;-sin(θ) 0 cos(θ) 0;0 0 0 1]
@inline tfrotz(θ) = @SMatrix [cos(θ) -sin(θ) 0 0;sin(θ) cos(θ) 0 0;0 0 1 0;0 0 0 1]
@inline rotx(θ) = @SMatrix [1 0 0;0 cos(θ) -sin(θ);0 sin(θ) cos(θ)];
@inline roty(θ) = @SMatrix [cos(θ) 0 sin(θ);0 1 0;-sin(θ) 0 cos(θ)];
@inline rotz(θ) = @SMatrix [cos(θ) -sin(θ) 0;sin(θ) cos(θ) 0;0 0 1];
@inline tfxyz((x,y,z)) = @SMatrix [1.0 0 0 x;0 1 0 y;0 0 1 z;0 0 0 1]
@inline mcross(x) = @SMatrix [0 -x[3] x[2];x[3] 0 -x[1];-x[2] x[1] 0]
@inline dmcross(ω̂) = @SMatrix [ω̂[3,2] ω̂[1,3] ω̂[2,1]];
@inline invT(T) = vcat(hcat((T[SOneTo(3),SOneTo(3)])',-(T[SOneTo(3),SOneTo(3)])'*T[SOneTo(3),4]),@SMatrix [0.0 0.0 0.0 1.0]);
@inline tfm2R(T) = T[SOneTo(3),SOneTo(3)];
@inline function SE3Adjoint(T)
	R,P = T2rt(T);
	Ad = vcat(hcat(R,(@SMatrix zeros(3,3))),hcat(mcross(P)*R,R))
	return Ad
end
@inline function SE3AdjointInv(T)
	R,P = T2rt(T);
	AdInv = vcat(hcat(R',(@SMatrix zeros(3,3))),hcat(-mcross(R'*P)*R',R'))
end


"""
	ProjectToSO3(R)
Takes mat: A matrix `R` near SO(3) to project to SO(3).
Returns R representing the closest rotation matrix that is in SO(3).
This function uses singular-value decomposition (see [this link](http://hades.mech.northwestern.edu/index.php/Modern_Robotics_Linear_Algebra_Review))
and is only appropriate for matrices close to SO(3).
"""
function ProjectToSO3(R::T) where {T}
	if !(T<: SMatrix)
		@warn "R is not an SMatrix, use SMatrix for best performance"
	end
	U,Σ,V = svd(R);
	R2 = U*V';
	if det(R2) < 0
		R2 = hcat(R2[:,SOneTo(2)],R2[:,3]);
	end
	return R2
end


function ProjectToSE3(T::SMatrix{4,4,Float64,16})
	R,t = T2rt(T);
	R_ = ProjectToSO3(R);
	return rt2T(R_,t)
end


function wraptoπ(a)
	b = a
	if b>0
		n = floor(a/2π);
		b = a-n*2π
		if b > π
			b = b - 2π
		end
	else
		n = floor(a/(-2π))
		b = a-n*(-2π)
		if b < -π
			b = b+2π
		end
	end
	return b
end

function wrapto1(a)
	if a > 1
		return 1
	elseif a < -1
		return -1
	else
		return a
	end
end
"""
	MatrixExp3(so3)
Takes a 3×3 so(3) representation of exponential coordinates.
Returns `R` ∈ SO(3) that is achieved by rotating about omghat by theta 
from an initial orientation `R = I`.

*From Modern Robotics*
"""
function MatrixExp3(so3::T) where {T}
	if !(T<: SMatrix)
		# @warn "so3 is not an SMatrix, use SMatrix for best performance"
	end
	ωθ = dmcross(so3);
	if NearZero(norm(ωθ))
		R = eye(3);
	else
		θ = norm(ωθ); ω = ωθ / θ;
		matω = so3 / θ;
		R = eye(3) .+ sin(θ) .* matω .+ (1 - cos(θ)) .* (matω * matω);
	end
	return R
end
"""
	MatrixExp6(se3)
Takes a se(3) representation of exponential coordinates.
Returns a `T` ∈ SE(3) that is achieved by traveling along/about the 
screw axis S for a distance theta from an initial configuration `T = I`.

*From Modern Robotics*
"""
function MatrixExp6(se3::T₁) where {T₁}
	if !(T₁<: SMatrix)
		# @warn "se3 is not an SMatrix, use SMatrix for best performance"
	end
	ωθ = dmcross(se3[SOneTo(3), SOneTo(3)]);
	if NearZero(norm(ωθ))
		R = eye(3);
		t = se3[SOneTo(3),4];
		T = rt2T(R,t);
	else
		θ = norm(ωθ); ω = ωθ / θ;
		matω = se3[SOneTo(3), SOneTo(3)] ./ θ;
		R = eye(3) .+ sin(θ) .* matω .+ (1 - cos(θ)) .* (matω * matω);
		t = (eye(3) .* θ + (1 - cos(θ)) .* matω + (θ - sin(θ)) .* (matω * matω)) * se3[SOneTo(3),4] ./ θ;
		T = rt2T(R,t);
	end
	return T
end

"""
	MatrixLog6(T)
Takes a transformation matrix `T` ∈ SE(3).
Returns the corresponding se(3) representation of exponential 
coordinates.

*From Modern Robotics*
"""
function MatrixLog6(T::T₁) where {T₁}
	if !(T₁<: SMatrix)
		# @warn "T is not an SMatrix, use SMatrix for best performance"
	end
	R, p = T2rt(T);
	omgmat = MatrixLog3(R);
	if isequal(omgmat, zeros(3,3))
		expmat = vcat(hcat((@SMatrix zeros(3,3)), T[SOneTo(3), 4]), @SMatrix zeros(1,4));
	else
		theta = acos(wrapto1((trace(R) - 1) / 2));
		expmat = vcat(hcat(omgmat, (eye(3) .- omgmat ./ 2 .+ (1 / theta - cot(theta / 2) / 2) .* (omgmat * omgmat) ./ theta) * p), @SMatrix zeros(1,4))
	end
	return expmat
end


"""
	T2rt(T)
Extract Rotation matrix `R` and Translation Vector `t` from Homogenous Transformation Matrix `T`
"""
function T2rt(T::T₁) where {T₁}
	if !(T₁<: SMatrix)
		# @warn "T is not an SMatrix, use SMatrix for best performance"
	end
	R = T[SOneTo(3),SOneTo(3)];
	P = T[SOneTo(3),4];
	return R,P
end

"""
	MatrixLog3(R)
Takes `R`` (rotation matrix).
Returns the corresponding so(3) representation of exponential 
coordinates.

*From Modern Robotics*
"""
function MatrixLog3(R::T) where {T}
	if !(T<: SMatrix)
		# @warn "R is not an SMatrix, use SMatrix for best performance"
	end
	acosinput = (trace(R) - 1) / 2;
	if acosinput >= 1
		so3mat = zeros(3,3);
	elseif acosinput <= -1
		if ~NearZero(1 + R[3, 3])
			omg = (1 / sqrt(2 * (1 + R[3, 3])))* [R[1,3]; R[2,3]; 1 + R[3,3]];
		elseif ~NearZero(1 + R[2, 2])
			omg = (1 / sqrt(2 * (1 + R[2, 2])))* [R[1, 2]; 1 + R[2, 2]; R[3, 2]];
		else
			omg = (1 / sqrt(2 * (1 + R[1, 1])))* [1 + R[1, 1]; R[2, 1]; R[3, 1]];
		end
		so3mat = mcross(pi * omg);
	else
		theta = acos(acosinput);
		so3mat = theta * (1 / (2 * sin(theta))) * (R .- R');
	end
	return so3mat
end

function rm2Quat(R::T) where {T}
	if !(T<: SMatrix)
		# @warn "R is not an SMatrix, use SMatrix for best performance"
	end
	R_ = t2r(R)
	Quat = QuatRotation(R_);
	return @SVector [Quat.q.s,Quat.q.v1,Quat.q.v2,Quat.q.v3]
end

"""
	vex(T)
Takes se3mat a 3×3 skew matrix or 4×4 se(3) matrix
Returns the corresponding 6-vector (representing spatial velocity).
"""
function vex(T::T₁) where{T₁}
	if !(T₁<: SMatrix)
		# @warn "T is not an SMatrix, use SMatrix for best performance"
	end
	if size(T,1) == 3		## so3 to vec(so3)
		return dmcross(T)
	else
		return hcat(dmcross(T),T[SOneTo(3),4]')
	end
end

function vee(twist::SMatrix{6,1,Float64,6})
	 return vcat(hcat(mcross(twist[SOneTo(3),1]),twist[StaticArrays.SUnitRange(4,6),1]),@SMatrix zeros(1,4));
end


Quat2rm(q) = QuatRotation(q[1],q[2],q[3],q[4])
t2r(T) = T[SOneTo(3),SOneTo(3)]
t2T(t) = vcat(hcat(eye(3),(@SMatrix [t[1];t[2];t[3]])), @SMatrix [0.0 0.0 0.0 1.0])
r2t(R) = vcat(hcat(R,(@SMatrix [0;0;0])),@SMatrix [0.0 0.0 0.0 1.0])
rt2T(R,t) = vcat(hcat(R,@SMatrix [t[1];t[2];t[3]]),@SMatrix [0.0 0.0 0.0 1.0])
angvec2r(θ,r) = exp(mcross(r)*θ)