# This function implements the AXB=YCZ algorithm based on the paper "Simultaneous 
# calibration of multicoordinates for a dual-robot system by solving the AXB = YCZ problem" 
# by Wang et al. (2021).

#########################################################################
####                    FUNCTIONS FOR AXB=YCZ WANG                   ####
#########################################################################
"""
	AXB_YCZ_Wang(A,B,C,conf)

Solve the AXB=YCZ problem using the algorithm proposed by Wang et al. (2021).

# Argument
- `A`: Vector of Transformation Matrices in the form of Vector{SMatirx{4,4,Float64,16}}
- `B`: Vector of Transformation Matrices in the form of Vector{SMatirx{4,4,Float64,16}}
- `C`: Vector of Transformation Matrices in the form of Vector{SMatirx{4,4,Float64,16}}
- `conf`: Configuration dictionary. The following keys are supported:
	- `:max_iter::Int`: Maximum number of iterations. Must be provided.
	- `:err::String`: Error type. Must be provided.
		- "Ext": External error. The error is calculated by the error function provided by the user.
		- "ls": Least square error. ||AXB-YCZ||²
		- "NAN": Disable error calculation.
	- `:errfunc::Function`: Error function. Should be provided if conf.err is "Ext".
	- `:svd::Bool`: Use SVD to solve the translation part instead of derived from the iteration. Default is false.
"""
function AXB_YCZ_Wang(A,B,C,conf)
	Xᵢₙᵢₜ, Yᵢₙᵢₜ, Zᵢₙᵢₜ = AXB_YCZ_Close(A,B,C,conf)
	# Xᵢₙᵢₜ[1:3,4] = zeros(1,3);Yᵢₙᵢₜ[1:3,4] = zeros(1,3);Xᵢₙᵢₜ[1:3,4] = zeros(1,3);
	V = toSM.([Xᵢₙᵢₜ, Yᵢₙᵢₜ, Zᵢₙᵢₜ]); P = [A,B,C];
	Result,err = iter_solve_Katyusha(V,P,conf)
	if haskey(conf,:svd)
		if conf.svd
			n = size(A,1)
			RA = [A[i][SOneTo(3),SOneTo(3)] for i = 1:n];
			RB = [B[i][SOneTo(3),SOneTo(3)] for i = 1:n];
			RC = [C[i][SOneTo(3),SOneTo(3)] for i = 1:n];
			tA = [A[i][SOneTo(3),4] for i = 1:n];
			tB = [B[i][SOneTo(3),4] for i = 1:n];
			tC = [C[i][SOneTo(3),4] for i = 1:n];
			Rx_sln = ProjectToSO3(Result[SOneTo(3),SOneTo(3)]);
			Ry_sln = ProjectToSO3(Result[SOneTo(3),StaticArrays.SUnitRange(4,6)]);
			Rz_sln = ProjectToSO3(Result[SOneTo(3),StaticArrays.SUnitRange(7,9)]);
			Result = T_SVD(Rx_sln,Ry_sln,Rz_sln,RA,RB,RC,tA,tB,tC)
			Xᵣ = Result[1]; Yᵣ = Result[2]; Zᵣ = Result[3];
			return Xᵣ,Yᵣ,Zᵣ,err
		end
	end
	Xᵣ,Yᵣ,Zᵣ = post_process_Wang(Result);
	return Xᵣ,Yᵣ,Zᵣ,err
end

"""
	AXB_YCZ_Close(A,B,C,conf)

The closed-form part of Wang's solution to the AXB=YCZ problem.

# Argument
- `A`: Vector of Transformation Matrices in the form of Vector{SMatirx{4,4,Float64,16}}
- `B`: Vector of Transformation Matrices in the form of Vector{SMatirx{4,4,Float64,16}}
- `C`: Vector of Transformation Matrices in the form of Vector{SMatirx{4,4,Float64,16}}
- `conf`: Not used in this function, left here for consistency with the other functions.
"""
function AXB_YCZ_Close(A,B,C,conf)
	n = size(A,1)
	RA = [A[i][SOneTo(3),SOneTo(3)] for i = 1:n];
	RB = [B[i][SOneTo(3),SOneTo(3)] for i = 1:n];
	RC = [C[i][SOneTo(3),SOneTo(3)] for i = 1:n];
	tA = [A[i][SOneTo(3),4] for i = 1:n];
	tB = [B[i][SOneTo(3),4] for i = 1:n];
	tC = [C[i][SOneTo(3),4] for i = 1:n];
	Rxᵢₙᵢₜ,Ryᵢₙᵢₜ,Rzᵢₙᵢₜ = solveRS(RA,RB,RC);
	# return T_SVD(Rxᵢₙᵢₜ,Ryᵢₙᵢₜ,Rzᵢₙᵢₜ,RA,RB,RC,tA,tB,tC)
	Rxᵢₙᵢₜ,Txᵢₙᵢₜ,Ryᵢₙᵢₜ,Tyᵢₙᵢₜ,Rzᵢₙᵢₜ,Tzᵢₙᵢₜ = solveBestT(RA,tA,RB,tB,RC,tC,Rxᵢₙᵢₜ,Ryᵢₙᵢₜ,Rzᵢₙᵢₜ)
	return rt2T(Rxᵢₙᵢₜ,Txᵢₙᵢₜ),rt2T(Ryᵢₙᵢₜ,Tyᵢₙᵢₜ),rt2T(Rzᵢₙᵢₜ,Tzᵢₙᵢₜ)
end

#########################################################################
####                  CLOSED-FORM HELPER FUNCTIONS                   ####
#########################################################################
"""
	solveRS(RA,RB,RC)

Solve the rotation matrices of the AXB=YCZ problem in closed-form.

# Argument
- `RA`: Vector of Rotation Matrices in the form of Vector{SMatirx{3,3,Float64,9}}
- `RB`: Vector of Rotation Matrices in the form of Vector{SMatirx{3,3,Float64,9}}
- `RC`: Vector of Rotation Matrices in the form of Vector{SMatirx{3,3,Float64,9}}
"""
function solveRS(RA,RB,RC)
	R₁ = solveR(RA,RB,RC);
	R₂ = solveR(transpose.(RA),RC,RB)
	R₃ = solveR(RC,transpose.(RB),RA)
	# return ProjectToSO3(R₁),ProjectToSO3(R₂),ProjectToSO3(R₃)
	return R₁,R₂,R₃
end

"""
	solveR(A,B,C)

Solve the rotation matrix of the AXB=YCZ problem in closed-form with kronecker product.
"""
function solveR(A,B,C)
	n = size(A,1)
	M = zeros(9n,90)
	for i = 1:n
		M[9i-8:9i,:] = [kron(transpose(B[i]),A[i]) -toMc(tovec(C[i]))]
	end
	U,Σ,V = svd(M'*M);
	R = reshape(V[1:9,90]/norm(V[1:9,90]),(3,3))*√3
	R = SMatrix{3,3,Float64,9}(R);
	# R = reshape(V[1:9,90],(3,3))*2*√3
end

function toMc(vecRc)
	Mc = zeros(9,81);
	for i = 1:9
		Mc[i,(i-1)*9+1:i*9] = vecRc
	end
	return Mc
end

function tovec(R)
	reshape(R,(size(R,1)*size(R,2),1))
end

"""
	solveBestT(RA,TA,RB,TB,RC,TC,Rxᵢₙᵢₜ,Ryᵢₙᵢₜ,Rzᵢₙᵢₜ)

Solve the translation vectors of the AXB=YCZ problem in closed-form with SVD method.

# Argument
- `RA`: Vector of Rotation Matrices in the form of Vector{SMatirx{3,3,Float64,9}}
- `TA`: Vector of Translation Vectors in the form of Vector{SVector{3,Float64}}
- `RB`: Vector of Rotation Matrices in the form of Vector{SMatirx{3,3,Float64,9}}
- `TB`: Vector of Translation Vectors in the form of Vector{SVector{3,Float64}}
- `RC`: Vector of Rotation Matrices in the form of Vector{SMatirx{3,3,Float64,9}}
- `TC`: Vector of Translation Vectors in the form of Vector{SVector{3,Float64}}
- `Rxᵢₙᵢₜ`: Rotation Matrix of the AXB=YCZ problem in the form of SMatrix{3,3,Float64,9}
- `Ryᵢₙᵢₜ`: Rotation Matrix of the AXB=YCZ problem in the form of SMatrix{3,3,Float64,9}
- `Rzᵢₙᵢₜ`: Rotation Matrix of the AXB=YCZ problem in the form of SMatrix{3,3,Float64,9}
"""
function solveBestT(RA,TA,RB,TB,RC,TC,Rxᵢₙᵢₜ,Ryᵢₙᵢₜ,Rzᵢₙᵢₜ)
	Rx = sign(det(Rxᵢₙᵢₜ))*Rxᵢₙᵢₜ
	Ry = sign(det(Ryᵢₙᵢₜ))*Ryᵢₙᵢₜ
	Rz = sign(det(Rzᵢₙᵢₜ))*Rzᵢₙᵢₜ
	Tx,Ty,Tz,err = solveTs(RA,TA,RB,TB,RC,TC,Rx,Ry,Rz)
	return Rx, Tx, Ry, Ty, Rz, Tz
end

function solveTs(RA,TA,RB,TB,RC,TC,Rx,Ry,Rz)
	n = size(RA,1)
	J = zeros(3n,9)
	b = zeros(3n,1)
	for i = 1:n
		J[3i-2:3i,:] = [RA[i] -eye(3) -Ry*RC[i]]
		b[3i-2:3i,1] = Ry*TC[i]-TA[i]-RA[i]*Rx*TB[i]
	end
	t = (J'*J)\(J'*b);
	Tx = t[1:3,1]; Ty = t[4:6,1]; Tz = t[7:9,1];
	Aₚ = [RA[1] TA[1];0.0 0.0 0.0 1.0]
	Bₚ = [RB[1] TB[1];0.0 0.0 0.0 1.0]
	Cₚ = [RC[1] TC[1];0.0 0.0 0.0 1.0]
	Xₚ = [Rx Tx;0.0 0.0 0.0 1.0]
	Yₚ = [Ry Ty;0.0 0.0 0.0 1.0]
	Zₚ = [Rz Tz;0.0 0.0 0.0 1.0]
	# err = local_error(RA,TA,RB,TB,RC,TC,Xₚ,Yₚ,Zₚ)
	err = norm(Aₚ*Xₚ*Bₚ-Yₚ*Cₚ*Zₚ,2)
	# println(err)
	return Tx,Ty,Tz,err
end

############################################################################
####                  ITERATIVE PART HELPER FUNCTIONS                   ####
############################################################################
"""
	get_full_grad_wangᵥ(V,P,conf)

Compute the gradient of the objective function of the AXB=YCZ problem with respect to the whole dataset.

See also [`get_one_grad_wangᵥ`](@ref).
"""
function get_full_grad_wangᵥ(V,P,conf)
	Aₛ = P[1]; Bₛ = P[2]; Cₛ = P[3];
	n = size(Aₛ,1)
	@assert size(Aₛ,1)==size(Bₛ,1)==size(Cₛ,1)
	grad = sum(get_one_grad_wangᵥ.(Aₛ,Bₛ,Cₛ,[V],[conf.μ]))/n
	return grad
end


"""
	get_i_grad_wangᵥ(V,P,i,conf)

Computes the gradient of the Wang objective function for a specific (ith) triple of poses A,B,C.

See also [`get_one_grad_wangᵥ`](@ref)
"""
function get_i_grad_wangᵥ(V,P,i,conf)
	A = P[1][i]; B = P[2][i]; C = P[3][i];
	get_one_grad_wangᵥ(A,B,C,V,conf.μ)
end

"""
	get_one_grad_wangᵥ(A,B,C,V,μ)

Computes the gradient of the Wang objective function for a single triple of poses A,B,C.

# Arguments
- `A::SMatrix{4,4,Float64,16}`
- `B::SMatrix{4,4,Float64,16}`
- `C::SMatrix{4,4,Float64,16}`
- `V::SMatrix{3,12,Float64,16}` the current value of optimization variables
- `μ::Float64` the weights of the translation part of the objective function
"""
function get_one_grad_wangᵥ(A,B,C,V,μ)
	Rx = V[:,SOneTo(3)]; Ry = V[:,StaticArrays.SUnitRange(4,6)]; Rz = V[:,StaticArrays.SUnitRange(7,9)];
	Tx = V[:,10]; Ty = V[:,11]; Tz = V[:,12];
	Ra = A[SOneTo(3),SOneTo(3)]; Rb = B[SOneTo(3),SOneTo(3)]; Rc = C[SOneTo(3),SOneTo(3)];
	Ta = A[SOneTo(3),4]; Tb = B[SOneTo(3),4]; Tc = C[SOneTo(3),4];
	row1 = μ[1]*(Rx*Rb*Rb'-Ra'*Ry*Rc*Rz*Rb')
	      +μ[2]*(Rx*Tb*Tb'+Tx*Tb'+Ra'*Ta*Tb'-Ra'*Ry*Rc*Tz*Tb'-Ra'*Ry*Tc*Tb'-Ra'*Ty*Tb')+2*μ[3]*(Rx*Rx'*Rx-Rx);
	row2 = μ[1]*(Ry*Rc*Rz*Rz'*Rc'-Ra*Rx*Rb*Rz'*Rc')
	      +μ[2]*(Ry*Tc*Tc'+Ry*Rc*Tz*Tz'*Rc'-Ra*Rx*Tb*Tz'*Rc'-Ra*Rx*Tb*Tc'-Ra*Tx*Tz'*Rc'-Ra*Tx*Tc'-Ta*Tz'*Rc'-Ta*Tc'+Ry*Tc*Tz'*Rc'+Ry*Rc*Tz*Tc'+Ty*Tz'*Rc'+Ty*Tc')+2*μ[4]*(Ry*Ry'*Ry-Ry);
	row3 = μ[1]*(Rc'*Ry'*Ry*Rc*Rz-Rc'*Ry'*Ra*Rx*Rb)+2μ[5]*(Rz*Rz'*Rz-Rz)
	row4 = μ[2]*(Tx+Rx*Tb+Ra'*Ta-Ra'*Ry*Rc*Tz-Ra'*Ry*Tc-Ra'*Ty)
	row5 = μ[2]*(Ty-Ra*Rx*Tb-Ra*Tx-Ta+Ry*Rc*Tz+Ry*Tc)
	row6 = μ[2]*(Rc'*Ry'*Ry*Rc*Tz-Rc'*Ry'*Ra*Rx*Tb-Rc'*Ry'*Ra*Tx-Rc'*Ry'*Ta+Rc'*Ry'*Ry*Tc+Rc'*Ry'*Ty)
	grad = hcat(row1,row2,row3,row4,row5,row6)
	return grad
end

"""
	errfunc_wang(V,P,∇ᵥgṼₜ,ηₜ,conf)

Error function for Wang's Algoirthm based on the Frobenius norm of the error matrix.
"""
function errfunc_wang(V,P,∇ᵥgṼₜ,ηₜ,conf)
	X,Y,Z = post_process_Wang(V)
	A = P[1]; B = P[2]; C = P[3];
	return sum(err.(A,B,C,[X],[Y],[Z]))/size(P[1],1)
	# return norm(A*X*B-Y*C*Z,2)
end

"""
	errN(V,P,∇ᵥgṼₜ,ηₜ,conf)

Null error function for testing purposes
"""
function errN(V,P,∇ᵥgṼₜ,ηₜ,conf)
	return 0.0
end

"""
	err_Wang_Full(V,P,∇ᵥgṼₜ,ηₜ,conf)

Complete error function for Wang's Algoirthm. 
Returns the error in the rotation and translation matrices, the norm of the gradient and the error function.
"""
function err_Wang_Full(V,P,∇ᵥgṼₜ,ηₜ,conf)
	Xᵣ,Yᵣ,Zᵣ = post_process_Wang(V)
	A = P[1]; B = P[2]; C = P[3];
	Errx = get_error(Xᵣ,X);
	Erry = get_error(Yᵣ,X);
	Errz = get_error(Zᵣ,Z);
	norme = norm(∇ᵥgṼₜ,2)
	error = sum(err.(A,B,C,[X],[Y],[Z]))/size(P[1],1)
	return [Errx[1] Errx[2] Erry[1] Erry[2] Errz[1] Errz[2] ηₜ norme error]
end

err(A,B,C,X,Y,Z) = norm(A*X*B-Y*C*Z,2)

"""
	get_ηₜ(Ṽₜ,Ṽₜ₁,∇ᵥgṼₜ,∇ᵥgṼₜ₁,conf)

Computes the step size of SVRG-BB for Wang's Algoirthm.
"""
function get_ηₜ(Ṽₜ,Ṽₜ₁,∇ᵥgṼₜ,∇ᵥgṼₜ₁,conf)
	m = conf.m
	ηₜ = 1/m*norm((Ṽₜ-Ṽₜ₁),2)^2/norm((Ṽₜ-Ṽₜ₁)'*(∇ᵥgṼₜ-∇ᵥgṼₜ₁),2)
	return ηₜ
end

"""
	post_process_Wang(Ṽₜ)

Decode the result vector of Wang's Algoirthm into the transformation matrices.
Noted that the rotation matrices are projected to the SO(3) manifold.
"""
function post_process_Wang(Ṽₜ)
	Rx = Ṽₜ[:,SOneTo(3)]; Ry = Ṽₜ[:,StaticArrays.SUnitRange(4,6)]; Rz = Ṽₜ[:,StaticArrays.SUnitRange(7,9)];
	Rx = ProjectToSO3(Rx);Ry = ProjectToSO3(Ry);Rz = ProjectToSO3(Rz);
	Tx = Ṽₜ[:,10]; Ty = Ṽₜ[:,11]; Tz = Ṽₜ[:,12];
	return rt2T(Rx,Tx), rt2T(Ry,Ty), rt2T(Rz,Tz)
end

"""
	iter_solve(V,P,conf)

Solve the problem with SVRG-BB for Wang's Algoirthm (aligned with the method in Wang's paper).

# Arguments
- `V::Vector{SMatrix{4,4,Float64,16}}`: Initial guess of the transformation matrices.
- `P::Vector{Vector{SMatrix{4,4,Float64,16}}}`: The data sets as `{{A},{B},{C}}`
"""
function iter_solve(V,P,conf)
	Vᵢₙᵢₜ = toV(V)
	if conf.err == "NAN"
		errfunc = errN
	elseif conf.err == "LS"
		errfunc = errfunc_wang
	elseif conf.err == "Ext"
		errfunc = conf.errfunc
	end
	Result,Err = SVRG_BB(Vᵢₙᵢₜ,P,get_i_grad_wangᵥ,get_full_grad_wangᵥ,get_ηₜ,errfunc,conf)
	return Result,Err
end

"""
	iter_solve(V,P,conf)

Solve the problem with SVRG for Wang's Algoirthm.

# Arguments
- `V::Vector{SMatrix{4,4,Float64,16}}`: Initial guess of the transformation matrices.
- `P::Vector{Vector{SMatrix{4,4,Float64,16}}}`: The data sets as `{{A},{B},{C}}`
"""
function iter_solve_SVRG(V,P,conf)
	Vᵢₙᵢₜ = toV(V)
	if conf.err == "NAN"
		errfunc = errN
	elseif conf.err == "LS"
		errfunc = errfunc_wang
	elseif conf.err == "Ext"
		errfunc = conf.errfunc
	end
	Result,Err = SVRG(Vᵢₙᵢₜ,P,get_i_grad_wangᵥ,get_full_grad_wangᵥ,errfunc,normfunc,conf)
	return Result,Err
end

"""
	iter_solve(V,P,conf)

Solve the problem with KatyushaX for Wang's Algoirthm.

# Arguments
- `V::Vector{SMatrix{4,4,Float64,16}}`: Initial guess of the transformation matrices.
- `P::Vector{Vector{SMatrix{4,4,Float64,16}}}`: The data sets as `{{A},{B},{C}}`
"""
function iter_solve_Katyusha(V,P,conf)
	Vᵢₙᵢₜ = toV(V)
	if conf.err == "NAN"
		errfunc = errN
	elseif conf.err == "LS"
		errfunc = errfunc_wang
	elseif conf.err == "Ext"
		errfunc = conf.errfunc
	end
	Result,Err = katyushaX(Vᵢₙᵢₜ,P,get_i_grad_wangᵥ,get_full_grad_wangᵥ,errfunc,normfunc,conf)
	return Result,Err
end


function normfunc(V)
	V[1:3,1:3] = ProjectToSO3(V[1:3,1:3])
	V[1:3,4:6] = ProjectToSO3(V[1:3,4:6])
	V[1:3,7:9] = ProjectToSO3(V[1:3,7:9])
	return V
end

function toV(Init)
	return hcat(Init[1][SOneTo(3),SOneTo(3)],Init[2][SOneTo(3),SOneTo(3)],Init[3][SOneTo(3),SOneTo(3)],Init[1][SOneTo(3),4],Init[2][SOneTo(3),4],Init[3][SOneTo(3),4])
end

function toSM(M)
	return @SMatrix [M[1,1] M[1,2] M[1,3] M[1,4]; M[2,1] M[2,2] M[2,3] M[2,4]; M[3,1] M[3,2] M[3,3] M[3,4]; 0.0 0.0 0.0 1.0]
end