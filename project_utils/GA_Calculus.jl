#########################################################################
####                      FUNCTIONS FOR AXB=YCZ                      ####
#########################################################################

## Iteration Solve
function iter_solve_GCᵥ(V,P,conf)
	scale = conf.scale
	V = toVᵥ(scale_motorS.(V,scale))
	P = [(scale_motorS.(P[i],scale)) for i∈eachindex(P)]
	if conf.err == "NAN"
		errfunc = Errₙ
	elseif conf.err == "LS"
		errfunc = Errᵥ₂
	elseif conf.err == "Ext"
		errfunc = conf.errfunc;
	end
	if conf.reg
		V,err = SVRG(V,P,get_i_grad_reg_AXBYCZᵥ₂,get_full_grad_regᵥ₂, errfunc,normfunc2ᵥ, conf)
	else
		V,err = SVRG(V,P,get_i_grad_AXBYCZᵥ₂,get_full_gradᵥ₂, errfunc ,normfuncᵥ₂, conf)
	end
	V = motorNormS.(toXYZS(V))
	V = scale_motorS.(V,1/scale)
	return V,err
end

function iter_solve_GC_katyushaᵥ(V,P,conf)
	scale = conf.scale
	V = toVᵥ(scale_motorS.(V,scale))
	P = [(scale_motorS.(P[i],scale)) for i∈eachindex(P)]
	if conf.err == "NAN"
		errfunc = Errₙ
	elseif conf.err == "LS"
		errfunc = Errᵥ₂
	elseif conf.err == "Ext"
		errfunc = conf.errfunc;
	end
	if conf.reg
		V,err = katyushaX(V,P,get_i_grad_reg_AXBYCZᵥ₂,get_full_grad_regᵥ₂, errfunc,normfuncᵥ₂, conf)
	else
		V,err = katyushaX(V,P,get_i_grad_AXBYCZᵥ₂,get_full_gradᵥ₂, errfunc ,normfuncᵥ₂, conf)
	end
	V = motorNormS.(toXYZS(V))
	V = scale_motorS.(V,1/scale)
	return V,err
end

function iter_solve_BB_GCᵥ(V,P,conf)
	scale = conf.scale
	V = toVᵥ(motor2vectorS.(scale_motor.(V,scale)))
	P = [motor2vectorS.(scale_motor.(P[i],scale)) for i∈eachindex(P)]
	if conf.err == "NAN"
		errfunc = Errₙ
	elseif conf.err == "LS"
		errfunc = Errᵥ₂
	elseif conf.err == "Ext"
		errfunc = conf.errfunc;
	end
	V,err = SVRG_BB(V,P,get_i_grad_reg_AXBYCZᵥ₂, get_full_grad_regᵥ₂, get_ηₜᵥ₂, errfunc, conf)
	V = motorNorm.(toXYZᵥ₂(V))
	V = scale_motor.(V,1/scale)
	return V,err
end

## Get Gradient
function get_i_grad_AXBYCZᵥ₂(V,P,i,conf)
    A = P[1][i]; B = P[2][i]; C = P[3][i];
	return get_one_gradᵥ₂(A,B,C,V)
end

function get_one_gradᵥ₂(A,B,C,V)
	E = E_ABCXYZᵥ₂(A,B,C,V)
	Y = V[StaticArrays.SUnitRange(9,16)];
	Z = V[StaticArrays.SUnitRange(17,24)];
	gradx = ecmul(ereversion(B),E,ereversion(A));
	grady = -ecmul(ereversion(ecmul(Z,C)),E);
	gradz = -ecmul(E,ereversion(ecmul(C,Y)));
	return vcat(gradx,grady,gradz)
end

function get_one_grad_regᵥ₂(A,B,C,V,μ)
	return get_one_gradᵥ₂(A,B,C,V).+get_regulatorᵥ₂(V,μ)
end

function get_i_grad_reg_AXBYCZᵥ₂(V,P,i,conf)
	return get_i_grad_AXBYCZᵥ₂(V,P,i,conf).+get_regulatorᵥ₂(V,conf.μ)
end

function get_regulatorᵥ₂(V,μ)
	X = V[StaticArrays.SUnitRange(1,8)];
	Y = V[StaticArrays.SUnitRange(9,16)];
	Z = V[StaticArrays.SUnitRange(17,24)];
	E₁ = E_normᵥ(X); E₂ = E_normᵥ(Y); E₃ = E_normᵥ(Z)
	gradx = μ.*ecmul((ereversion(E₁).+E₁),X);
	grady = μ.*ecmul((ereversion(E₂).+E₂),Y);
	gradz = μ.*ecmul((ereversion(E₃).+E₃),Z);
	return vcat(gradx,grady,gradz)
end

function get_full_gradᵥ₂(V,P,conf)
	grad = sum(get_one_gradᵥ₂.(P[1],P[2],P[3],[V]))/size(P[1],1)
	return grad
end

function get_full_grad_regᵥ₂(V,P,conf)
	μ = conf.μ
	grad = sum(get_one_grad_regᵥ₂.(P[1],P[2],P[3],[V],[μ]))/size(P[1],1)
	return grad
end

## Variable Step size
function get_ηₜᵥ₂(Ṽₜ::SVector,Ṽₜ₁::SVector,∇ᵥgṼₜ::SVector,∇ᵥgṼₜ₁::SVector,conf)
	m = conf.m
	a1 = Ṽₜ.-Ṽₜ₁; b1 = ∇ᵥgṼₜ.-∇ᵥgṼₜ₁;
	n1 = norm(a1);
	n2 = sum(a1.*b1)
	ηₜ = 1/m*n1^2/n2;
	return ηₜ
end

## Err FUNCTIONS

# function Errᵥ₂(V,P,∇ᵥgṼₜ,ηₜ,conf)
# 	A = P[1][1]; B = P[2][1]; C = P[3][1];
# 	return norm(E_ABCXYZᵥ₂(A,B,C,V))
# end

function Errᵥ₂(V,P,∇ᵥgṼₜ,ηₜ,conf)
	Vₙ = normfuncᵥ₂(V)
	return norm(E_ABCXYZᵥ₂.(P[1],P[2],P[3],[Vₙ]))
end
function Errₙ(V,P,∇ᵥgṼₜ,ηₜ,conf)
	return 0.0
end
# function Err_Fullᵥ₂(V,P,∇ᵥgṼₜ,ηₜ,conf)
# 	  Result = scale_motor.(toXYZ(V),1/conf.scale);
#     Errx = get_error(Result[1],Xₜ);
#     Erry = get_error(Result[2],Yₜ);
#     Errz = get_error(Result[3],Zₜ);
#     norme = norm(∇ᵥgṼₜ)
#     err = Errᵥ₂(V,P,∇ᵥgṼₜ,ηₜ,conf)
#     return [Errx[1] Errx[2] Erry[1] Erry[2] Errz[1] Errz[2] ηₜ norme err];
# end

E_ABCXYZᵥ₂(A,B,C,V) = ecmul(B,V[SOneTo(8)],A).-ecmul(V[StaticArrays.SUnitRange(17,24)],C,V[StaticArrays.SUnitRange(9,16)])

function E_normᵥ(X)
	v1 = @SVector [1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
	return (ecmul(X,ereversion(X)).-v1);
end
err_normᵥ(X) = norm(E_normᵥ(X));

## Other Functions
function normfuncᵥ₂(V)
	Vₘ = [motorNormS(V[StaticArrays.SUnitRange(8*(i-1)+1,8*i)]) for i = 1:3];
	return toVᵥ(Vₘ)
end

toVᵥ(Init) = vcat(Init[1],Init[2],Init[3])

function toXYZS(V)
	Result = [(V[StaticArrays.SUnitRange(8*(i-1)+1,8*i)]) for i = 1:3];
	return Result
end

# function toXYZᵥ₂(V)
#     Result = [vector2motor(V[StaticArrays.SUnitRange(8*(i-1)+1,8*i)]) for i = 1:3]
#     return Result
# end

@inline function normfunc2ᵥ(V)
	return V
end