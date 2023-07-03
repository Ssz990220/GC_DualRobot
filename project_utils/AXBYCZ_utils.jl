"""
	get_ABC(n::Int,X,Y,Z,unit="mm")

Generate a set of A,B,Cs in matrix form.
**Motions are not random!!**
When unit is "mm", last column should in order of 1 or larger.
When unit is "m", last column should be in order of 0.001 and no more than 10(Unless magnitude of your problem is beyond 10m)
"""
function get_ABC(n::Int,X,Y,Z,unit="mm")
	# X,Y,Z should be given in meter
	(unit == "mm") ? scale = 1 : scale = 1/1000
	Aₛ = Array{Array{Float64,2}}(undef,n);
	Bₛ = Array{Array{Float64,2}}(undef,n);
	Cₛ = Array{Array{Float64,2}}(undef,n);
	@inbounds for i = 1:n
		xyz = ([900. 0 500]+(rand(1,3).-0.5*200))*scale;
		rots = tfrotx((rand()-0.5)*π/3)*tfrotz((rand()-0.5)*π/6)*tfroty((rand()-0.5)*π/8)
		xyzb = ((rand(1,3).-0.5)*50)*scale;
		rotb = tfrotx((rand()-0.5)*π/36)*tfrotz((rand()-0.5)*π/36)*tfroty((rand()-0.5)*π/36)
		Mid = tfxyz(xyz)*rots
		Aₛ[i] = Mid*(@SMatrix [0. 0 1 0;1 0 0 0;0 1 0 0;0 0 0 1])/X;
		Bₛ[i] = tfxyz(xyzb)*rotb;
		Cₛ[i] = Y\((Aₛ[i]*X*Bₛ[i])/Z)
	end
	return Aₛ,Bₛ,Cₛ
end

"""
	get_ABCSₘ(n::Int,Xₘ,Yₘ,Zₘ,unit="mm")

Generate a set of A,B,Cs in MotorS form.
**Motions are random!!**
When unit is "mm", dual part should in order of 1 or larger.
When unit is "m", dual part should be in order of 0.001 and no more than 10(Unless magnitude of your problem is beyond 10m)
"""

function get_ABCSₘ(n::Int,Xₘ,Yₘ,Zₘ,unit="mm")
	Q₁ =  [SVector{6,Float64}(rand(Uniform(-π,π),1,6)) for i = 1:n]
	Q₂ =  [SVector{6,Float64}(rand(Uniform(-π,π),1,6)) for i = 1:n]
	Aₘₛ = UR10.fkineGAS.(Q₁,unit);
	Cₘₛ = UR10.fkineGAS.(Q₂,unit);
	Bₘₛ = ecmul.([Zₘ],ecmul.(Cₘₛ,ecmul.([Yₘ],ecmul.(ereversion.(Aₘₛ),[ereversion(Xₘ)]))))
	return Aₘₛ,Bₘₛ,Cₘₛ,Q₁,Q₂
end

"""
	get_ABC_noiseSₘ(n::Int,Xₘ::SVector,Yₘ::SVector,Zₘ::SVector,n_angle,n_dis,unit="mm", fix=false, dist = "uniform")

Generate a set of A,B,Cs in MotorS form. Configurations for A, Cs are also returned.

check `add_noiseSₘ` for details on noise addition.
"""
function get_ABC_noiseSₘ(n::Int,Xₘ::SVector,Yₘ::SVector,Zₘ::SVector,n_angle,n_dis,unit="mm", fix=false, dist = "uniform")
	Aₘₛ, Bₘₛ, Cₘₛ,Q₁,Q₂ = get_ABCSₘ(n,Xₘ,Yₘ,Zₘ,unit);
	Aₙₘₛ = add_noiseSₘ.(Aₘₛ,n_angle, n_dis,fix,dist)
	Bₙₘₛ = add_noiseSₘ.(Bₘₛ,n_angle, n_dis,fix,dist)
	Cₙₘₛ = add_noiseSₘ.(Cₘₛ,n_angle, n_dis,fix,dist)
	return Aₙₘₛ,Bₙₘₛ,Cₙₘₛ,Q₁,Q₂
end

"""
	match_signS(True_Value,Init_Value)

Match `Init_Value`'s sign to `True_Value`'s sign.
The sign of largest absolute value in rotation part is matched.
"""
function match_signS(True_Value::T,Init_Value::T) where {T<:SVector{8,Float64}}
	Xᵥ = motorS2rotorS(True_Value); Xᵢᵥ = motorS2rotorS(Init_Value);
	index = argmax(abs.(Xᵥ));
	s = sign(Xᵥ[index])*sign(Xᵢᵥ[index]);
	Init_Valueₛ = Init_Value*s;
	return Init_Valueₛ
end

function match_signS(Xᵥ::T,Xᵢᵥ::T) where {T<:SVector{4,Float64}}
	index = argmax(abs.(Xᵥ));
	s = sign(Xᵥ[index])*sign(Xᵢᵥ[index]);
	Xᵢᵥₛ = Xᵢᵥ*s;
	return Xᵢᵥₛ
end

"""
	get_error_XYZ(X,Y,Z,Ans)

Evaluate error of X,Y,Z in angle and displacement
"""
function get_error_XYZ(X,Y,Z,Ans)
	return get_error_XYZ(X,Y,Z,Ans[1],Ans[2],Ans[3])
end

function get_error_XYZ(Result,Ans)
	return get_error_XYZ(Result[1],Result[2], Result[3], Ans[1], Ans[2], Ans[3])
end

function get_error_XYZ(X,Y,Z,A,B,C)
	error = zeros(6)
	xerr = get_error(A,X);
	yerr = get_error(B,Y);
	zerr = get_error(C,Z);
	error[1] = xerr[1]; error[2] = norm(xerr[2]);
	error[3] = yerr[1]; error[4] = norm(yerr[2]);
	error[5] = zerr[1]; error[6] = norm(zerr[2]);
	return error
end


##########################################
##        Close Form SVD solve T        ##
##########################################


"""
	T_SVD(RX_sln,RY_sln,RZ_sln, RA,RB,RC,tA,tB,tC)
Solve Tranlsational Part of X,Y,Z by solving and LSQ problem with SVD.
RX_sln, RY_sln, RZ_sln provides rotational part of X,Y,Z
RA,RB,RC,tA,tB,tC come from the datasets.
"""
function T_SVD(RX_sln,RY_sln,RZ_sln, RA,RB,RC,tA,tB,tC)
	J = GetJ(RA, RY_sln, RC);
	p = Getp(RA, RX_sln, RY_sln, tA, tB, tC);

	t_sln=J\p; tX_sln=t_sln[1:3]; tY_sln=t_sln[4:6]; tZ_sln=t_sln[7:9];

	
	Xₛ = SMatrix{4,4,Float64,16}([RX_sln tX_sln;
    	0.0 0.0 0.0 1.0]);
	Yₛ = SMatrix{4,4,Float64,16}([RY_sln tY_sln;
	    0.0 0.0 0.0 1.0]);
	Zₛ = SMatrix{4,4,Float64,16}([RZ_sln tZ_sln;
	    0.0 0.0 0.0 1.0]);
	return Xₛ,Yₛ,Zₛ
end
"""
	GetJ(RA_noise,RY_sln, RC_noise)
Helper function of T_SVD
"""
function GetJ(RA_noise, RY_sln, RC_noise)
	M = size(RA_noise,1);
	J = zeros(3*M,9);
	
	for i=1:M
	    J[3*i-2:3*i,:]=[RA_noise[i] -eye(3) -RY_sln*RC_noise[i]];
	end
	return J
end
"""
	Getp(RA_noise, RX_sln, RY_sln, tA_noise, tB_noise, tC_noise)
Helper function of T_SVD
"""
function Getp(RA_noise, RX_sln, RY_sln, tA_noise, tB_noise, tC_noise)
	
	M = size(RA_noise,1);
	p = zeros(3*M,1);
	
	for i=1:M
	    p[3*i-2:3*i,1]=-tA_noise[i] - RA_noise[i]*RX_sln*tB_noise[i] + RY_sln*tC_noise[i];
	end
	return p
end

function T_EKF(RX_sln,RY_sln,RZ_sln, RA,RB,RC,tA,tB,tC)
	
	J = GetJ(RA, RY_sln, RC);
	p = Getp(RA, RX_sln, RY_sln, tA, tB, tC);

	t_sln=EKF_SVD(J,p); tX_sln=t_sln[1:3]; tY_sln=t_sln[4:6]; tZ_sln=t_sln[7:9];

	
	Xₛ = SMatrix{4,4,Float64,16}([RX_sln tX_sln;
    	0.0 0.0 0.0 1.0]);
	Yₛ = SMatrix{4,4,Float64,16}([RY_sln tY_sln;
	    0.0 0.0 0.0 1.0]);
	Zₛ = SMatrix{4,4,Float64,16}([RZ_sln tZ_sln;
	    0.0 0.0 0.0 1.0]);
	return Xₛ,Yₛ,Zₛ
end

function EKF_SVD(A,b)
	r1 = 1e-3
	r2 = 5
	@show "EKF" r1 r2
	neqPerSample = 3;
	n_GeoError = size(A,2);
	R=Sdiag((@SVector ones(neqPerSample))*r1);
	Pk=Sdiag((@SVector ones(n_GeoError))*r2);
	Xk=@SVector zeros(n_GeoError);
	for i=1:size(b,1)÷neqPerSample
		Jk=SMatrix{neqPerSample,n_GeoError,Float64,neqPerSample*n_GeoError}(A[(neqPerSample*i-neqPerSample+1):neqPerSample*i,:]);
		Yk=SVector{neqPerSample,Float64}(b[(neqPerSample*i-neqPerSample+1):neqPerSample*i,:]);
		Kk=Pk*Jk'/(Jk*Pk*Jk'+R);
		Xk=Xk+Kk*(Yk-Jk*Xk);
		Pk=(eye(n_GeoError)-Kk*Jk)*Pk;
	end
	return Xk
end