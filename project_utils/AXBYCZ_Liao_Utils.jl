# This function implements the AXB=YCZ algorithm based on the paper "Simultaneous 
# Hand-Eye, Tool-Flange, and Robot-Robot Calibration for Comanipulation by Solving the 
# AXB = YCZ Problem" by Liao Wu, (2016). A matlab version published by Wu is available
# at https://github.com/drliaowu/AXBYCZ
"""
	AXBYCZ_Liao(A,B,C,Aₘ,Bₘ,Cₘ,conf)

A modified version of Liao Wu 2016's paper on AXB=YCZ.
A,B,C is in matrix form. Aₘ,Bₘ,Cₘ is in motor's form, which gives a precies sign of A, B, C in Quaternion
**This version cannot be used in physical experiment.**
**The sign fixing process of B is removed for performance consideration**
"""
function AXBYCZ_Liao(A::T1,B::T1,C::T1,Aₘ::T2,Bₘ::T2,Cₘ::T2,conf=()) where {T1<:Vector{SMatrix{4,4,Float64,16}},T2<:Vector{SVector{8,Float64}}}
	n = size(A,1);
	RA = [@view A[i][SOneTo(3),SOneTo(3)] for i = 1:n];
	RB = [@view B[i][SOneTo(3),SOneTo(3)] for i = 1:n];
	RC = [@view C[i][SOneTo(3),SOneTo(3)] for i = 1:n];
	tA = [@view A[i][SOneTo(3),4] for i = 1:n];
	tB = [@view B[i][SOneTo(3),4] for i = 1:n];
	tC = [@view C[i][SOneTo(3),4] for i = 1:n];	
	RX01,RY01,RZ01 = liao_close(motorS2rotorG3S.(Aₘ),motorS2rotorG3S.(Bₘ),motorS2rotorG3S.(Cₘ));
	RX_sln,RY_sln,RZ_sln = liao_iter(RA,RB,RC,RX01,RY01,RZ01,conf);
	return T_SVD(RX_sln,RY_sln,RZ_sln, RA,RB,RC,tA,tB,tC)
end


"""
	liao_close(Aₘ,Bₘ,Cₘ)

Solve the rotational component of X,Y,Z in closed form.

# Argument
- `Aₘ`: Vector of rotors in the form of Vector{SVector{4,Float64}}
- `Bₘ`: Vector of rotors in the form of Vector{SVector{4,Float64}}
- `Cₘ`: Vector of rotors in the form of Vector{SVector{4,Float64}}

Noted that A,B,C in Aₘ,Bₘ,Cₘ should satisfies AXB=YCZ. 
Any combination that does not satisfy this condition (AXB=-YCZ) will result in a wrong solution.
This is the signed problem mentioned in Wu's and this paper.
Sign check should be carried out before the data is injected.

# return
- `RX01`: Rotation matrix of X
- `RY01`: Rotation matrix of Y
- `RZ01`: Rotation matrix of Z
"""
function liao_close(Aₘ::T,Bₘ::T,Cₘ::T) where {T<:Vector{SVector{4,Float64}}}
	M_data = size(Aₘ,1);
	WAB = zeros(4*M_data,4);
	WC = zeros(4*M_data,16);
	for i = 1:M_data
		WAB[4*i-3:4*i,1:4] = RQ2M(rotorG32quat(Bₘ[i]))*LQ2M(rotorG32quat(Aₘ[i]));
    	WC[4*i-3:4*i,1:16] = GetWC(rotorG32quat(Cₘ[i]));
	end
	WABC = [WAB -WC];
	F = eigen(WABC'*WABC)
	V1 = F.vectors; Dtemp = F.values;
	
	RX01, RY01, RZ01 = qXYZ2RXYZ(V1);
	return RX01,RY01,RZ01
end

"""
	liao_iter(RA,RB,RC,RX_init,RY_init,RZ_init,conf)

Iterative part of Wu's algorithm.

# Argument
- `RA`: Vector of rotation matrix in the form of Vector{SMatrix{3,3,Float64,9}}
- `RB`: Vector of rotation matrix in the form of Vector{SMatrix{3,3,Float64,9}}
- `RC`: Vector of rotation matrix in the form of Vector{SMatrix{3,3,Float64,9}}
- `RX_init`: Initial guess of rotation matrix of X
- `RY_init`: Initial guess of rotation matrix of Y
- `RZ_init`: Initial guess of rotation matrix of Z
- `conf`: Configuration of the algorithm. 
	- `conf.max_iter`: Maximum iteration number. Must be provied.
	- `conf.err`: Error type. Must be provied.
		- "Ext": External error. The error is calculated by the error function provided by the user.
	- `conf.errfunc`: Error function. Should be provied if conf.err is "Ext".
	- `conf.errdim`: Error dimension. Should be provied if conf.err is "Ext".
"""
function liao_iter(RA,RB,RC,RX_init,RY_init,RZ_init,conf)
	if haskey(conf,:err)
		if conf.err == "Ext"
			errfunc = conf.errfunc;
			err = zeros(conf.max_iter,conf.errdim)
		else
			err = 0;
		end
	else 
		err = 0;
	end
	iter = 0;
	while iter < conf.max_iter
	    F=GetF(RA,RB,RC,RX_init,RY_init,RZ_init);
	    c=Getc(RA,RB,RC,RX_init,RY_init,RZ_init);
	    r = SVector{9,Float64}(F\c);
		
	    # RX_init = angvec2r(norm(r[1:3]),r[1:3]/norm(r[1:3]))*RX_init;
	    # RY_init = angvec2r(norm(r[4:6]),r[4:6]/norm(r[4:6]))*RY_init;
	    # RZ_init = angvec2r(norm(r[7:9]),r[7:9]/norm(r[7:9]))*RZ_init;
		RX_init = RotationVec(r[SOneTo(3)]...) * RX_init;
		RY_init = RotationVec(r[StaticArrays.SUnitRange(4,6)]...)*RY_init;
		RZ_init = RotationVec(r[StaticArrays.SUnitRange(7,9)]...)*RZ_init;
	    iter = iter+1;
		if haskey(conf,:err)
			if conf.err == "Ext"
				err[iter,:] = errfunc([RX_init,RY_init,RZ_init]);
			end
		end
		# if norm(r)<1e-15	# Break In case numerical failure
		# 	break
		# end
	end
	RX_sln = SMatrix{3,3,Float64,9}(RX_init);
	RY_sln = SMatrix{3,3,Float64,9}(RY_init);
	RZ_sln = SMatrix{3,3,Float64,9}(RZ_init);
	return RX_sln,RY_sln,RZ_sln,err
end

RQ2M = Qₗ
LQ2M = Qᵣ

function GetWC(qC)
	qC0 = qC[1]; qC1 = qC[2]; qC2 = qC[3]; qC3 = qC[4]
	WC = @SMatrix [qC0     -qC1    -qC2    -qC3    -qC1    -qC0    qC3     -qC2    -qC2    -qC3    -qC0    qC1     -qC3    qC2     -qC1    -qC0;
     qC1     qC0     -qC3    qC2     qC0     -qC1    -qC2    -qC3    qC3     -qC2    qC1     qC0     -qC2    -qC3    -qC0    qC1;
     qC2     qC3     qC0     -qC1    -qC3    qC2     -qC1    -qC0    qC0     -qC1    -qC2    -qC3    qC1     qC0     -qC3    qC2;
     qC3     -qC2    qC1     qC0     qC2     qC3     qC0     -qC1    -qC1    -qC0    qC3     -qC2    qC0     -qC1    -qC2    -qC3 ];
end

function qXYZ2RXYZ(V)

	qX = V[1:4,1]/norm(V[1:4,1]);
	
	qYZ = V[5:20,1]/norm(V[1:4,1]);
	
	RX0 = Quat2rm(qX)
	
	qYZ = qYZ*sign(qYZ[1]);
	
	SignY0 = 1; SignY1 = sign(qYZ[5]); SignY2 = sign(qYZ[9]); SignY3 = sign(qYZ[13]);
	qY0 = norm(qYZ[collect(1:4)])*SignY0; qY1 = norm(qYZ[collect(5:8)])*SignY1; 
	qY2 = norm(qYZ[collect(9:12)])*SignY2; qY3 = norm(qYZ[collect(13:16)])*SignY3;
	qY = [qY0;qY1;qY2;qY3]; qY = qY/norm(qY);
	
	SignZ0 = 1; SignZ1 = sign(qYZ[2]); SignZ2 = sign(qYZ[3]); SignZ3 = sign(qYZ[4]);
	qZ0 = norm(qYZ[[1,5,9,13]])*SignZ0; qZ1 = norm(qYZ[[2,6,10,14]])*SignZ1; 
	qZ2 = norm(qYZ[[3,7,11,15]])*SignZ2; qZ3 = norm(qYZ[[4,8,12,16]])*SignZ3;
	qZ = [qZ0;qZ1;qZ2;qZ3]; 
	qZ = qZ/norm(qZ);
	
	RY0 = Quat2rm(qY)
	RZ0 = Quat2rm(qZ)

	return RX0, RY0, RZ0
end
angvec2r(θ,r) = MatrixExp3(mcross(r)*θ)

function bitget(n,i)
	D = digits(n, base = 2)
	if i > length(D)
		return 0
	else
		return D[i]
	end
end

function GetF(RA_noise,RB_noise, RC_noise, RX_init, RY_init,RZ_init)
	M = size(RA_noise,1);
	F = zeros(9*M,9);
	for i = 1:M
		RA = RA_noise[i];
		RB = RB_noise[i];
		RC = RC_noise[i];
	    RXB = RX_init*RB; 
	    RYCZ = RY_init*RC*RZ_init; 
	    
	    F[9*i-8:9*i-6,:] = hcat(-RA*mcross(RXB[:,1]),mcross(RYCZ[:,1]),RY_init*RC*mcross(RZ_init[:,1]));
	    F[9*i-5:9*i-3,:] = hcat(-RA*mcross(RXB[:,2]),mcross(RYCZ[:,2]),RY_init*RC*mcross(RZ_init[:,2]));
	    F[9*i-2:9*i,:] = hcat(-RA*mcross(RXB[:,3]),mcross(RYCZ[:,3]),RY_init*RC*mcross(RZ_init[:,3]));
	end
	return F
end

function Getc(RA_noise,RB_noise,RC_noise,RX_init,RY_init,RZ_init)
	M=size(RA_noise,1);
	c=zeros(9*M,1);
	for i=1:M
	    RA=RA_noise[i]; 
	    RB=RB_noise[i]; 
	    RC=RC_noise[i];
	    RAXBYCZ=-RA*RX_init*RB + RY_init*RC*RZ_init;
	    c[9*i-8:9*i]=reshape(RAXBYCZ,(9,1));
	end
	return c
end

