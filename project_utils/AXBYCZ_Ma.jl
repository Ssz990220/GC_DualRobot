# This function implements the AXB=YCZ algorithm based on the paper "Probabilistic approaches 
# to the AXB= YCZ calibration problem in multi-robot systems" by Ma et al., (2018). A matlab 
# version published by Sipu Ruan (Second author of their paper) is available at https://github.com/ruansp/axbycz_calibration
"""
	axbyczProb1(A1,B1,C1,A2,B2,C3,opt,nstd1,nstd2)
This function implements the Prob1 method in the paper
Prerequisites on the input

  A1 is constant with B1 and C1 free

  C2 is constant with A2 adn B2 free

**Authors**: Qianli Ma, qianli.ma622@gmail.com; 
         Zachariah Goh, zach_goh@yahoo.com
**Modifications**: Sipu Ruan, ruansp@jhu.edu
"""
function axbyczProb1(A1::T,B1::T,C1::T,A2::T,B2::T,C2::T,opt=0,nstd1=0.0,nstd2=0.0) where {T<:Vector{SMatrix{4,4,Float64,16}}}
	@assert size(A1,1)==size(B1,1)==size(C1,1)
	@assert size(A2,1)==size(B2,1)==size(C2,1)
	A1 = A1[1]; C2 = C2[1];

	## Solve Z
	# A1 fixed, B1 and C1 free
	Zg, _, MeanC1, MeanB1, _, _ = batchSolveXY(C1, B1, opt, nstd1, nstd2);

	Z_index = det.(Zg) .> 0; s_Z = sum(Z_index);
	Z = Zg[Z_index];

	## Solve X
	A2ᵢₙᵥ = [eye(4)]./A2; B2ᵢₙᵥ = [eye(4)]./B2;
	Xg, _, MeanA2, _, _, _ = batchSolveXY(A2, B2ᵢₙᵥ, opt, nstd1, nstd2);

	_, _, _, MeanB2, _, _ = batchSolveXY(A2ᵢₙᵥ, B2, opt, nstd1, nstd2);

	X_index = det.(Xg) .> 0; s_X = sum(X_index);
	X = Xg[X_index];

	## Solve Y
	Y = [(A1*x*MeanB1/z)/MeanC1 for x∈X for z∈Z];
	Y2 = [(MeanA2*x*MeanB2) / z / C2 for x∈X for z∈Z];
	Y = vcat(Y, Y2);
	s_Y = size(Y,1);

	## Find Optimal X,Y,Z
	w = 1.5
	counter = 1;
	cost = zeros(s_X*s_Y*s_Z);
	@inbounds for (x,y,z) ∈ Iterators.product(X,Y,Z)
		left1 = A1*x*MeanB1; right1 = y*MeanC1*z;
		diff1 = roterror(left1,right1) + w * tranerror(left1, right1);
		left2 = MeanA2 * x * MeanB2; right2 = y * C2 * z;
		diff2 = roterror(left2, right2) + w * tranerror(left2, right2);
			
		cost[counter] = norm(diff1) + norm(diff2);
		counter += 1;
	end

	## Decode index
	I1 = argmin(cost);
	xid = mod(I1,s_X) == 0 ? s_X : mod(I1,s_X);
	yid = Int(ceil(mod(I1,(s_Y*s_X))/s_X)) == 0 ? s_Y : Int(ceil(mod(I1,(s_Y*s_X))/s_X));
	zid = Int(ceil(I1/(s_Y*s_X)));
	
	return X[xid], Y[yid], Z[zid]
end
"""
	axbyczProb3(N,n,A,B,C,Xinit,Yinit,Zinit)
This function implements the Iterative method in the paper

N: number of Pairs of data

n: total number of each pair of Data

## Data format
In each *pair* of data, n/2 of it holds A static, the other n/2 holds C static.

data is stored in the following fashion.

A = [A11(static),A12,A21(static),A22,...,AN1(static),AN2]
B = [B11,B12,B21,B22,...,BN1,BN2]
C = [C11,C12(static),C21,C22(static),...,CN1,CN2(static)]

**Authors**: Qianli Ma, qianli.ma622@gmail.com; 
         Zachariah Goh, zach_goh@yahoo.com
**Modifications**: Sipu Ruan, ruansp@jhu.edu
"""
function axbyczProb3(N,n,A,B,C,Xinit,Yinit,Zinit)
	@assert size(A,1) == size(B,1) == size(C,1) == N*n;
	X_cal = Xu = Xinit; Y_cal = Yu = Yinit; Z_cal = Zu = Zinit;

	xi = @SVector ones(18);
	max_num = 500;
	tol = 1e-15;
	num = 1;
	
	A1s,B1s,C1s,A2s,B2s,C2s = dataSpliter(N,n,A,B,C);
	
	MSA1 = meanCov.(A1s); 
	A1ms = [MSA1[i][1] for i ∈ eachindex(MSA1)]
	SigA1s = [MSA1[i][2] for i ∈ eachindex(MSA1)]
	MSB1 = meanCov.(B1s);
	B1ms = [MSB1[i][1] for i ∈ eachindex(MSB1)]
	SigB1s = [MSB1[i][2] for i ∈ eachindex(MSB1)]
	MSC1 = meanCov.(C1s);
	C1ms = [MSC1[i][1] for i ∈ eachindex(MSC1)]
	SigC1s = [MSC1[i][2] for i ∈ eachindex(MSC1)]
	
	MSA2 = meanCov.(A2s); 
	A2ms = [MSA2[i][1] for i ∈ eachindex(MSA2)]
	SigA2s = [MSA2[i][2] for i ∈ eachindex(MSA2)]
	MSB2 = meanCov.(B2s);
	B2ms = [MSB2[i][1] for i ∈ eachindex(MSB2)]
	SigB2s = [MSB2[i][2] for i ∈ eachindex(MSB2)]
	MSC2 = meanCov.(C2s);
	C2ms = [MSC2[i][1] for i ∈ eachindex(MSC2)]
	SigC2s = [MSC2[i][2] for i ∈ eachindex(MSC2)]
	
	B2is = [invT.(B2s[i]) for i ∈ eachindex(B2s)];
	MSB2i = meanCov.(B2is);
	B2ims = [MSB2i[i][1] for i ∈ eachindex(MSB2i)]
	SigB2is = [MSB2i[i][2] for i ∈ eachindex(MSB2i)]
	
	
	diff = 1;
	while norm(xi) >= tol && diff > tol && num <= max_num
		Mb1 = MbMat_1.(A1ms,[Xu],B1ms,[Yu],C1ms,[Zu],SigB1s,SigC1s);
		MM1 = [Mb1[i][1] for i ∈ eachindex(Mb1)];
		bb1 = [Mb1[i][2] for i ∈ eachindex(Mb1)];
		Mb2 = MbMat_2.(C2ms,[Zu],B2ims,[invT(Yu)],A2ms,[Xu],SigB2s,SigA2s,B2ms)
		MM2 = [Mb2[i][1] for i ∈ eachindex(Mb2)];
		bb2 = [Mb2[i][2] for i ∈ eachindex(Mb2)];
		
		M = vcat(MM1...,MM2...);
		b = vcat(bb1...,bb2...);
		
		xi = M\b;
		
		w_X = SVector{3,Float64}(xi[1:3]);v_X = SVector{3,Float64}(xi[4:6]);
		w_Y = SVector{3,Float64}(xi[7:9]);v_Y = SVector{3,Float64}(xi[10:12]);
		w_Z = SVector{3,Float64}(xi[13:15]);v_Z = SVector{3,Float64}(xi[16:18]);

		X̂ = vcat(hcat(mcross(w_X),v_X),(@SMatrix zeros(1,4)));
		Ŷ = vcat(hcat(mcross(w_Y),v_Y),(@SMatrix zeros(1,4)));
		Ẑ = vcat(hcat(mcross(w_Z),v_Z),(@SMatrix zeros(1,4)));

		X_cal = Xu * MatrixExp6(X̂);
		Y_cal = Yu * MatrixExp6(Ŷ);
		Z_cal = Zu * MatrixExp6(Ẑ);

		Xu = X_cal;
		Yu = Y_cal;
		Zu = Z_cal;
		num = num + 1;
		diff = sum(norm.(A.*[Xu].*B .- [Yu].*C.*[Zu]))/(N*n);
	end
	return Xu, Yu, Zu
	
end

function batchSolveXY(A, B, opt, nstd_A, nstd_B)
	
	MeanA, SigA = meanCov(A);
	MeanB, SigB = meanCov(B);

	if opt==1
		SigA = SigA - nstd_A * eye(6);
		SigB = SigB - nstd_B * eye(6);
	end

	VA = eigen(SigA[SOneTo(3),SOneTo(3)]).vectors; VA = SMatrix{3,3,Float64,9}(VA);
	VB = eigen(SigB[SOneTo(3),SOneTo(3)]).vectors; VB = SMatrix{3,3,Float64,9}(VB);

	Rx = get_Rx(VA,VB)
	tx = -Rx.*so3_vec.(transpose.((transpose.(Rx).*[SigA[SOneTo(3),SOneTo(3)]].*Rx).\([SigB[SOneTo(3),StaticArrays.SUnitRange(4,6)]].-transpose.(Rx).*[SigA[SOneTo(3),StaticArrays.SUnitRange(4,6)]].*Rx)))
	
	Xs = vcat.(hcat.(Rx,tx),[@SMatrix [0.0 0.0 0.0 1.0]]);
	Ys = [MeanA] .* Xs ./ [MeanB];

	return Xs, Ys, MeanA, MeanB, SigA, SigB
end
@inline so3_vec(X) = @SMatrix [-X[2,3];X[1,3];-X[1,2]]

"""
	get_Rx(VA,VB)
Helper function of batchSolveXY.
"""
function get_Rx(VA, VB)
	Q1 = eye(3);
	Q2 = @SMatrix [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0];
	Q3 = @SMatrix [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0];
	Q4 = @SMatrix [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0];

	Rx_solved = Vector{SMatrix{3,3,Float64,9}}(undef,8);
	Rx_solved[1] = VA*Q1*VB';
	Rx_solved[2]= VA*Q2*VB';
	Rx_solved[3] = VA*Q3*VB';
	Rx_solved[4] = VA*Q4*VB';
	Rx_solved[5] = VA*-Q1*VB';
	Rx_solved[6] = VA*-Q2*VB';
	Rx_solved[7] = VA*-Q3*VB';
	Rx_solved[8] = VA*-Q4*VB';
	return Rx_solved
end

"""
	meanCov(X)

Iteratively solve mean of set {X} , where X ∈ SE(3).
> Yunfeng Wang and Chirikjian, “Nonparametric Second-Order Theory of Error Propagation on Motion Groups.”
"""
function meanCov(X::Vector{SMatrix{4,4,Float64,16}},
				options = (max_iter = 100, tol = 1e-5))
	N = size(X,1);
	
	## Initial Mean Approximation
	Mean = MatrixExp6(sum(MatrixLog6.(X)) ./ N);

	## Iterative Update Mean
	count = 0
	diff = @SMatrix ones(4,4)
	while norm(diff*N) >= options.tol && count <= options.max_iter
		diff = sum(MatrixLog6.([Mean] .\ X)) / N
		Mean = Mean * MatrixExp6(diff)
		count = count + 1
	end

	## Covariance
	diffs = MatrixLog6.([Mean].\X);
	difftwist = vex.(diffs);
	Cov = sum(transpose.(difftwist) .* difftwist) / N;

	return Mean, Cov
end

begin
	roterror(T1::SMatrix{4,4,Float64,16},T2::SMatrix{4,4,Float64,16}) = norm(dmcross(MatrixLog3(T1[SOneTo(3),SOneTo(3)]'*T2[SOneTo(3),SOneTo(3)])));
	tranerror(T1::SMatrix{4,4,Float64,16},T2::SMatrix{4,4,Float64,16}) = norm(T1[SOneTo(3),4].-T2[SOneTo(3),4]);
end
"""
	MbMat_1(A,X,B,Y,C,Z,SigB,SigC)
Helper function of axbyczProb3
"""
function MbMat_1(A,X,B,Y,C,Z,SigB,SigC)
	e1 = @SVector [1,0,0]; e2 = @SVector [0,1,0]; e3 = @SVector [0,0,1];
	RA,PA = T2rt(A); RB,PB = T2rt(B); RC,PC = T2rt(C);
	RX,PX = T2rt(X); RY,PY = T2rt(Y); RZ,PZ = T2rt(Z);
	
	M11 = -RA * RX * mcross(RB*e1);
	M13 = RY * mcross(RC*RZ*e1);
	M15 = RY * RC * RZ * mcross(e1);
	
	M21 = -RA * RX * mcross(RB*e2);
	M23 = RY * mcross(RC*RZ*e2);
	M25 = RY * RC * RZ * mcross(e2);
	
	M31 = -RA * RX * mcross(RB*e3);
	M33 = RY * mcross(RC*RZ*e3);
	M35 = RY * RC * RZ * mcross(e3);
	
	# Translation part
	M41 = -RA * RX * mcross(PB);
	M42 = RA * RX;
	M43 = RY * mcross(RC*Z[SOneTo(3),4] + PC);
	M44 = -RY;
	M46 = -RY * RC * RZ;
	
	M = [M11 zeros(3,3) M13 zeros(3,3)        M15 zeros(3,3);
	    M21 zeros(3,3) M23 zeros(3,3)        M25 zeros(3,3);
	    M31 zeros(3,3) M33 zeros(3,3)        M35 zeros(3,3);
	    M41        M42 M43 M44 zeros(3,3)       M46];
	
	# RHS
	RHS = - A * X * B + Y * C * Z;
	b = vcat(RHS[SOneTo(3),1], RHS[SOneTo(3),2], RHS[SOneTo(3),3], RHS[SOneTo(3),4]);
	## SigBi = Ad^{-1}(Z) * SigCi * Ad^{-T}(Z)
	#First block
	SigB1313 = SigB[SOneTo(3),SOneTo(3)];
	SigB1346 = SigB[SOneTo(3),StaticArrays.SUnitRange(4,6)];
	SigB4646 = SigB[StaticArrays.SUnitRange(4,6),StaticArrays.SUnitRange(4,6)];
	SigB4613 = SigB[StaticArrays.SUnitRange(4,6),SOneTo(3)];
	M55 = -mcross(SigB[SOneTo(3),1]) + SigB1313 * mcross(e1);
	M56 = zeros(3,3);
	M65 = -mcross(SigB[SOneTo(3),2]) + SigB1313 * mcross(e2);
	M66 = zeros(3,3);
	M75 = -mcross(SigB[SOneTo(3),3]) + SigB1313 * mcross(e3);
	M76 = zeros(3,3);
	
	# Second block
	M85 = -mcross(SigB[SOneTo(3),4]) + SigB1346 * mcross(e1);
	M86 = SigB1313 * mcross(e1);
	M95 = -mcross(SigB[SOneTo(3),5]) + SigB1346 * mcross(e2);
	M96 = SigB1313 * mcross(e2);
	M105 = -mcross(SigB[SOneTo(3),6]) + SigB1346 * mcross(e3);
	M106 = SigB1313 * mcross(e3);
	
	# Third block
	M115 = -mcross(SigB[StaticArrays.SUnitRange(4,6),1]) + SigB4613 * mcross(e1);
	M116 = -mcross(SigB[SOneTo(3),1]);
	M125 = -mcross(SigB[StaticArrays.SUnitRange(4,6),2]) + SigB4613 * mcross(e2);
	M126 = -mcross(SigB[SOneTo(3),2]);
	M135 = -mcross(SigB[StaticArrays.SUnitRange(4,6),3]) + SigB4613 * mcross(e3);
	M136 = -mcross(SigB[SOneTo(3),3]);
	
	# Fourth block
	M145 = -mcross(SigB[StaticArrays.SUnitRange(4,6),4]) + SigB4646 * mcross(e1);
	M146 = -mcross(SigB[SOneTo(3),4]) + SigB4613 * mcross(e1);
	M155 = -mcross(SigB[StaticArrays.SUnitRange(4,6),5]) + SigB4646 * mcross(e2);
	M156 = -mcross(SigB[SOneTo(3),5]) + SigB4613 * mcross(e2);
	M165 = -mcross(SigB[StaticArrays.SUnitRange(4,6),6]) + SigB4646 * mcross(e3);
	M166 = -mcross(SigB[SOneTo(3),6]) + SigB4613 * mcross(e3);
	
	M = [M;
	    zeros(3,12)  M55  M56;
	    zeros(3,12)  M65  M66;
	    zeros(3,12)  M75  M76;
	    zeros(3,12)  M85  M86;
	    zeros(3,12)  M95  M96;
	    zeros(3,12) M105 M106;
	    zeros(3,12) M115 M116;
	    zeros(3,12) M125 M126;
	    zeros(3,12) M135 M136;
	    zeros(3,12) M145 M146;
	    zeros(3,12) M155 M156;
	    zeros(3,12) M165 M166];
	
	RHS2 = SE3AdjointInv(Z) * SigC * SE3AdjointInv(Z)' - SigB;
	RHS2 = reshape(RHS2, 3, 12);
	b = [b; RHS2[:]];
	return M,b
end
"""
	MbMat_2(C,Z,Binv,Yinv,A,X,SigB,SigA,B)
Helper function of axbyczProb3
"""
function MbMat_2(C,Z,Binv,Yinv,A,X,SigB,SigA,B)
	e1 = @SVector [1,0,0]; e2 = @SVector [0,1,0]; e3 = @SVector [0,0,1];
	
	Binv = invT(B);
	SigBinv = SE3Adjoint(B) * SigB * SE3Adjoint(B)';
	
	RYinv = Yinv[SOneTo(3),SOneTo(3)];
	RA = A[SOneTo(3),SOneTo(3)];
	RX = X[SOneTo(3),SOneTo(3)];
	RC = C[SOneTo(3),SOneTo(3)];
	RZ = Z[SOneTo(3),SOneTo(3)];
	RBinv = Binv[SOneTo(3),SOneTo(3)];
	M11 = RYinv * RA * RX * mcross(e1);
	M13 = -mcross(RYinv * RA * RX * e1);
	M15 = -RC * RZ * mcross(RBinv*e1);
	
	M21 = RYinv * RA * RX * mcross(e2);
	M23 = -mcross(RYinv * RA * RX * e2);
	M25 = -RC * RZ * mcross(RBinv*e2);
	
	M31 = RYinv * RA * RX * mcross(e3);
	M33 = -mcross(RYinv * RA * RX * e3);
	M35 = -RC * RZ * mcross(RBinv*e3);
	
	# Translation part
	M42 = -RYinv * RA * RX;
	M43 = -mcross(RYinv*RA*X[SOneTo(3),4] + RYinv*A[SOneTo(3),4] + Yinv[SOneTo(3),4]);
	M44 = eye(3);
	M45 = -RC * RZ * mcross(Binv[SOneTo(3),4]);
	M46 = RC * RZ;
	
	M = [       M11 zeros(3,3) M13 zeros(3,3) M15 zeros(3,3);
	            M21 zeros(3,3) M23 zeros(3,3) M25 zeros(3,3);
	            M31 zeros(3,3) M33 zeros(3,3) M35 zeros(3,3);
	     zeros(3,3)        M42 M43        M44 M45       M46];
	
	# RHS
	RHS = - C * Z * Binv + Yinv * A * X;
	
	b = vcat(RHS[SOneTo(3),1], RHS[SOneTo(3),2], RHS[SOneTo(3),3], RHS[SOneTo(3),4]);
	## SigBi^{-1} = Ad^{-1}(X) * SigAi * Ad^{-T}(X)
	# First block
	SigBinv1313 = SigBinv[SOneTo(3),SOneTo(3)];
	SigBinv1346 = SigBinv[SOneTo(3),StaticArrays.SUnitRange(4,6)];
	SigBinv4646 = SigBinv[StaticArrays.SUnitRange(4,6),StaticArrays.SUnitRange(4,6)]
	SigBinv4613 = SigBinv[StaticArrays.SUnitRange(4,6),SOneTo(3)];
	M51 = -mcross(SigBinv[SOneTo(3),1]) + SigBinv1313 * mcross(e1);
	M52 = zeros(3,3);
	M61 = -mcross(SigBinv[SOneTo(3),2]) + SigBinv1313 * mcross(e2);
	M62 = zeros(3,3);
	M71 = -mcross(SigBinv[SOneTo(3),3]) + SigBinv1313 * mcross(e3);
	M72 = zeros(3,3);
	
	# Second block
	M81 = -mcross(SigBinv[SOneTo(3),4]) + SigBinv1346 * mcross(e1);
	M82 = SigBinv1313 * mcross(e1);
	M91 = -mcross(SigBinv[SOneTo(3),5]) + SigBinv1346 * mcross(e2);
	M92 = SigBinv1313 * mcross(e2);
	M101 = -mcross(SigBinv[SOneTo(3),6]) + SigBinv1346 * mcross(e3);
	M102 = SigBinv1313 * mcross(e3);
	
	# Third block
	M111 = -mcross(SigBinv[StaticArrays.SUnitRange(4,6),1]) + SigBinv4613 * mcross(e1);
	M112 = -mcross(SigBinv[SOneTo(3),1]);
	M121 = -mcross(SigBinv[StaticArrays.SUnitRange(4,6),2]) + SigBinv4613 * mcross(e2);
	M122 = -mcross(SigBinv[SOneTo(3),2]);
	M131 = -mcross(SigBinv[StaticArrays.SUnitRange(4,6),3]) + SigBinv4613 * mcross(e3);
	M132 = -mcross(SigBinv[SOneTo(3),3]);
	
	# Fourth block
	M141 = -mcross(SigBinv[StaticArrays.SUnitRange(4,6),4]) + SigBinv4646 * mcross(e1);
	M142 = -mcross(SigBinv[SOneTo(3),4]) + SigBinv4613 * mcross(e1);
	M151 = -mcross(SigBinv[StaticArrays.SUnitRange(4,6),5]) + SigBinv4646 * mcross(e2);
	M152 = -mcross(SigBinv[SOneTo(3),5]) + SigBinv4613 * mcross(e2);
	M161 = -mcross(SigBinv[StaticArrays.SUnitRange(4,6),6]) + SigBinv4646 * mcross(e3);
	M162 = -mcross(SigBinv[SOneTo(3),6]) + SigBinv4613 * mcross(e3);
	
	M = [M;
	     M51  M52 zeros(3,12);
	     M61  M62 zeros(3,12);
	     M71  M72 zeros(3,12);
	     M81  M82 zeros(3,12);
	     M91  M92 zeros(3,12);
	    M101 M102 zeros(3,12);
	    M111 M112 zeros(3,12);
	    M121 M122 zeros(3,12);
	    M131 M132 zeros(3,12);
	    M141 M142 zeros(3,12);
	    M151 M152 zeros(3,12);
	    M161 M162 zeros(3,12)];
	
	RHS2 = SE3AdjointInv(X) * SigA * SE3AdjointInv(X)' - SigBinv;
	RHS2 = reshape(RHS2, 3, 12);
	
	b = [b; RHS2[:]];
	return M,b
end
"""
	dataSpliter(N,n,A,B,C)
Helper function of axbyczProb3
"""
function dataSpliter(N,n,A,B,C)
	k = n>>1;
	A1s = [A[(i-1)*2k+1:(i-1)*2k+k] for i = 1:N];
	B1s = [B[(i-1)*2k+1:(i-1)*2k+k] for i = 1:N];
	C1s = [C[(i-1)*2k+1:(i-1)*2k+k] for i = 1:N];
	A2s = [A[(i-1)*2k+k+1:i*2k] for i = 1:N];
	B2s = [B[(i-1)*2k+k+1:i*2k] for i = 1:N];
	C2s = [C[(i-1)*2k+k+1:i*2k] for i = 1:N];
	return A1s,B1s,C1s,A2s,B2s,C2s
end
###############################################
##              Data Generation              ##
###############################################

function generateSetsOfABC(num, optPDF, Mean, Cov, XActual, YActual, ZActual,unit="m")
	opt = 1;
	A1, B1, C1 = generateABC(num, opt, optPDF, Mean, Cov, XActual, YActual, ZActual,unit);
	opt = 3;
	A2, B2, C2 = generateABC(num, opt, optPDF, Mean, Cov, XActual, YActual, ZActual,unit);
	opt = 2;
	A3, B3, C3 = generateABC(num, opt, optPDF, Mean, Cov, XActual, YActual, ZActual,unit);
	return A1, B1, C1, A2, B2, C2, A3, B3, C3
end

function myGenerateSetsOfABC(num, optPDF, Mean, Cov, Xₘ, Yₘ, Zₘ,unit="m")
	opt = 1;
	A1, B1, C1 = myGenerateABC(num, opt, optPDF, Mean, Cov, Xₘ, Yₘ, Zₘ,unit);
	opt = 3;
	A2, B2, C2 = myGenerateABC(num, opt, optPDF, Mean, Cov, Xₘ, Yₘ, Zₘ,unit);
	return A1, B1, C1, A2, B2, C2
end

"""
	generateABC(num, optFix, optPDF, M, Sig, X, Y, Z)
Data generation for AXB = YCZ problem
"""
function generateABC(m, optFix, optPDF, M, Sig, X, Y, Z,unit="m")
    a = @SMatrix randn(6,1); a = a./norm(a); A_initial = MatrixExp6(vee(a));
    b = @SMatrix randn(6,1); b = b./norm(b); B_initial = MatrixExp6(vee(b));
    c = @SMatrix randn(6,1); c = c./norm(c); C_initial = MatrixExp6(vee(c));
	
	if optFix == 1 # Fix A, randomize B and C
		if optPDF == 1
			se3s = [mvg(M, Sig, 1) for i = 1:m]
			B = MatrixExp6.(vee.(se3s)) .* [B_initial];
		elseif optPDF == 2
			se3s = [mvg(M, Sig, 1) for i = 1:m]
			B = [B_initial] .* MatrixExp6.(vee.(se3s));
		elseif optPDF == 3
			gmean = [(@SMatrix zeros(6,1)) for i = 1:m]
			B = sensorNoise([B_initial], gmean, [Sig(1)], 1)
		end
		C =  [Y] .\ ([A_initial * X] .* B ./ [Z])
		A = [A_initial for i = 1:m]
	elseif optFix == 2  # Fix B, randomize A and C
		if optPDF == 1
			se3s = [mvg(M, Sig, 1) for i = 1:m]
			A = MatrixExp6.(vee.(se3s)) .* [A_initial];
		elseif optPDF == 2
			se3s = [mvg(M, Sig, 1) for i = 1:m]
			A = [A_initial] .* MatrixExp6.(vee.(se3s));
		elseif optPDF == 3
			gmean = [(@SMatrix zeros(6,1)) for i = 1:m]
			A = sensorNoise([A_initial], gmean, [Sig(1)], 1)
		end
		C = [Y] .\ (A .* [X * B_initial] ./ [Z])
		B = B_initial;
	elseif optFix == 3 # Fix C, randomize A and B
		if optPDF == 1
			se3s = [mvg(M, Sig, 1) for i = 1:m]
			B = MatrixExp6.(vee.(se3s)) .* [B_initial];
		elseif optPDF == 2
			se3s = [mvg(M, Sig, 1) for i = 1:m]
			B = [B_initial] .* MatrixExp6.(vee.(se3s));
		elseif optPDF == 3
			gmean = [(@SMatrix zeros(6,1)) for i = 1:m]
			B = sensorNoise([B_initial], gmean, [Sig(1)], 1)
		end
		A = ([Y*C_initial*Z] ./ B) ./ [X];
		C = [C_initial for i = 1:m]
	end

	return A,B,C
			
end


function myGenerateABC(m, optFix, optPDF, M, Sig, Xₘ, Yₘ, Zₘ,unit="m")
	qa = @SVector rand(6); A_initial = UR10.fkine(qa,unit); A_initialₘ = UR10.fkineGAS(qa,unit);
    b = @SMatrix randn(6,1); b = b./norm(b); B_initial = MatrixExp6(vee(b));
	B_initialₘ = tfm2motorS(B_initial);
	qc = @SVector rand(6); C_initial = UR10.fkine(qc,unit); C_initialₘ = UR10.fkineGAS(qa,unit);

	if optFix == 1 # Fix A, randomize B and C
		if optPDF == 1
			se3s = [mvg(M, Sig, 1) for i = 1:m]
			B = MatrixExp6.(vee.(se3s)) .* [B_initial];
		elseif optPDF == 2
			se3s = [mvg(M, Sig, 1) for i = 1:m]
			B = [B_initial] .* MatrixExp6.(vee.(se3s));
		elseif optPDF == 3
			gmean = [(@SMatrix zeros(6,1)) for i = 1:m]
			B = sensorNoise([B_initial], gmean, [Sig(1)], 1)
		end
		if unit == "m"
			B_initial = tometer(B_initial)
		end
		A = [A_initialₘ for i = 1:m]
		B = match_signS.([B_initialₘ], tfm2motorS.(B))
		C = ecmul.([ereversion(Zₘ)],ecmul.(B,[Xₘ],A),[ereversion(Yₘ)])
	elseif optFix == 3 # Fix C, randomize A and B
		if optPDF == 1
			se3s = [mvg(M, Sig, 1) for i = 1:m]
			B = MatrixExp6.(vee.(se3s)) .* [B_initial];
		elseif optPDF == 2
			se3s = [mvg(M, Sig, 1) for i = 1:m]
			B = [B_initial] .* MatrixExp6.(vee.(se3s));
		elseif optPDF == 3
			gmean = [(@SMatrix zeros(6,1)) for i = 1:m]
			B = sensorNoise([B_initial], gmean, [Sig(1)], 1)
		end
		if unit == "m"
			B_initial = tometer(B_initial)
		end
		C = [C_initialₘ for i = 1:m]
		B = match_signS.([B_initialₘ], tfm2motorS.(B))
		A = ecmul.(ereversion.(ecmul.(B,[Xₘ])),ecmul.([Zₘ],C,[Yₘ]))
	end

	return A,B,C
end


"""
	mvg(mu,Sigma, N)
Multivariate Gaussian random number generator.
"""
function mvg(mu, Sigma, N::Int)
	@assert size(mu,1)==size(Sigma,1) "Length(mu) must equal size(Sigma,1)"
	@assert size(Sigma,1)==size(Sigma,2) "Sigma must be square"
	@assert norm(Sigma-Sigma')<1e-15 "Sigma must be symmetric."
	_,R = cholesky(Sigma)
	n = length(mu);
	y = R'* (@SMatrix randn(n, N)) + hcat([mu for i = 1:N]...)
	return y
end

"""
	get_data_Maₘ(n, n_angle, n_dis)
Generate n {A,B,C} in motor that satisfies AXB=YCZ based on Ma's requirement.
First half of the data shares a static `A` and the other half a static `C`. 

*Noted that `Xₘ, Yₘ, Zₘ` which are the motor form of `X,Y,Z` respctively, and `unit` should be provided in workspace.*
"""
function get_data_Maₘ(n, Xₘ::T,Yₘ::T,Zₘ::T, unit) where {T<:SVector{8,Float64}}
	# Fix A
	Q1 = rand(Uniform(-π,π),1,6); Q1a = [Q1 for i =1:n>>1];
	A1ₘ = UR10.fkineGAS(Q1,unit); A1ₘₛ = [A1ₘ for i = 1:n>>1];

	Q1c = [rand(Uniform(-π,π),1,6) for i = 1:n>>1]; 	# Random C
	C1ₘₛ = UR10.fkineGAS.(Q1c,unit);
	
	B1ₘₛ = ecmul.([Zₘ],ecmul.(C1ₘₛ,ecmul.([Yₘ],ecmul.(ereversion.(A1ₘₛ),[ereversion(Xₘ)]))))
	
	# Fix C
	Q2 = rand(Uniform(-π,π),1,6); Q2c = [Q2 for i =1:n>>1];
	C2 = UR10.fkineGAS(Q2,unit); C2ₘₛ = [C2 for i = 1:n>>1];
	
	Q2a = [rand(Uniform(-π,π),1,6) for i = 1:n>>1]; 		# Random A
	A2ₘₛ = UR10.fkineGAS.(Q2a,unit);
	
	B2ₘₛ = ecmul.([Zₘ],ecmul.(C2ₘₛ,ecmul.([Yₘ],ecmul.(ereversion.(A2ₘₛ),[ereversion(Xₘ)]))))
	
	Q₁ = vcat(Q1a,Q2a); Q₂ = vcat(Q1c, Q2c);
	Aₘₛ = vcat(A1ₘₛ, A2ₘₛ); Bₘₛ = vcat(B1ₘₛ, B2ₘₛ); Cₘₛ = vcat(C1ₘₛ, C2ₘₛ);

	return Aₘₛ, Bₘₛ, Cₘₛ, Q₁, Q₂
end

function get_mvg_dataSₘ(n,angle_noise, dis_noise,Xₘ,Yₘ,Zₘ,fix,unit="m",dist="normal")
	m = n>>1;
	gmean = @SMatrix zeros(6,1);
	k = 0.02;
	Cov = eye(6);
	A1,B1,C1,A2,B2,C2 = myGenerateSetsOfABC(m,1,gmean,k*Cov,Xₘ,Yₘ,Zₘ,unit)
	## Add Noise for Set 1
	a1 = add_noiseSₘ(A1[1],angle_noise,dis_noise,fix,dist)
	A1 = [a1 for i = 1:m]
	B1 = add_noiseSₘ.(B1,angle_noise,dis_noise,fix,dist);
	C1 = add_noiseSₘ.(C1,angle_noise,dis_noise,fix,dist);
	## Add Noise for Set 2
	A2 = add_noiseSₘ.(A2,angle_noise,dis_noise,fix,dist);
	B2 = add_noiseSₘ.(B2,angle_noise,dis_noise,fix,dist);
	c2 = add_noiseSₘ(C2[1],angle_noise,dis_noise,fix,dist)
	C2 = [c2 for i = 1:m]
	## Collect 
	Aₘₛ = vcat(A1, A2); Bₘₛ = vcat(B1, B2); Cₘₛ = vcat(C1, C2);

	Q₁ = zeros(6,1); Q₂ = zeros(6,1);
	return Aₘₛ, Bₘₛ, Cₘₛ, Q₁, Q₂
end

function get_N_mvg_dataSₘ(N,n,angle_noise, dis_noise, Xₘ, Yₘ, Zₘ, fix,unit,dist="normal")
	func(n) = get_mvg_dataSₘ(n,angle_noise, dis_noise, Xₘ, Yₘ, Zₘ, fix,unit,dist)
	ABCQ = func.([n for i = 1:N])
	Aₘₛ = vcat([ABCQ[i][1] for i = 1:N]...)
	Bₘₛ = vcat([ABCQ[i][2] for i = 1:N]...)
	Cₘₛ = vcat([ABCQ[i][3] for i = 1:N]...)
	Q₁ = vcat([ABCQ[i][4] for i = 1:N]...)
	Q₂ = vcat([ABCQ[i][5] for i = 1:N]...)
	return Aₘₛ,Bₘₛ,Cₘₛ,Q₁,Q₂
end

function add_mvg_NoiseS(A::SMatrix{4,4,Float64,16},angle_noise, dis_noise)
	M = @SMatrix zeros(6,1); Sig = SMatrix{6,6,Float64,36}(diagm(vcat([1,1,1]*angle_noise,[1,1,1]*dis_noise)))
	se3s = mvg(M, Sig, 1)
	An = A * MatrixExp6(vee(se3s))
end

function selectMaData(NN,nn,ns,A,B,C)
	AS = vcat([A[(i-1)*(nn>>1)+1:(i-1)*(nn>>1)+ns>>1] for i = 1:2*NN]...)
	BS = vcat([B[(i-1)*(nn>>1)+1:(i-1)*(nn>>1)+ns>>1] for i = 1:2*NN]...)
	CS = vcat([C[(i-1)*(nn>>1)+1:(i-1)*(nn>>1)+ns>>1] for i = 1:2*NN]...)
	return AS,BS,CS
end

function selectMaDataR(NN,nn,ns,A,B,C)
	AS = Vector{SMatrix{4,4,Float64,16}}(undef,NN*ns);
	BS = Vector{SMatrix{4,4,Float64,16}}(undef,NN*ns);
	CS = Vector{SMatrix{4,4,Float64,16}}(undef,NN*ns);
	for i = 1:NN
		s1 = sample(1:nn>>1,ns>>1,replace=false)
		s2 = sample(1:nn>>1,ns>>1,replace=false)
		AS[(i-1)*ns+1:(i-1)*ns+ns>>1] = A[s1.+(i-1)*nn]; 
		AS[(i-1)*ns+ns>>1+1:i*ns] = A[s2.+((i-1)*nn+nn>>1)]; 
		BS[(i-1)*ns+1:(i-1)*ns+ns>>1] = B[s1.+(i-1)*nn]; 
		BS[(i-1)*ns+ns>>1+1:i*ns] = B[s2.+((i-1)*nn+nn>>1)]; 
		CS[(i-1)*ns+1:(i-1)*ns+ns>>1] = C[s1.+(i-1)*nn]; 
		CS[(i-1)*ns+ns>>1+1:i*ns] = C[s2.+((i-1)*nn+nn>>1)]; 
	end
	return AS,BS,CS
end


###############################################
##               Hybrid Method               ##
###############################################
"""
	AXBYCZ_Wang2(A,B,C,X,Y,Z)
Function that solves AXB=YCZ with method given in 
"towards simultaneous coordinates calibrations for cooperative multiple robots"

"""
function AXBYCZ_Wang2(A::T1,B::T1,C::T1,X,Y,Z) where {T1<:Vector{SMatrix{4,4,Float64,16}}}
	n = size(A,1);
	RA,tA = T2rt.(A);
	RB,tB = T2rt.(B);
	RC,tC = T2rt.(C);
	RA = [A[i][SOneTo(3),SOneTo(3)] for i = 1:n]
	RB = [B[i][SOneTo(3),SOneTo(3)] for i = 1:n]
	RC = [C[i][SOneTo(3),SOneTo(3)] for i = 1:n]
	tA = [A[i][SOneTo(3),4] for i = 1:n]
	tB = [B[i][SOneTo(3),4] for i = 1:n]
	tC = [C[i][SOneTo(3),4] for i = 1:n]

	# e = π/5 .* @SVector ones(3,1);
	# RX_init = MatrixExp6(mcross(e))*X[SOneTo(3),SOneTo(3)];
	# RY_init = MatrixExp6(mcross(e))*X[SOneTo(3),SOneTo(3)];
	# RX_init = MatrixExp6(mcross(e))*X[SOneTo(3),SOneTo(3)];
	RX_init,_ = T2rt(X);RY_init,_ = T2rt(Y);RZ_init,_ = T2rt(Z);
	delR = 1e4 * @SVector ones(9);
	local n_step = 0;
	while (norm(delR) > 1e-4 && n_step < 500)
		q = zeros(n*9,1);
		F = zeros(n*9,9);
		@inbounds for i = 1:n
			tmp1 = RX_init * RB[i];
			tmp2 = RY_init * RC[i] * RZ_init;
			qq = -RA[i]*tmp1 + tmp2;
			q[(i-1)*9+1:i*9] = qq[:];

			F11 = -RA[i] * mcross(tmp1[:,1]);
			F21 = -RA[i] * mcross(tmp1[:,2]);
			F31 = -RA[i] * mcross(tmp1[:,3]);
			F12 = mcross(tmp2[:,1]);
			F22 = mcross(tmp2[:,2]);
			F32 = mcross(tmp2[:,3]);
			F13 = RY_init * RC[i] * mcross(RZ_init[:,1])
			F23 = RY_init * RC[i] * mcross(RZ_init[:,2])
			F33 = RY_init * RC[i] * mcross(RZ_init[:,3])
			F[(i-1)*9+1:i*9,:] = [F11 F12 F13; F21 F22 F23; F31 F32 F33];
		end
		delR = SVector{9,Float64}(F\q);
		RX_init = MatrixExp3(mcross(delR[SOneTo(3)]))*RX_init;
		RY_init = MatrixExp3(mcross(delR[StaticArrays.SUnitRange(4,6)]))*RY_init;
		RZ_init = MatrixExp3(mcross(delR[StaticArrays.SUnitRange(7,9)]))*RZ_init;
		n_step = n_step + 1;
	end

	J = zeros(3*n,9); p = zeros(3*n,1);
	for i = 1:n
		J[(i-1)*3+1:i*3,:] = [RA[i] -eye(3) -RY_init*RC[i]];
		p[(i-1)*3+1:i*3] = -tA[i] - RA[i]*RX_init*tB[i] + RY_init*tC[i];
	end
	translation = SVector{9,Float64}(J\p);
	tX = translation[SOneTo(3)];
	tY = translation[StaticArrays.SUnitRange(4,6)];
	tZ = translation[StaticArrays.SUnitRange(7,9)];
	return rt2T(RX_init,tX), rt2T(RY_init,tY), rt2T(RZ_init,tZ)
end