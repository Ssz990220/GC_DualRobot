using Printf
using Distributed
using SharedArrays
using JLD


function benchmark(noise_angle, noise_dis, n_sample, n_time::Int,get_XYZ,get_error_func, conf)
    @assert size(noise_angle,1) == size(noise_dis,1);
    err = zeros(size(noise_angle,1),size(n_sample,1),6)
	@inbounds for i = 1:size(n_sample,1)
		@inbounds for j = 1:size(noise_angle,1)
		    n_angle = noise_angle[j]; n_dis = noise_dis[j];
            local_err = [zeros(1,n_time) for i = 1:6]
		    @inbounds for k = 1:n_time
                X_,Y_,Z_ = get_XYZ(n_sample[i],n_angle,n_dis)
                errs = get_error_func(X_,Y_,Z_)
                for m = 1:6
                    local_err[m][k] = errs[m];
                end
			end
			# println("noise level $i, sample level $j, is done...")
            @inbounds for m = 1:6
                err[j,i,m] = Statistics.mean(abs.(local_err[m]))
            end
		end
	end
	return err
end

"""
    benchmark_p(noise_angle, noise_dis, n_sample, n_time::Int,get_data,get_XYZs,get_error_func, conf)
**Deprecated Function**

Parallel Experiment With Multi-Processing.
"""
function benchmark_p(noise_angle, noise_dis, n_sample, n_time::Int,get_data,get_XYZs,get_error_func, err_dim, conf; disp = "sample")
    # @show get_XYZs
    nfunc = size(get_XYZs,1)
    @assert size(noise_angle,1) == size(noise_dis,1);
    err = SharedArray{Float64,5}(size(noise_angle,1),size(n_sample,1),err_dim,n_time,nfunc)
    if disp=="sample"
        nsum = sum(n_sample);
    elseif disp == "noise"
        nsum = size(noise_angle,1)
    end
    finised = 0;
    totalTime = 0.0;
	for i ∈ eachindex(n_sample)
        if disp=="sample"
            start = time();
        end
		for j ∈ eachindex(noise_angle)
            if disp=="noise"
                start = time();
            end
		    n_angle = noise_angle[j]; n_dis = noise_dis[j];
		    @sync @distributed for k = 1:n_time
                ValidationData = get_data(100,n_angle,n_dis)
                for s ∈ eachindex(get_XYZs)
                    data = get_data(n_sample[i],n_angle,n_dis)
                    X,Y,Z = get_XYZs[s](data,conf[s])
                    errs = get_error_func(X,Y,Z,ValidationData)
                    for m = 1:err_dim
                        err[j,i,m,k,s] = errs[m];
                    end
                end
			end
            if disp=="noise"
                finished = finised + 1
                etime = time();
                batchTime = etime - start;
                totalTime = totalTime + batchTime
                ETA = (nsum - finished) * batchTime;
                @printf "Batch %d finished, \tbatch time is %.4f s, \ttotal elipsed time %.4fs, \tETA %.4f s\n" j batchTime totalTime ETA;
            end
		end
        if disp=="sample"
            finished = finised + n_sample[i]
            etime = time();
            batchTime = etime - start;
            totalTime = totalTime + batchTime
            ETA = (nsum - finished) * batchTime / n_sample[i];
            @printf "Batch %d finished, \tbatch time is %.4f s, \ttotal elipsed time %.4fs, \tETA %.4f s\n" n_sample[i] batchTime totalTime ETA;
        end
    end
	return err
end


begin
    function error_func(X,Y,Z,ValidationData,Ans)
        if size(X,1) == 8
            XMatrix = motorS2tfm(X);YMatrix = motorS2tfm(Y);ZMatrix = motorS2tfm(Z);
        else
            XMatrix = X; YMatrix = Y; ZMatrix = Z;
        end
        A = motorS2tfm.(ValidationData.P[1]); B = motorS2tfm.(ValidationData.P[2]); C = motorS2tfm.(ValidationData.P[3]);
        errs = @MVector zeros(30); errts = @MVector zeros(30);
        @inbounds for i = 1:30
            try
                errs[i],errts[i] = get_error(A[i]*XMatrix*B[i],YMatrix*C[i]*ZMatrix)
            catch e
                @show XMatrix, YMatrix, ZMatrix
            end
        end
        ErrorToTrue = get_error_XYZ(XMatrix, YMatrix, ZMatrix, Ans)
        return vcat(ErrorToTrue,[Statistics.mean(errs), Statistics.mean(errts)])
    end
    
    function liao(data,conf)
        A = data.P[1]; B = data.P[2]; C = data.P[3];
        Xᵣ,Yᵣ,Zᵣ = get_XYZ_rotor_Close(A,B,C)
        Xᵣ = r2t(rotorG32rm(Xᵣ));Yᵣ = r2t(rotorG32rm(Yᵣ));Zᵣ = r2t(rotorG32rm(Zᵣ));
        Aₘ = data.P[1]; Bₘ = data.P[2]; Cₘ = data.P[3];
        A = motorS2tfm.(Aₘ); B = motorS2tfm.(Bₘ); C = motorS2tfm.(Cₘ);
        A = tommeter.(A); B = tommeter.(B); C = tommeter.(C);
        n = size(A,1)
        RA = [A[i][SOneTo(3),SOneTo(3)] for i = 1:n];
        RB = [B[i][SOneTo(3),SOneTo(3)] for i = 1:n];
        RC = [C[i][SOneTo(3),SOneTo(3)] for i = 1:n];
        TA = [A[i][SOneTo(3),4] for i = 1:n];
        TB = [B[i][SOneTo(3),4] for i = 1:n];
        TC = [C[i][SOneTo(3),4] for i = 1:n];
        RX_sln,RY_sln,RZ_sln = liao_iter(RA,RB,RC,Xᵣ[SOneTo(3),SOneTo(3)],Yᵣ[SOneTo(3),SOneTo(3)],Zᵣ[SOneTo(3),SOneTo(3)],conf);
        Xᵣ,Yᵣ,Zᵣ = T_SVD(RX_sln,RY_sln,RZ_sln, RA,RB,RC,TA,TB,TC)
        Xᵣ = tometer(Xᵣ); Yᵣ = tometer(Yᵣ); Zᵣ = tometer(Zᵣ);
        return Xᵣ,Yᵣ,Zᵣ
    end
    
    # function wang_full(data,conf)
    #     A = data.P[1]; B = data.P[2]; C = data.P[3];
    #     A = motor2tfm.(A); B = motor2tfm.(B); C = motor2tfm.(C);
    #     A = tommeter.(A); B = tommeter.(B); C = tommeter.(C);
    #     Xᵣ,Yᵣ,Zᵣ = AXB_YCZ_Close(A,B,C,conf)
    #     V = toSM.([Xᵣ,Yᵣ,Zᵣ])
    #     Result,err = iter_solve_Katyusha(V,[A,B,C],conf)
    #     Xᵣ,Yᵣ,Zᵣ = post_process_Wang(Result)
    #     Xᵣ = tometer(Xᵣ); Yᵣ = tometer(Yᵣ); Zᵣ = tometer(Zᵣ);
    #     return Xᵣ,Yᵣ,Zᵣ
    # end
    
    # function close_wang(data,conf)
    #     A = data.P[1]; B = data.P[2]; C = data.P[3];
    #     A = motor2tfm.(A); B = motor2tfm.(B); C = motor2tfm.(C);
    #     Xᵣ,Yᵣ,Zᵣ = AXB_YCZ_Close(A,B,C,conf)
    #     return Xᵣ,Yᵣ,Zᵣ
    # end
    
    function GA_Full(data,conf)
        A = data.P[1]; B = data.P[2]; C = data.P[3];
        Q₁ = data.Q[1]; Q₂ = data.Q[2]
        Xᵣ,Yᵣ,Zᵣ = GAClose(A,B,C)
        # A_ = motor2tfm.(A); B_ = motor2tfm.(B); C_ = motor2tfm.(C);
        # X2ᵣ,Y2ᵣ,Z2ᵣ = AXB_YCZ_Close(A_,B_,C_,conf)
        # Xᵣ = match_sign(Xᵣ,tfm2motor(X2ᵣ));Yᵣ = match_sign(Yᵣ,tfm2motor(Y2ᵣ));Zᵣ = match_sign(Zᵣ,tfm2motor(Z2ᵣ))
        result,err = iter_solve_GC_katyushaᵥ([Xᵣ,Yᵣ,Zᵣ],[A,B,C],conf)
        return result[1], result[2], result[3]
    end
    
    # function GA_Close(data,conf)
    #     A = data.P[1]; B = data.P[2]; C = data.P[3];
    #     Q₁ = data.Q[1]; Q₂ = data.Q[2]
    #     Xᵣ,Yᵣ,Zᵣ = GAClose(A,B,C,Xₘ,Yₘ,Zₘ,Q₁,Q₂,unit)
    #     return Xᵣ,Yᵣ,Zᵣ
    # end
    
    """
        G3_Close(data,conf)
    Wrapper for proposed G3 initial value finder in benchmark
    """
    function G3_Close(data,conf)
        A = data.P[1]; B = data.P[2]; C = data.P[3];
        Xᵣ,Yᵣ,Zᵣ = get_XYZ_rotor_Close(A,B,C)
        Rxᵢₙᵢₜ = rotorG32rm(Xᵣ); Ryᵢₙᵢₜ = rotorG32rm(Yᵣ); Rzᵢₙᵢₜ = rotorG32rm(Zᵣ); 
        A = motorS2tfm.(A); B = motorS2tfm.(B); C = motorS2tfm.(C);  
        n = size(A,1)
        RA = [A[i][SOneTo(3),SOneTo(3)] for i = 1:n];
        RB = [B[i][SOneTo(3),SOneTo(3)] for i = 1:n];
        RC = [C[i][SOneTo(3),SOneTo(3)] for i = 1:n];
        TA = [A[i][SOneTo(3),4] for i = 1:n];
        TB = [B[i][SOneTo(3),4] for i = 1:n];
        TC = [C[i][SOneTo(3),4] for i = 1:n];
        return T_SVD(Rxᵢₙᵢₜ,Ryᵢₙᵢₜ,Rzᵢₙᵢₜ,RA,RB,RC,TA,TB,TC)
    end
    
    """
        G3(data,conf)
    Wrapper for proposed G3 AXB=YCZ solver in benchmark
    """
    function G3(data,conf)
        A = data.P[1]; B = data.P[2]; C = data.P[3];
        # Q₁ = data.Q[1]; Q₂ = data.Q[2]
        # Xᵣ,Yᵣ,Zᵣ = GAClose(A,B,C,Xₘ,Yₘ,Zₘ,Q₁,Q₂,unit)
        Xᵣ,Yᵣ,Zᵣ = get_XYZ_rotor_Close(A,B,C)       ## A closed-form solver to find initial value
        result,err = iter_solve_Katyusha_G3([Xᵣ,Yᵣ,Zᵣ],data.P,conf)     ## Iteratively update
    
        ## Solve Translation with SVD
        Rxᵢₙᵢₜ = rotorG32rm(result[1]); Ryᵢₙᵢₜ = rotorG32rm(result[2]); Rzᵢₙᵢₜ = rotorG32rm(result[3]);
        A = motorS2tfm.(A); B = motorS2tfm.(B); C = motorS2tfm.(C);  
        n = size(A,1)
        RA = [@view A[i][SOneTo(3),SOneTo(3)] for i = 1:n];
        RB = [@view B[i][SOneTo(3),SOneTo(3)] for i = 1:n];
        RC = [@view C[i][SOneTo(3),SOneTo(3)] for i = 1:n];
        TA = [@view A[i][SOneTo(3),4] for i = 1:n];
        TB = [@view B[i][SOneTo(3),4] for i = 1:n];
        TC = [@view C[i][SOneTo(3),4] for i = 1:n];
        return T_SVD(Rxᵢₙᵢₜ,Ryᵢₙᵢₜ,Rzᵢₙᵢₜ,RA,RB,RC,TA,TB,TC)
        # return T_EKF(Rxᵢₙᵢₜ,Ryᵢₙᵢₜ,Rzᵢₙᵢₜ,RA,RB,RC,TA,TB,TC)
    end
    """
        wang_cic(data,conf)
    Wrapper for Wang's AXB=YCZ solver in benchmark.
    Notice that the translation part of the final solution comes from an addtional SVD step.
    """
    function wang_cic(data,conf)
        A = data.P[1]; B = data.P[2]; C = data.P[3];
        A = motorS2tfm.(A); B = motorS2tfm.(B); C = motorS2tfm.(C);
        A = tommeter.(A); B = tommeter.(B); C = tommeter.(C);
        Xᵣ,Yᵣ,Zᵣ = AXB_YCZ_Wang(A,B,C,conf)
        return tometer(Xᵣ),tometer(Yᵣ),tometer(Zᵣ)
    end
    
    function fu_bayro(data,conf)
        n = size(data.P[1],1);
        k = n; k2 = n;
        A = data.P[1][1:k2]; B = data.P[2][1:k2]; C = data.P[3][1:k2];
        Ax = data.axxb[1][1:k]; Bx = data.axxb[2][1:k];
        Xᵢₙᵢₜ = AXXBs(Ax,Bx);
        Yᵢ,Zᵢ = AXB_YCZ_FUS(A,B,C,Xᵢₙᵢₜ)
        return Xᵢₙᵢₜ,Yᵢ,Zᵢ
    end
    
    function fu_bayro13(data,conf)
        n = size(data.P[1],1);
        k = n÷3; k2 = n*2÷3;
        A = data.P[1][1:k2]; B = data.P[2][1:k2]; C = data.P[3][1:k2];
        Ax = data.axxb[1][1:k]; Bx = data.axxb[2][1:k];
        Xᵢₙᵢₜ = AXXBs(Ax,Bx);
        Yᵢ,Zᵢ = AXB_YCZ_FUS(A,B,C,Xᵢₙᵢₜ)
        return Xᵢₙᵢₜ,Yᵢ,Zᵢ
    end

    function Ma(data,conf)
        A = motorS2tfm.(data.P[1]); B = motorS2tfm.(data.P[2]); C = motorS2tfm.(data.P[3]);
        N = conf.N;
        n = size(data.P[1],1)÷N;
        A1 = A[1:n>>1];B1 = B[1:n>>1];C1 = C[1:n>>1];
        A2 = A[n>>1+1:n];B2 = B[n>>1+1:n];C2 = C[n>>1+1:n];
        Xₙ,Yₙ,Zₙ = axbyczProb1(A1,B1,C1,A2,B2,C2,0,0,0)
        AMa_mean = [meanCov(A[(i-1)*(n>>1)+1:i*(n>>1)])[1] for i = 1:2*N]
        BMa_mean = [meanCov(B[(i-1)*(n>>1)+1:i*(n>>1)])[1] for i = 1:2*N]
        CMa_mean = [meanCov(C[(i-1)*(n>>1)+1:i*(n>>1)])[1] for i = 1:2*N]
        Xₙ,Yₙ,Zₙ = axbyczProb3(N,n,A,B,C,Xₙ,Yₙ,Zₙ)
        return Xₙ,Yₙ,Zₙ
    end

    
    function MaOne(data,conf)
        A = motorS2tfm.(data.P[1]); B = motorS2tfm.(data.P[2]); C = motorS2tfm.(data.P[3]);
        n = size(A,1);
        A1 = A[1:n>>1];B1 = B[1:n>>1];C1 = C[1:n>>1];
        A2 = A[n>>1+1:n];B2 = B[n>>1+1:n];C2 = C[n>>1+1:n];
        Xₙ,Yₙ,Zₙ = axbyczProb1(A1,B1,C1,A2,B2,C2,0,0,0)
        return Xₙ,Yₙ,Zₙ
    end

end