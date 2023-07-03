
function add_noiseSₘ(M::SVector,n_angle,n_dis,fix=false,dist="uniform")
	T = rand_translatorS(n_dis,fix,dist)
	R = rand_rotorS(n_angle,fix,dist)
	return ecmul(R,ecmul(M,T))
end

function get_error(A::SVector{8,Float64},B::SVector{8,Float64})
	get_error(motorS2tfm(A),motorS2tfm(B))
end

function get_error(M,Mₜ)
	rₑᵣᵣ = get_rot_error(M[SOneTo(3),SOneTo(3)],Mₜ[SOneTo(3),SOneTo(3)])
	tₑᵣᵣ = norm(M[SOneTo(3),4]-Mₜ[SOneTo(3),4])
	return rₑᵣᵣ,tₑᵣᵣ
end

function get_error(M::SVector,Mₜ::Union{Matrix,SMatrix})
	Mₘ = motorS2tfm(M)
	return get_error(Mₘ,Mₜ)
end

function get_error(M::Union{Matrix,SMatrix},Mₜ::SVector)
	Mmₜ = motorS2tfm(Mₜ)
	return get_error(M,Mmₜ)
end

function get_error_result(Zₛ,Yₛ,Zₜ,Yₜ)
	errᵣ₁,errₜ₁=get_error(Zₛ,Zₜ)
	errᵣ₂,errₜ₂=get_error(Yₛ,Yₜ)
	return [errᵣ₁,errₜ₁,errᵣ₂,errₜ₂]
end

function get_rot_error(R::Union{Matrix,SMatrix},Rₜ::Union{Matrix,SMatrix})
	Rₙ = ProjectToSO3(R)
	err = AngleAxis(Rₙ*Rₜ').theta
	# err = opnorm((R'*Rₜ)-diagm([1.0,1.0,1.0]))
end

function tometer(X)
	return @SMatrix [X[1,1] X[1,2] X[1,3] X[1,4]/1000.0;
					X[2,1] X[2,2] X[2,3] X[2,4]/1000.0;
					X[3,1] X[3,2] X[3,3] X[3,4]/1000.0;
					0.0 0.0 0.0 1.0]
end

function tommeter(X)
	return @SMatrix [X[1,1] X[1,2] X[1,3] X[1,4]*1000.0;
					X[2,1] X[2,2] X[2,3] X[2,4]*1000.0;
					X[3,1] X[3,2] X[3,3] X[3,4]*1000.0;
					0.0 0.0 0.0 1.0]
end

function get_report(Init,Ans,Result)
	dis = zeros(3); Angle = zeros(3); rdis =zeros(3); rangle =zeros(3);
	odis = zeros(3); oangle =zeros(3)
	for i = 1:3
		Origin_err = get_error(Ans[i],Init[i]);
		Result_err = get_error(Ans[i],Result[i]);
		dis[i] = norm(Result_err[2]); odis[i] = norm(Origin_err[2]);
		rdis[i] = round((norm(Origin_err[2])-norm(Result_err[2]))/norm(Origin_err[2])*100; sigdigits=4);
		Angle[i] = Result_err[1]; oangle[i] = Origin_err[1];
		rangle[i] = round(((Origin_err[1]-Result_err[1])/Origin_err[1])*100; sigdigits=4);
	end
	dis = round.(dis; sigdigits = 4);
	odis = round.(odis; sigdigits = 4);
	Angle = round.(Angle; sigdigits = 4);
	oangle = round.(oangle; sigdigits = 4);
	return [dis,odis,rdis, Angle,oangle,rangle]
end;