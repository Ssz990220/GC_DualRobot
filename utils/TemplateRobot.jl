export TEMPLATEROBOT
module TEMPLATEROBOT
	using StaticArrays
	using Rotations
	include("Alge.jl")
	include("GA.jl")
	include("Robotics.jl")
	include("GA_robotics.jl")
     			  #theta a(x)  	 d(z)	 alpha(x)
	DH = @SMatrix [];			# Fill in DH for robot, length in mm, angle in rad.
	function fkine(qs,unit="mm")
		if unit=="mm"
			dh = DH;
		else
			dh = [DH[:,1] DH[:,2]/1000 DH[:,3]/1000 DH[:,4]]
		end
		FkineDH(qs,dh);
	end
	# function fkineGA(qs,unit="mm")
	# 	if unit=="mm"
	# 		dh = DH;
	# 	else
	# 		dh = [DH[:,1] DH[:,2]/1000 DH[:,3]/1000 DH[:,4]]
	# 	end
	# 	FkineGA(qs,dh);
	# end

	function fkineGAS(qs, unit="mm")
		if unit=="mm"
			dh = DH;
		else
			dh = [DH[:,1] DH[:,2]/1000 DH[:,3]/1000 DH[:,4]]
		end
		FkineGAS(qs,dh);
	end

	function fkineG3S(qs, unit="mm")
		if unit=="mm"
			dh = DH;
		else
			dh = [DH[:,1] DH[:,2]/1000 DH[:,3]/1000 DH[:,4]]
		end
		FkineG3S(qs,dh);
	end

	function jacobe(qs,unit="mm")
		if unit=="mm"
			dh = DH;
		else
			dh = [DH[:,1] DH[:,2]/1000 DH[:,3]/1000 DH[:,4]]
		end
		Jacobe(qs,dh);
	end
	ikine(Te,q0=rand(1,6),Tol=1e-9) = ikine_iter(Te,DH,q0,Tol)
end