export UR10
module UR10
	using StaticArrays
	using Rotations
	include("Alge.jl")
	include("GA.jl")
	include("Robotics.jl")
	include("GA_robotics.jl")
	 		  	   #theta 	a(x)  	d(z)		alpha(x)
	DH = @SMatrix [0 		0 		127.3 		π/2;
			       0 		-612 	0 			0;
			       0 		-572.3 	0 			0;
			       0 		0 		163.941 	π/2;
			       0 		0 		115.7 		-π/2
			       0 		0 		92.2 		0]
	function fkine(qs,unit="mm")
		if unit=="mm"
			dh = DH;
		else
			dh = [DH[:,1] DH[:,2]/1000 DH[:,3]/1000 DH[:,4]]
		end
		FkineDH(qs,dh);
	end
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
