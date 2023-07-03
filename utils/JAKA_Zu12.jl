export JAKA
module JAKA
	using StaticArrays
	using Rotations
	include("Alge.jl")
	include("GA.jl")
	include("Robotics.jl")
	include("GA_robotics.jl")
				 # theta a(x)  	 d(z)	alpha(x)
	DH = @SMatrix [0.0   0.0     140.6   π/2;
				   0.0   595.0   0.0     0.0;
				   0.0   574.5   0.0     0.0;
				   0.0   0.0     -128.5  π/2;
				   0.0   0.0     112.0   -π/2;
				   0.0   0.0     100.5   0.0];
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