include("../utils/GARobotics.jl")
using .GARobotics, StaticArrays, Random, Test, LinearAlgebra
include("../utils/UR10.jl")

@info "Random sample in configuration space of UR10 robot."
q = @SVector rand(6);
T = UR10.fkine(q)
RG3,tG3 = UR10.fkineG3S(q)
TG3 = rt2T(rotorG32rm(RG3),tG3)
UR10.fkineGAS(q)
TGA = motorS2tfm(UR10.fkineGAS(q))
@info "Test G3 Forward Kinematics"
@show @test norm(TG3 - T) <1e-10

@info "Test PGA Forward Kinematics"
@show @test norm(TGA - T) <1e-10