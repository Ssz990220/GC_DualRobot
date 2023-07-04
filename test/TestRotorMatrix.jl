## Test if Rand Rotor and Rand Matrix works identically
using LinearAlgebra, StaticArrays, Rotations
using Distributions
using Test
include("../utils/Alge.jl")
include("../utils/Robotics.jl")
include("../utils/GA.jl")
include("../utils/GA_robotics.jl")
## Read rand_rotorS in GA.jl for random rotation generation


axis = @SVector rand(3); axis = axis ./ norm(axis)
θ = rand(Uniform(-π,π));
###################
## Test Rotation ##
###################
RotationError = motorS2tfm(rotorS(θ,axis))[SOneTo(3), SOneTo(3)] - AngleAxis(θ,axis...)
@info "Rotation Test"
@test norm(RotationError)<1e-15

######################
## Test Translation ##
######################
dis = rand(Uniform(-10,10));
TranslationError = motorS2tfm(translatorS(dis*axis)) - t2T(dis*axis)
@info "Translation Error Test"
@test norm(TranslationError)<1e-15

#################
## Test Hybrid ##
#################
TrueValue = rt2T(rand(RotMatrix{3}),rand(3)*3)
TrueValueMotor = tfm2motorS(TrueValue);
begin
    axisR1 = @SVector rand(3);  axisR1 = axisR1 ./ norm(axisR1)
    axisR2 = @SVector rand(3);  axisR2 = axisR2 ./ norm(axisR2)
    θ₁ = rand(Uniform(-π,π));   θ₂ = rand(Uniform(-π,π));
    axist1 = @SVector rand(3);  axist1 = axist1 ./ norm(axist1)
    axist2 = @SVector rand(3);  axist2 = axist2 ./ norm(axist2)
    dis₁ = rand(Uniform(-10,10));dis₂ = rand(Uniform(-10,10));

    rotor1 = rotorS(θ₁,axisR1); rotor2 = rotorS(θ₂,axisR2); 
    translator1 = translatorS(dis₁*axist1);
    translator2 = translatorS(dis₂*axist2);
    R1 = r2t(AngleAxis(θ₁,axisR1...)); R2 = r2t(AngleAxis(θ₂,axisR2...));
    t1 = t2T(dis₁*axist1); t2 = t2T(dis₂*axist2);
end
MotorHybrid = ecmul(ecmul(rotor1,TrueValueMotor),translator1); motorS2tfm(MotorHybrid)
MotorHybrid2 = ecmul(ecmul(rotor2,MotorHybrid),translator2); motorS2tfm(MotorHybrid2)
@info "Hybrid Rotation Translation Test"
@test norm(t1*TrueValue*R1 - motorS2tfm(MotorHybrid))<1e-14
@test norm(t2*t1*TrueValue*R1*R2 - motorS2tfm(MotorHybrid2))<1e-14