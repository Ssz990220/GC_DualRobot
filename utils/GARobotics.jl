__precompile__()
module GARobotics
    using Reexport
    @reexport using StaticArrays
    @reexport using LinearAlgebra
    @reexport using Rotations
    @reexport using Random
    @reexport using Rotations
    @reexport using Statistics
    @reexport using Distributions

    include("Alge.jl")
    include("GA.jl")
    include("Robotics.jl")
    include("GA_robotics.jl")

## Alge.jl
export eye, bitget, NearZero, trace, Sdiag, cvtSMatrix

## GA.jl
export rotorS, translatorS, elineS, emomentS, motorS2tₛ, Qᵣ, Qₗ, approx, quat2rotor, 
       rotorS2quat, rotorG32quat, rotorS2vec, drotorS2quat, motorS2rotorS, motordS2rotorS,
       Mₗ, Mᵣ, rand_rotorS, rand_translatorS, ecmul, ereversion, eMₗ, motorS2rotorG3S, rotorG3S2motorS, ecmulG3,
       ereversionG3, Rᵣ, Rₗ

## Robotics.jl
export FkineDH, tfrotx, tfroty, tfrotz, rotx, roty, rotz, tfxyz, mcross, dmcross, invT, tfm2R,
       SE3Adjoint, SE3AdjointInv, ProjectToSO3, ProjectToSE3, wraptoπ, wrapto1, MatrixExp3, MatrixExp6, MatrixLog6, T2rt,
       MatrixLog3, rm2Quat, vex, vee, Quat2rm, t2r, t2T, r2t, rt2T, angvec2r

## GA_robotics.jl
export scale_motorS, tfm2motorS, motorS2tfm, rotorS2rm, rotorG32rm, rm2rotorS, MᵢⱼS, FkineGAS, rotorG3Norm,
       motorNormS, rm2rotorG3S

end