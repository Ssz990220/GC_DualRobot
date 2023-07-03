## 
conf = ((max_iter = 200, μ = [1,1e-5,1,1,1],m=3, η=0.1, err = "NAN", stop = true,τ=0.3,reg=true,svd=false))
confG3 = ((max_iter = 200, m = 3, η=0.3,μ=0.1,scale=1e-5, err = "NAN", stop = true, reg = true, τ=0.3))
confₘ₂ = ((max_iter = 200, m = 3, η=0.1,μ=0.1,scale=1, err = "NAN", stop = true, reg = true, τ=0.3))		# GA-Katyusha
confMa = ((N=N,))
## Prepare XYZ
X = (@SMatrix [1. 0 0 0;0 1 0 0;0 0 1 270;0 0 0 1])*tfrotx(π/3)*tfroty(π/6);
Z = tfxyz([0.,0.,100])*tfrotx(π)*tfroty(0.1)
Y = tfxyz([1700.,0,200])*tfrotz(π)*tfroty(-π/6)*tfrotx(π/12);
if UNIT == "mm"
else
    X = tometer(X); Y = tometer(Y); Z = tometer(Z)
end
Ans = [X,Y,Z];
Xₘ = tfm2motorS(X); Yₘ = tfm2motorS(Y); Zₘ = tfm2motorS(Z);

## simulation Config
n_angles = [0.5,1.2,1.8]./180.0.*π  # in radius
n_dises = [0.5,1.2,1.8]./1000   #in mm
n_sample = collect(10:8:98).*3
DIM_ERROR = 8

n_sampleMa = collect(10:8:98);      # Actual Dataset Size is n_sampleMa*N 
confMa = ((N = N,))

## Sim Data Generator
function get_data(n,n_angle,n_dis)
    fix=true
    A,B,C,Q₁,Q₂ = get_ABC_noiseSₘ(n,Xₘ,Yₘ,Zₘ,n_angle,n_dis,UNIT,fix,NOISE_DISTRIBUTION)
    A_axxb, B_axxb = axxb_movementsₛ(n,Xₘ,Yₘ,Zₘ,UNIT)
    A_axxb = add_noiseSₘ.(A_axxb,n_angle,n_dis,fix,NOISE_DISTRIBUTION);
    B_axxb = add_noiseSₘ.(B_axxb,n_angle,n_dis,fix,NOISE_DISTRIBUTION);
    return ((P = [A,B,C],Q = [Q₁,Q₂],axxb=[A_axxb,B_axxb]))
end

function get_dataMa(n,n_angle,n_dis)
    fix=true
    A,B,C,_,_ = get_N_mvg_dataSₘ(N,n,n_angle,n_dis,Xₘ,Yₘ,Zₘ,true,UNIT,NOISE_DISTRIBUTION);
    return ((P = [A,B,C],))
end

function get_dataMaOne(n,n_angles,n_dis)
    fix = true
    A,B,C,_,_ = get_mvg_dataSₘ(n,n_angles, n_dis, Xₘ, Yₘ, Zₘ, fix,UNIT,NOISE_DISTRIBUTION);
    return ((P = [A,B,C],))
end


## Error Evaluator
err_func(X,Y,Z,ValidationData) = error_func(X,Y,Z,ValidationData,Ans)