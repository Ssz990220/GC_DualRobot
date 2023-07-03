__precompile__()
module DualRobotCalibrator
    using Reexport

    include("../utils/GARobotics.jl")
    @reexport using .GARobotics

    include("utils.jl")
    export add_noiseSₘ, get_error, get_error_result, get_rot_error, tometer, tommeter, get_report

    include("SVRG.jl")
    export katyushaX, SVRG, SVRG_BB

    include("../utils/UR10.jl")

    include("AXBYCZ_utils.jl")
    export get_ABC, get_ABCSₘ, get_ABC_noiseSₘ, match_signS, get_error_XYZ, T_SVD

    include("AXBYCZ_Fu_Utils.jl")
    export AXB_YCZ_FUS

    include("G3_Close.jl")
    export match_sign_ABC, get_XYZ_rotor_Close
    
    include("AXBYCZ_Liao_Utils.jl")
    export AXBYCZ_Liao, liao_close, liao_iter

    include("AXBYCZ_Wang_Utils.jl")
    export AXB_YCZ_Wang, AXB_YCZ_Close

    include("AXXB_utils.jl")
    export AXXBs, axxb_movements_noiseₛ, axxb_movementsₛ

    include("G3_Calculus.jl")
    export AXBYCZ_G3, iter_solve_Katyusha_G3

    include("AXBYCZ_Ma.jl")
    export axbyczProb1, axbyczProb3, meanCov, get_N_mvg_dataSₘ, AXBYCZ_Wang2, get_mvg_dataSₘ

    ## GC
    include("GA_Calculus.jl")
    export iter_solve_GC_katyushaᵥ,get_i_grad_reg_AXBYCZᵥ₂,get_full_grad_regᵥ₂,normfunc2ᵥ

    include("GA_Close.jl")
    export GAClose

    # Benchmark Functions and Wrappers
    include("Benchmark_AXBYCZ.jl")
    export error_func, liao, G3, wang_cic, fu_bayro, fu_bayro13, Ma, benchmark_p, MaOne, GA_Full

    # Experiment Data IO
    include("../Experiments/pre_process.jl")
    export read_Data, clean_up, load_dataS, unselect_dataset, load_data
end