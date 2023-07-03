# DualRobotCalibrator
Algorithms for calibrating dual robot by solving $\mathbf{AXB=YCZ}$ problem.

* **DualRobotCalibrator**: Module Exporter. Import this file and `using .DualRobotCalibrator` to use this package.
    * **AXBYCZ_Fu_Utils**: reproduction of Fu's [A Dual Quaternion-Based Approach for Coordinate Calibration of Dual Robots in Collaborative Motion](https://ieeexplore.ieee.org/document/9072583).
    * **AXBYCZ_Liao_Utils**: Reproduction of Wu's [Simultaneous Hand–Eye, Tool–Flange, and Robot–Robot Calibration for Comanipulation by Solving the AXB=YCZ Problem](https://ieeexplore.ieee.org/document/7425217)
    * **AXBYCZ_Ma**: Reproduction of Ma's [Probabilistic approaches to the AXB= YCZ calibration problem in multi-robot systems](https://rpk.lcsr.jhu.edu/wp-content/uploads/2018/05/axb_auro_publish.pdf)
    * **AXBYCZ_Wang**: Reproduction of Wang's [Simultaneous calibration of multicoordinates for a dual-robot system by solving the AXB = YCZ problem](https://ieeexplore.ieee.org/document/9319563)
    * **AXXB_utils**: Reproduction of Bayro's [Motor Algebra for 3D Kinematics: The Case of the Hand-Eye Calibration](https://link.springer.com/article/10.1023/A:1026567812984)
    * **Benchmark_AXBYCZ**: BenchmarkUtils that wraps functions for benchmark
    * **G3_Calculus**: Proposed method that solves $\mathbf{AXB=YCZ}$ with $\mathbb{G}^3$.
    * **G3_Close**: Script that solve $\mathbf{AXB=YCZ}$ closed from and correct the signs of $\mathbf{A,B,C}$
    * **GA_Calculus**: Unpublished proposed method that solves $\mathbf{AXB=YCZ}$ with PGA and GC.
    * **GA_Close**: Unpublished closed-form method that solve $\mathbf{AXB=YCZ}$ closed from and correct the signs of $\mathbf{A,B,C}$ with PGA.
    * **SVRG**: SVRG based optimizors included SVRG, SVRG-BB, KatyushaX.
    * **utils**: some utility functions.