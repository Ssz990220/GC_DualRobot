# utils
Most fundamental utilites for this project.

* **GARobotics**: Module Exporter. Import this file and `using .GARobotics` to use the defined functions.
    * **GA.jl**: Basic functions for **Geometric Algebra**.
            *NOTED* that both PGA & G3 algebra are included. G3 is a subalgebra of PGA.
    * **GA_robotics.jl**: Basic functions for **Robot kinematics with GA & G3**.
    *NOTED* that forward kinematics is implemented in PGA & G3.
    * **Robotics.jl**: Basic functions for **Robotics**. Forward kinematics based on DH
    * **Alge.jl**: Modified algebric functions from MATLAB 

## Robot Models
* TemplatRobot
* UR10
* ABB 4600-60-2.05
* JAKA Zu 12
* Puma 560