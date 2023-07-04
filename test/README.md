# test
Scripts in this folder verifies some functions in the project.

Kinematics of this package is represented in both Transformation Matrices (primarily), $\mathbb{G}^3$ (partially, for completeness) and PGA (mostly). 
The correctness of the kinematics is verified by comparing the results from these two representations.

## Tests
* `TestFkine.jl`: Equivalence of the forward kinematics in PGA and Transformation Matrices.
* `TestRotorMatrix.jl`: Equivalence of the rotation matrix and rotor in $\mathbb{G}^3$.
* `TestMinimal.jl`: Test the minimal number of samples required by each algorithm.