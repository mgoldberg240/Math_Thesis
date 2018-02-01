# Math_Thesis
This repository includes MATLAB code used to study Numerical Methods for solving PDEs.
As of Feb 1, codes detail solutions to the 1D Diffusion Equation.  Future codes will be applied towards other problems and will feature in depth analysis and comparisons of truncation errors as well as accuracy.

diffusionFTCS.m
---------------
-Solve the 1D Diffusion Equation using the Forward in Time, Central in Space (FTCS) method
-Outputs satisfaction of stability condition
-Generates plot of approach from initial conditions to final solution
-Verified using analytical solution
-If stability condition is not satisfied, plot demonstrates numerical ringing effects

diffusionCN.m
-------------
-Solve the 1D Diffusion Equation using the Crank Nicholson method
-Employs recursive solution to tridiagonal matrix system (see tridiag.m)

tridiag.m
---------
-Construct and solve a tridiagonal matrix system Au=b
