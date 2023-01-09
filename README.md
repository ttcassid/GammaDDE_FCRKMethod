# GammaDDE_FCRKMethod
README:
Code from the paper: 

Numerical methods and hypoexponential approximations for gamma distributed delay differential equations, IMA Journal of Applied Mathematics, 2022; https://doi.org/10.1093/imamat/hxac027, 
Tyler Cassidy, Peter Gillich, Antony R Humphries, Christiaan H van Dorp. 

Contact Tyler Cassidy for article access and questions

These Matlab scripts replicate the convergence plot and simulation comparisions between the Functional Continuous Runge-Kutta (FCRK) method, hypoexponential approximation, and Erlang Approximation for a generic Gamma distributed DDE. The code is written in Matlab 2019b. 
 
Code for Convergence plots:

FCRKMethod_ConvergencePlots.m: Main driver file to produce the log-log plots showing the convergence rate of the FCRK method as a function of the fixed FCRK method stepsize. 

odefunNL.m & odefun: Helper functions to solve the equivalent ODE problem using the built-in ODE solvers in Matlab to provide the reference solutions when the underlying delay distribution is Erlang. The equivalent systems are derived using the Linear Chain Technique. 

ConvergencePlots_SimulationDynamics.m: Simulation script to plot the time-evolution dynamics of the system solved in the convergence plots; outputs a figure of of the solution of the DDE (y(t)) vs simulation time.

Code for FCRK Method simulations:

HypoODEVsGammaDDE_SimulationScript.m: This file computes a visual comparison of the solution of a generic Gamma distributed DDE calculated using the FCRK method or two
approximations. The two approximations are either the Hypoexponential approximation or the Erlang approximation derived in the main text. 

FCRK Method files:

ddef4.m: The implementation of the 4th order fixed stepsize FCRK method. This function is similar to ode45/dde23 Matlab solvers but adapted to the gamma distributed DDE.

herm3terp.m : Hermite interpolation used as a helper function in the FCRK method for the stage interpolation step

linTerp.m : Linear interpolation used as a helper function in the FCRK method for the stage interpolation step

