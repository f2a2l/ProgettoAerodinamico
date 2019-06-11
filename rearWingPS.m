%% Particle Swarm Optimization Procedure for F1 Rear Wing Design
clear all
close all
clc

%% Add Paths
addPaths;


%% Optimization
optiPAR; %as output: lb and ub. Here it is possibile to define boundaries for optimization

fun = @lapTime_objFun;
nvars = 21;

%Max simulation time
days = 1;
MaxTime = 60 * 60 * 24 * days;

options = optimoptions('particleswarm','SwarmSize',50,'Display','Iter',...
                       'FunctionTolerance',1e-6,'MaxIterations',200*nvars,...
                       'MaxTime',MaxTime,'PlotFcn','pswplotbestf',...
                       'UseParallel',true,'UseVectorized',false);
                                      
[optimalSolution,fval,exitflag,output] = particleswarm(fun,nvars,lb,ub,options);

%% Generate .dat file for CFD Meshing
geometryCFD;