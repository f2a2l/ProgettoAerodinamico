%% Particle Swarm Optimization Procedure for F1 Rear Wing Design
clear all
close all
clc

%% Add Paths
addPaths;


%% Optimization
optiPAR_v5; %as output: lb and ub. Here it is possibile to define boundaries for optimization

fun = @lapTime_objFun;
nvars = 21;

%Starts ParallelPool for Parallel Computation
parpool;

%Max simulation time
days = 7;
MaxTime = 60 * 60 * 24 * days;

options = optimoptions('particleswarm','SwarmSize',50,'Display','Iter',...
                       'FunctionTolerance',1e-6,'MaxIterations',200*nvars,...
                       'MaxTime',MaxTime,'PlotFcn','pswplotbestf',...
                       'UseParallel',true,'UseVectorized',false);
                                      
[optimalSolution,fval,exitflag,output] = particleswarm(fun,nvars,lb,ub,options);

%% Save data and generate .dat files for CFD Meshing
saveProfiles;
geometryCFD;
save workspace;
