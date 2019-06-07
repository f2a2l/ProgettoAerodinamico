%% Particle Swarm Optimization Procedure for F1 Rear Wing Design


optiPAR; %as output: lb and ub. Here it is possibile to define boundaries for optimization
fun = @lapTime_objFun;
nvars = 21;

options = optimoptions('particleswarm','SwarmSize',50,'Display','Iter','FunctionTolerance',1e-6,...
                        'MaxIterations',200*nvars,'MaxTime',60*60*24,'PlotFcn','pswplotbestf',...
                        'UseParallel',1);
                    
[optimalSolution] = particleswarm(fun,nvars,lb,ub);