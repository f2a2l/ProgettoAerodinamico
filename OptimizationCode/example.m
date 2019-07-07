fun = @dejong5fcn;
nvars = 2;
rng default % For reproducibility
lb = [-50;-50];
ub = -lb;
options = optimoptions('particleswarm','SwarmSize',5000,'PlotFcn','pswplotbestf');
[x,fval,exitflag] = particleswarm(fun,nvars,lb,ub,options)
