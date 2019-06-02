fun = @dejong5fcn;
nvars = 2;
rng default % For reproducibility
lb = [-50;-50];
ub = -lb;
options = optimoptions('particleswarm','SwarmSize',100);
[x,fval,exitflag] = particleswarm(fun,nvars,lb,ub,options)
