%% Parameters

%clear all
%clc
%close all

%c1 = 0.3; %coeff camber-line-abscissa parameter equation [0.01,0.96]
%c2 = 0.6; %coeff of camber-line-abscissa parameter equation [0.02,0.97]
%c3 = 0; %coeff of camber-line-ordinate parameter equation [-0.074,0.247]
%c4 = 0; %coeff of camber-line-ordinate parameter equation [-0.102,0.206]
%xt = 0.3; %chordwise location of maximum thickness [0.2002,0.4813]
%T = 0.12; %max thickness [0.0246,0.3227]
%%%%%rho0 = 1;
%%%%%rho0_bar = rho0 / (T/xt)^2; % dimensionless rho0 [0.175,1.4944]
%rho0_bar = 0.3;
%%%%beta_te = 10; %trailing edge boat-tail angle []
%%%%beta_te_bar = beta_te/(atan((T)/(1-xt))); %dimless beta_te [0.1452,4.8724]
%beta_te_bar = 1.5;


%param = [c1 c2 c3 c4 xt T rho0_bar,beta_te_bar];


%[x,y] = AirfoilShape(param);

%plot(x,y)
%xlim([0,1]) 
%axis equal

function [x,y] = AirfoilShape(param, npoint)
%AirfoilShape - Genera profili con IGP


    c1 = param(1);
    c2 = param(2);
    c3 = param(3);
    c4 = param(4);
    xt = param(5);
    T = param(6);
    rho0_bar = param(7);
    beta_te_bar = param(8);

    [t1, t2, t3, t4, t5] = getThickParam(xt,T,rho0_bar,beta_te_bar);

    xu = zeros(npoint,1);
    yu = xu;
    xl = xu;
    yl = xu;

    for i = 1 : npoint

        % linear - deprecated       k = (i - 1)/(npoint-1);
        k = 1 - 0.5 * (1 + cos(((i - 1) * pi) / (npoint - 1)));
        
        [xc,yc] = camberline(c1,c2,c3,c4,k);

        thickfun = @(x) t1 * sqrt(x) + t2 * x + t3 * x^2 + ...
                                           t4 * x^3 + t5 * x^4;
        t = thickfun(k);

        xu(i) = xc;
        yu(i) = yc + 0.5 * t;

        xl(i) = xc;
        yl(i) = yc - 0.5 * t;

        
    end

    x = [flip(xl); xu(2:end)];
    y = [flip(yl); yu(2:end)];

end
