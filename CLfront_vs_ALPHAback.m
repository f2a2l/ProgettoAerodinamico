close all
clear 
clc

%%
% Front airfoil geometry parameters
arflPar_single = [0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5];

% Tandem configuration geometry parameters
arflPar = [0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5;
           0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5];

alpha_f = -5:1:15;     % range of AoA front airfoil
alpha_b = 30 : 1 : 60; % range of AoA back airfoil


Cl_double = zeros(2,1);
Cl_front= zeros(length(alpha_f),length(alpha_b)); % front airfoil CL matrix
Cl_back = zeros(length(alpha_f),length(alpha_b)); % back airfoil CL matrix


%% ATTENZIONE

% scommenta per rifare i calcoli (più di 1 ora di calcoli)

% for j = 1 : length(alpha_f)
%     
%     for i = 1 : length(alpha_b)
%     
%         [Cl_double, ~, ~, ~] = solverHS(80, arflPar, [alpha_f(j), alpha_b(i)], [1.05, -0.05], 0.35);
%         Cl_front(j,i) = Cl_double(1);
%         Cl_back(j,i) = Cl_double(2);
%     
%     end 
% end

% altrimenti fai solo:     load('CLfront_vs_ALPHAback.mat')

%% Plot 

figure()
for j = 1 : length(alpha_f)
    plot(alpha_b,Cl_front(j,:),'-o','DisplayName',['\alpha_{front} = ',num2str(alpha_f(j)),'°'])
    hold on
    grid on
end
legend('show','Location','NW')
title ('Variazione del CL_{front} al variare dell'' \alpha_{back}')
xlabel ('\alpha_{back}')
ylabel ('CL_{front}')

%%

figure()
for i = 1 : length(alpha_b)
    plot(alpha_f,Cl_front(:,i),'-o','DisplayName',['\alpha_{back} = ',num2str(alpha_b(i)),'°'])
    hold on
    grid on
end
legend('show','Location','NW')
xlabel ('\alpha_{front}')
ylabel ('CL_{front}')

figure()
for i = 1 : length(alpha_f)
    plot(alpha_b,Cl_back(i,:),'-o','DisplayName',['\alpha_{front} = ',num2str(alpha_f(i)),'°'])
    hold on
    grid on
end
legend('show','Location','SW')
xlabel ('\alpha_{back}')
ylabel ('CL_{back}')


%% Approssimazione tramite parabola del CL FRONT

% Si seleziona un angolo front: es. alpha_b = 40° --> alpha_b(21)
index_back = 12;

XX = [ones(length(alpha_f),1), alpha_f', (alpha_f.^2)'];
YY = Cl_front(:,index_back);

bb = XX\YY;

Cl_front_approx = (bb(1)+bb(2).*alpha_f+bb(3).*alpha_f.^2)';

figure()
subplot(2,1,1)
plot(alpha_f,YY,'LineWidth',1)
hold on
plot(alpha_f,Cl_front_approx,'Linewidth',1)
grid on
xlabel('\alpha_{front}')
ylabel ('CL_{front}')
title (['Valutazione per \alpha_{back} = ',num2str(alpha_b(index_front)),'°'])

subplot(2,1,2)
plot(alpha_f,abs(YY-Cl_front_approx),'-o')
grid on
xlabel('\alpha_{front}')
ylabel 'errore'

B = zeros(3,length(alpha_b));
for j = 1 : length(alpha_b)
    
    B(:,j) = parabola(alpha_f,Cl_front(:,j));
    
end

figure()
subplot(3,1,1)
plot(alpha_b,B(1,:))
ylabel('b_0')
xlabel('\alpha_{back}')
grid on

subplot(3,1,2)
plot(alpha_b,B(2,:))
ylabel('b_1')
xlabel('\alpha_{back}')

grid on

subplot(3,1,3)
plot(alpha_b,B(3,:))
ylabel('b_2')
xlabel('\alpha_{back}')

grid on

%% Approssimazione tramite parabola del CL BACK
% Si seleziona un angolo front: es. alpha_f = 5° --> alpha_f(11)
index_front = 16;

XX = [ones(length(alpha_b),1), alpha_b', (alpha_b.^2)'];
YY = Cl_back(index_front,:)';

bb = XX\YY;

Cl_back_approx = (bb(1)+bb(2).*alpha_b+bb(3).*alpha_b.^2)';

figure()
subplot(2,1,1)
plot(alpha_b,YY,'LineWidth',1)
hold on
plot(alpha_b,Cl_back_approx,'Linewidth',1)
grid on
xlabel('\alpha_{back}')
ylabel ('CL_{back}')
title (['Valutazione per \alpha_{front} = ',num2str(alpha_f(index_front)),'°'])

subplot(2,1,2)
plot(alpha_b,abs(YY-Cl_back_approx),'-o')
grid on
xlabel('\alpha_{back}')
ylabel 'errore'



%%%%%%%%%%% VECCHIA ROBA %%%%%%%%%%%%%%%%

%% Linear regression of data

b = zeros(2,length(alpha_b));

figure()
for i = 1 : length(alpha_b)
    
    Y = Cl_front(:,i);
    X = [ones(length(alpha_f),1),alpha_f'];
    b(:,i) = X\Y;
    
    plot(alpha_f,alpha_f*b(2,i)+b(1,i),'DisplayName',['\alpha_{back} = ',num2str(alpha_b(i)),'°'])
    hold on
    grid on
end
legend('show','Location','NW')

figure()
subplot(2,1,1)
plot(alpha_b,b(1,:))
grid on
xlabel('\alpha_{back}')
ylabel('q')
title('C_{L,front} = m \cdot \alpha_{front} + q')

subplot(2,1,2)
plot(alpha_b,b(2,:))
grid on
xlabel('\alpha_{back}')
ylabel('m')






