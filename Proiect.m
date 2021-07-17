%% Date initiale
clc
clear all
format short
R = 2e6
L = 1e4
m = 0.1
K = 10
D = 10
epsilon = 10e-6 * 8.854
A = 2* 10e-2
alpha = epsilon * A
beta = 0.005

%% Date liniarizare
E = 0;
Fext = 0;
%% Spatiul starlior

syms x1 x2 x3 x4

x1p = x2;
x2p = 1/m * (Fext - D * x2 - K*x1 - x3^2/(2*alpha));
x3p = x4;
x4p = 1/L * (E - R*x4 - x1*x3/alpha);


Solutie = solve([x1p == 0, x2p == 0, x3p == 0, x4p ==0]);
Solutie.x1
Solutie.x2
Solutie.x3
Solutie.x4
%% Spatiul starilor
A = [0 1 0 0; -K/m -D/m 0 0;0 0 0 1;0 0 0 -R/L];
B = [0 0;0 1/m;0 0;1/L 0];
C = [1 0 0 0;0 0 1 0];
D1 = [0 0;0 0];
%%
Con = (ctrb(A,B));
rank(Con)
% K_control = ((place(A,B,[-1 -80 -0.5 -150])));
K_control = (place(A,B,[-1 -10 -0.5 -20]))
% K_control = lqr(A,B,eye(4),eye(2));
Obs = obsv(A,C);
rank(Obs)
% K_es =((place(A',C',[-1 -80 -0.5 -150])));
K_es =((place(A',C',[-1 -10 -0.5 -20]*5)));
Fx = inv(C*inv(-A + B*K_control)*B)
%%
Te = (1/200)/10;
sysC = ss(A,B,C,D1)
sysD = c2d(sysC,Te,'zoh')
%% Varianta 2 -Euler
Ad2 = eye(4) +Te*A;
Bd2 = Te*B
con2 = ctrb(Ad2,Bd2)
rank(con2)
pCont = [-1 -10 -0.5 -20]
pDis = exp(pCont.*Te)
K_control_dis2 = (place(Ad2,Bd2,pDis))
obs2 = obsv(Ad2,C)
rank(obs2)
% pCont_ds2 = [-1 -10 -0.5 -20]*5
% pDis_des2 = exp(pCont_ds2.*Te)
pDis_des2 = [0.5 0.7 0.6 0.8];
K_es_dis2 =((place(Ad2',C',pDis_des2)))
Fx_dis2 = inv(C*inv(-Ad2 + Bd2*K_control_dis2)*Bd2)
%% Liniar
x0 = [-1 -0.5 0.5 1];
x(:,1) = x0;

x0_es = [0 0 0 0];
x_es(:,1) = x0_es;

Te = 1/2000;
K_timp = 2/Te;


u1(1,1) = E;
u1(2,1) = Fext;

for k = 1:K_timp
    
    dx = [x(2,k),...
        (u1(2,k)-D*x(2,k) - K*x(1,k))/m,...
        x(4,k),...
        (u1(1,k)-R*x(4,k))/L]';
    
    x(:,k+1) = x(:,k)+Te*dx;
    u1(:,k+1) = -K_control_dis2 * x(:,k+1);
    
    dx_es = [x_es(2,k),...
        (u1(2,k)-D*x_es(2,k) - K*x_es(1,k))/m,...
        x_es(4,k),...
        (u1(1,k)-R*x_es(4,k))/L]';
    
    x_es(:,k+1) = x_es(:,k)+Te*dx_es + K_es_dis2'*C*(x(:,k) - x_es(:,k));
    
    
end

plot(x')
legend('x1','x2','x3','x4');
figure
plot(x_es')
title('Estimat');
legend('x1','x2','x3','x4');
figure
plot(x' - x_es')
title('Eroare');
legend('x1','x2','x3','x4');
%% Neliniar discret
clear x
clear x_es
x0 = [-0.01 -0.005 0.005 0.01];
x(:,1) = x0;

x0_es = [0 0 0 0];
x_es(:,1) = x0_es;

Te = 1/2000;
K_timp = 0.1/Te;

K_control_dis3 = (place(Ad2,Bd2,[0.8 0.9 0.7 0.6]))
u1(1,1) = E;
u1(2,1) = Fext;

for k = 1:K_timp
    
    dx = [x(2,k),...
        (u1(2,k)-D*x(2,k) - K*x(1,k) - x(3,k)^2/(2*alpha))/m,...
        x(4,k),...
        (u1(1,k)-R*x(4,k)-x(1,k)*x(3,k)/alpha)/L]';
    
    x(:,k+1) = x(:,k)+Te*dx;
    
    
    dx_es = [x_es(2,k),...
        (u1(2,k)-D*x_es(2,k) - K*x_es(1,k) - x_es(3,k)^2/(2*alpha))/m,...
        x_es(4,k),...
        (u1(1,k)-R*x_es(4,k)-x_es(1,k)*x_es(3,k)/alpha)/L]';
    
    x_es(:,k+1) = x_es(:,k)+Te*dx_es + K_es_dis2'*C*(x(:,k) - x_es(:,k));
    u1(:,k+1) = -K_control_dis3 * x_es(:,k+1);
    
end

plot(x')
legend('x1','x2','x3','x4');
figure
plot(x_es')
title('Estimat');
legend('x1','x2','x3','x4');
figure
plot(x' - x_es')
title('Eroare');
legend('x1','x2','x3','x4');
%% Kalman
clear x
clear x_es
clear dx
clear dx_es

Q = diag([1e-3,1e-3,10e-4,10e-4])*10e-4;
R_Kal = diag([10e-6 10e-7])
P = eye(4) * 1e11;


x0 = [-1 -0.5 0.5 1];
x(:,1) = x0;

x0_es = [0 0 0 0];
x_es(:,1) = x0_es;

Te = 1/2000;
K_timp = 2/Te;


u1(1,1) = E;
u1(2,1) = Fext;

for k = 1:K_timp
    
    dx = [x(2,k),...
        (u1(2,k)-D*x(2,k) - K*x(1,k))/m,...
        x(4,k),...
        (u1(1,k)-R*x(4,k))/L]';
    
    x(:,k+1) = x(:,k)+Te*dx + cholcov(Q)*randn(4,1);
    y(:,k) = C*x(:,k)+ cholcov(R_Kal)*randn(2,1);
    
    u1(:,k+1) = -K_control_dis2 * x(:,k+1);
    
    dx_es = [x_es(2,k),...
        (u1(2,k)-D*x_es(2,k) - K*x_es(1,k))/m,...
        x_es(4,k),...
        (u1(1,k)-R*x_es(4,k))/L]';
    
    x_pred = x_es(:,k)+Te*dx_es;
    P_pred = Ad2*P*Ad2' + Q;
    K_es_dis4 = P_pred*C'*inv(C*P_pred*C'+R_Kal);
    x_es(:,k+1) = x_pred + K_es_dis4*(y(:,k) - C*x_pred);
    P = (eye(4) - K_es_dis4*C)*P_pred*(eye(4)-K_es_dis4*C)'+K_es_dis4*R_Kal*K_es_dis4';
    
    
    
end
plot(x')
legend('x1','x2','x3','x4');
figure
plot(x_es')
title('Estimat');
legend('x1','x2','x3','x4');
figure
plot(x' - x_es')
title('Eroare');
legend('x1','x2','x3','x4');
%% Kalman extins
% Punem dx de la cel neliniar
% Recalculam Ad2 la pasul k
clear x
clear x_es
clear dx
clear dx_es

x0 = [-0.01 -0.005 0.005 0.01];
x(:,1) = x0;

x0_es = [0 0 0 0];
x_es(:,1) = x0_es;

Te = 1/2000;
K_timp = 0.1/Te;

Q = diag([1e-3,1e-4,10e-3,10e-4])*10e-3;
R_Kal = diag([10e-10 10e-12])*1e3
P = eye(4) * 1e8;


K_control_dis3 = (place(Ad2,Bd2,[0.8 0.9 0.7 0.6]))
u1(1,1) = E;
u1(2,1) = Fext;
% Fara estimator
for k = 1:K_timp
    
    dx = [x(2,k),...
        (u1(2,k)-D*x(2,k) - K*x(1,k) - x(3,k)^2/(2*alpha))/m,...
        x(4,k),...
        (u1(1,k)-R*x(4,k)-x(1,k)*x(3,k)/alpha)/L]';
    
    
    x(:,k+1) = x(:,k)+Te*dx + cholcov(Q)*randn(4,1);
    y(:,k) = C*x(:,k)+ cholcov(R_Kal)*randn(2,1);
    
    dx_es = [x_es(2,k),...
        (u1(2,k)-D*x_es(2,k) - K*x_es(1,k) - x_es(3,k)^2/(2*alpha))/m,...
        x_es(4,k),...
        (u1(1,k)-R*x_es(4,k)-x_es(1,k)*x_es(3,k)/alpha)/L]';
    
    x_pred = x_es(:,k)+Te*dx_es;
    
    A_es = Te*[0,1,0,0; ...
        -K/m,-D/m,-x_es(3,k)/(alpha*m),0;...
        0,0,0,1;...
        -x_es(3,k)/(alpha*L),0,-x_es(1,k)/(alpha*L),-R/L];
    B_es =[0,0;0,1;0,0;1,0];
    
    P_pred = A_es*P*A_es' + Q;
    
    
    K_es_dis4 = P_pred*C'*inv(C*P_pred*C'+R_Kal);
    x_es(:,k+1) = x_pred + K_es_dis4*(y(:,k) - C*x_pred);
    P = (eye(4) - K_es_dis4*C)*P_pred*(eye(4)-K_es_dis4*C)'+K_es_dis4*R_Kal*K_es_dis4';
    
    u1(:,k+1) = -K_control_dis3 * x_es(:,k+1);
    
end

plot(x')
legend('x1','x2','x3','x4');
title('Starea')
figure
plot(x_es')
title('Estimat');
legend('x1','x2','x3','x4');
figure
plot(x' - x_es')
title('Eroare');
legend('x1','x2','x3','x4');
%% Estimarea unei intrari necunoscute
close all
clear x
clear x_es

d =1;
E_nec = [0;d;0;0;0];
A_nec = [[A;0,0,0,1],E_nec];
B_nec = [B;0 0]
C_nec = [C,[0;0]];

Ad2_nec = eye(5) +Te*A_nec;
Bd2_nec = Te*B_nec;

con2 = ctrb(Ad2_nec,Bd2_nec)
rank(con2)

% pCont = [-1 -10 -0.5 -20]
% pDis = exp(pCont.*Te)
% K_control_dis2 = (place(Ad2_nec,Bd2_nec,pDis))
obs2 = obsv(Ad2_nec,C_nec)
rank(obs2)
pDis_des_nec = [0.2 0.8 0.5 0.3 0.1];
K_es_dis_nec =((place(Ad2_nec',C_nec',pDis_des_nec)))

x0 = [-0.1 +0.05 -0.5 1 10];
x(:,1) = x0;

x0_es = [0 0 0 0 0];
x_es(:,1) = x0_es;

Te = 1/2000;
K_timp = 0.1/Te;


u1(1,1) = E;
u1(2,1) = Fext;

for k = 1:K_timp
    
    %Ax+BU
    dx_nec =[x(2,k),...
        (u1(2,k)-D*x(2,k) - K*x(1,k))/m,...
        x(4,k),...
        (u1(1,k)-R*x(4,k))/L,...
        0]';
    %Ax+Bu+Ed
    dx_nec = dx_nec + E_nec*x(5,k);
    
    x(:,k+1) = x(:,k)+Te*dx_nec;
    y(:,k) = C_nec*x(:,k);
    
    dx_es_nec = [x_es(2,k),...
        (u1(2,k)-D*x_es(2,k) - K*x_es(1,k))/m,...
        x_es(4,k),...
        (u1(1,k)-R*x_es(4,k))/L,...
        0]';
     dx_es_nec = dx_es_nec + E_nec*x_es(5,k);
    y_es(:,k) = C_nec*x_es(:,k);
    x_es(:,k+1) = x_es(:,k)+Te*dx_es_nec +K_es_dis_nec'*(y(:,k) - y_es(:,k));
%     u1(:,k+1) =  u1(:,k);
     u1(:,k+1) = -K_control_dis3 * x_es(1:4,k+1);
     u1(2,k+1) = u1(2,k+1) - x_es(5,k+1)*m;
end

plot(x')
legend('x1','x2','x3','x4','d');
figure
plot(x_es')
title('Estimat');
legend('x1','x2','x3','x4','d');
figure
plot(x' - x_es')
title('Eroare');
legend('x1','x2','x3','x4','d');

%% Decuplare
clear x
clear x_es
clear z
clear d_dec
clear u1
clear C
clear y

Te = 1/2000;
K_timp = 20/Te;
C = [1 0 0 0;0 0 1 0;0 1 0 0];
E_dec = Te*[0;1;0;0]
H = E_dec*pinv(C*E_dec)
T = eye(4)-H*C;
% T=[0,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1]
A_dec=T*Ad2;
rank(obsv(A_dec,C))
% syms k11 k12 k21 k22 k31 k32 k41 k42
% K_dec = [k11 k12;k21 k22;k31 k32;k41 k42];
% vpa(A_dec-K_dec*C)
% [0.5 0.7 0.6 0.8]
% F = diag([0.4,0.5,0.3,0.1]);

K_dec =place((A_dec)',C',[0.4,0.5,0.3,0.1])';
F = A_dec-K_dec*C;
K_dec2 = F*H;

x0 = [-0.1 +0.05 -0.5 1];
x(:,1) = x0;

y(:,1) = C*x(:,1);
z0 = [0 0 0 0];
z(:,1) = z0;
x_es(:,1) = z0'+H*y(:,1);




u1(1,1) = E;
u1(2,1) = Fext;
u1(:,1) = [0,0]';
d_dec(:,1) = 0;

for k = 1:K_timp
    
    %Ax+BU
%     dx_nec =[x(2,k),...
%         (u1(2,k)-D*x(2,k) - K*x(1,k))/m,...
%         x(4,k),...
%         (u1(1,k)-R*x(4,k))/L]';
%     %Ax+Bu+Ed
%     dx_nec = dx_nec + E_dec*d_dec(:,k);
%     
    
%     x(:,k+1) = x(:,k)+Te*dx_nec;
   x(:,k+1) = Ad2*x(:,k) + B*Te*u1(:,k)+ E_dec*d_dec(:,k);
    y(:,k+1) = C*x(:,k+1);
    d_dec(:,k+1) = rand()*110;
%     d_dec(:,k+1) = 0.1;
   
    
    z(:,k+1) = F*z(:,k)+T*B*Te*u1(:,k)+(K_dec+K_dec2)*y(:,k);
    x_es(:,k+1) = z(:,k+1)+H*y(:,k+1);
    
% %     if(k< K_timp*0.4)
%     u1(:,k+1) =  [1e-3*sin(20*Te*k);1e-3*cos(38*Te*k)];
% %     else
% %         u1(:,k+1) = [0,0];
% %     end

     u1(:,k+1) = u1(:,k);
end
figure
plot(u1')
legend('u1','u2')
figure
plot(x')
legend('x1','x2','x3','x4');
figure
plot(x_es')
title('Estimat');
legend('x1','x2','x3','x4');
figure
plot(x' - x_es')
title('Eroare');
legend('x1','x2','x3','x4');