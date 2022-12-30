% Drone

clear all, close all, clc

set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

sstdarkblue     = [0,73,219]/255;
sstblue         = [0,128,255]/255;
sstlightblue    = [48,208,216]/255;
sstdarkgreen    = [20,120,60]/255;
sstgreen        = [43,191,92]/255;
sstlightgreen   = [140,255,200]/255;
sstdarkgray     = [150,150,150]/255;
sstgray         = [190,190,190]/255;
sstlightgray    = [230,230,230]/255;

%% 1.1.8 2-D Model

clear X Y
    
g = 9.8065; % gravity force [m/s^2]
beta = 1/2*1.225*0.3*0.1; % 1/2*p*Cd*A where p is air density, Cd is drag coefficient, A is frontal Area
v1e = 3;
v1e_x = 1; % speed of the follower drone
v1e_y = sqrt(v1e^2-v1e_x^2);
v0e_x = v1e_x ; % speed of the leader drone
v0e_y = v1e_y;
p10 = 0.5;
p10_x = 0.1;
p10_y = sqrt(p10^2-p10_x^2);
m = 0.1; % mass in kg
we_x = 0; % typ. wind speed [m/s]
we_y = 0;
theta1e = atan((beta/(m*g))*(v0e_x-we_x)^2); % pitch angle in rad
theta1_e = rad2deg(theta1e),
phi1e = atan(-(beta/(m*g))*cos(theta1e)*(v0e_y-we_y)^2);% roll angle in rad
phi1_e = rad2deg(phi1e)
Al = [0 0 -1 0; 0 0 0 -1; 0 0 -2*beta*v1e_x/m 0;0 0 0 -2*beta*v1e_y/m];
Bl = [0 0; 0 0;(g*sec(theta1e)^2) 0;(-g*sec(theta1e)*tan(theta1e)*tan(phi1e)) (-g*sec(phi1e)^2*sec(theta1e))];
El = [0 0 0 1;0 0 0 1; 2*beta*we_x/m 0 0 0; 0 2*beta*we_y/m 0 0];
Cl = eye(4); % Cl = [1 0 0 0; 0 1 0 0] but we wanted to access speed
Dl = 0;

% Discrete-time state-space model approximation

Ts = 0.1; % sampling time [s]

Ad = [1 0 -Ts 0; 0 1 0 -Ts; 0 0 1-2*beta*v1e_x*Ts/m 0;0 0 0 1-2*beta*v1e_y*Ts/m];
Bd = [0 0; 0 0;Ts*(g*sec(theta1e)^2) 0;Ts*(-g*sec(theta1e)*tan(theta1e)*tan(phi1e)) Ts*(-g*sec(phi1e)^2*sec(theta1e))];
Ed = [0 0 0 Ts;0 0 0 Ts; 2*beta*we_x*Ts/m 0 0 0; 0 2*beta*we_y*Ts/m 0 0];
Cd = eye(4);  % Cd = [1 0 0 0; 0 1 0 0] but we wanted to access speed
Dd = 0;

%  Stability, controllability, and observability of the discrete time model.

[eigenvectors,eigenvalues] = eig(Ad); % stability -> eigenvalues <= 1
[similarity_transform,J] = jordan(Ad);
disp('System marginally stable because eigenvectors <= 1');
disp(eigenvalues);
disp('and all jordan blocks associated with eigenvalues = 1 are 1x1 blocks');
disp(J);

Co = ctrb(Ad,Bd); % controlable - rank of controlability matrix equals n = 2
if rank(Co) == size(Ad,1)
   disp('Controlable') 
else
    disp('Not Controlable') 
end    

Ob = obsv(Ad,Cd); % observable - rank of observability matrix equals n = 2
if rank(Ob) == size(Ad,1)
   disp('Observable') 
else
    disp('Not Observable') 
end   

% Simulate system response

% Continuous response
sys_c = ss(Al,Bl,Cl,Dl); 

tmax = 5; % duration of the simulation
t = 0:Ts:tmax;
u = 0.1*ones(length(t),2);
x0 = [0 0 0 0]'; % [p10x p10y v1x v1y]'
y_c = lsim(sys_c,u,t,x0);
y_c(:,1) = y_c(:,1) + p10_x;
y_c(:,2) = y_c(:,2) + p10_y;
y_c(:,3) = y_c(:,3) + v1e_x;
y_c(:,4) = y_c(:,4) + v1e_y;

% Discrete response
Ts = 0.1;

sys_d = ss(Ad,Bd,Cd,Dd,Ts);
y_d = lsim(sys_d,u,t,x0);

figure('Name','Pole Zero Map - 2-D Linearized and Discretized System','NumberTitle','off');
set(gcf,'defaultTextInterpreter','tex');
pzplot(sys_d,'b');
grid on;
axis equal;
% change marker size
a = findobj(gca,'type','line');
for i = 1:length(a)
    set(a(i),'markersize',12) %change marker size
    set(a(i), 'linewidth',2)  %change linewidth
    set(a(i),'color',sstblue) %change marker size
end
title('Pole Zero Map');

y_d(:,1) = y_d(:,1) + p10_x;
y_d(:,2) = y_d(:,2) + p10_y;
y_d(:,3) = y_d(:,3) + v1e_x;
y_d(:,4) = y_d(:,4) + v1e_y;

% sys_d_z = c2d(sys_c, Ts, 'zoh');
% y_d_z = lsim(sys_d_z,u,t,x0);
% y_d_z(:,1) = y_d_z(:,1) + p10_x;
% y_d_z(:,2) = y_d_z(:,2) + p10_y;
% y_d_z(:,3) = y_d_z(:,3) + v1e_x;
% y_d_z(:,4) = y_d_z(:,4) + v1e_y;
% 
% sys_d_t = c2d(sys_c, Ts, 'tustin');
% y_d_t = lsim(sys_d_t,u,t,x0);
% y_d_t(:,1) = y_d_t(:,1) + p10_x;
% y_d_t(:,2) = y_d_t(:,2) + p10_y;
% y_d_t(:,3) = y_d_t(:,3) + v1e_x;
% y_d_t(:,4) = y_d_t(:,4) + v1e_y;
Ny = length(y_d(:,1));
y_c_total = zeros(Ny,2);
y_d_total = zeros(Ny,2);
for i = 1:size(y_d(:,1))
    y_c_total(i,1) = sqrt(y_c(i,1)^2 + y_c(i,2)^2);
    y_c_total(i,2) = sqrt(y_c(i,3)^2 + y_c(i,4)^2);
    y_d_total(i,1) = sqrt(y_d(i,1)^2 + y_d(i,2)^2);
    y_d_total(i,2) = sqrt(y_d(i,3)^2 + y_d(i,4)^2);
end

figure('Name','Phase Plot - 2-D Linearized and Discretized System','NumberTitle','off');
subplot(2,1,1)
plot(y_d(:,1),y_d(:,3),'s-','Color',sstblue);
grid on;
hold off;
xlabel('$$p_{10_x}$$');
ylabel('$$v_{1_x}$$');
legend('trajectory');
title('Phase plot');

subplot(2,1,2)
plot(y_d(:,2),y_d(:,4),'s-','Color',sstblue);
grid on;
hold off;
xlabel('$$p_{10_y}$$');
ylabel('$$v_{1_y}$$');
legend('trajectory');
title('Phase plot');

figure('Name','2-D Simulation: 4 states','NumberTitle','off');
subplot(2,1,1);
plot(t,y_c(:,1),'Color',sstgray,'LineStyle','-');
hold on
stairs(t,y_d(:,1),'Color',sstlightgray,'LineStyle','-');
hold on
plot(t,y_c(:,2),'Color',sstgreen,'LineStyle','-');
hold on
stairs(t,y_d(:,2),'Color',sstblue,'LineStyle','-');
hold on
stairs(t,y_c_total(:,1),'Color','r','LineStyle','--');
hold on
stairs(t,y_d_total(:,1),'Color','r','LineStyle','-');
% stairs(t,y_d_z(:,1),'Color',sstgray,'LineStyle','-');
% hold on
% stairs(t,y_d_t(:,1),'Color',sstdarkblue,'LineStyle','-');
grid on;
xlabel('t [s]');
ylabel('$p_{10}$');
title('Position');
legend('$p_{10_x}$ cont.','$p_{10_x}$ disc.','$p_{10_y}$ cont.','$p_{10_y}$ disc.','$p_{10}$ cont.','$p_{10}$ disc.');

subplot(2,1,2);
plot(t,y_c(:,3),'Color',sstgray,'LineStyle','-');
hold on
stairs(t,y_d(:,3),'Color',sstlightgray,'LineStyle','--');
hold on
plot(t,y_c(:,4),'Color',sstgreen,'LineStyle','-');
hold on
stairs(t,y_d(:,4),'Color',sstblue,'LineStyle','--');
hold on
stairs(t,y_c_total(:,2),'Color','r','LineStyle','--');
hold on
stairs(t,y_d_total(:,2),'Color','r','LineStyle','-');
grid on;
xlabel('t [s]');
ylabel('$v_{1}$');
title('Velocity');
legend('$v_{1_x}$ cont.','$v_{1_x}$ disc.','$v_{1_y}$ cont.','$v_{1_y}$ disc.','$v_1$ cont.','$v_1$ disc.');
hold off;

%% 1.1 Linear 1-D Model 

clear X Y
    
g = 9.8065; % gravity force [m/s^2]
beta = 1/2*1.225*0.3*0.1; % 1/2*p*Cd*A where p is air density, Cd is drag coefficient, A is frontal Area
v1e = 3; % speed of the follower drone
v0e = 3; % speed of the leader drone
m = 0.1; % mass in kg
we = 0; % typ. wind speed [m/s]
theta1e = atan((beta/(m*g))*(v0e-we)^2); % pitch angle in rad


Al = [0 -1; 0 -2*beta*v1e/m];
Bl = [0 g*sec(theta1e)^2]';
El = [0 1; 2*beta*we/m 0];
Cl = eye(2); % Cl = [1 0] but we wanted to access speed
Dl = 0;

% Discrete-time state-space model approximation

Ts = 0.1; % sampling time [s]

Ad = [1 -Ts; 0 1-2*beta*v1e*Ts/m];
Bd = [0 Ts*g*sec(theta1e)^2]';
Ed = [0 Ts; 2*beta*we*Ts/m 0];
Cd = eye(2);  % Cd = [1 0] but we wanted to access speed
Dd = 0;

%  Stability, controllability, and observability of the discrete time model.

[eigenvectors,eigenvalues] = eig(Ad); % stability -> eigenvalues <= 1
[similarity_transform,J] = jordan(Ad);
disp('System marginally stable because eigenvectors <= 1');
disp(eigenvalues);
disp('and all jordan blocks associated with eigenvalues = 1 are 1x1 blocks');
disp(J);

Co = ctrb(Ad,Bd); % controlable - rank of controlability matrix equals n = 2
if rank(Co) == size(Ad,1)
   disp('Controlable') 
else
    disp('Not Controlable') 
end    

Ob = obsv(Ad,Cd); % observable - rank of observability matrix equals n = 2
if rank(Ob) == size(Ad,1)
   disp('Observable') 
else
    disp('Not Observable') 
end   

% Simulate system response

% Continuous response
sys_c = ss(Al,Bl,Cl,Dl); 

tmax = 5; % duration of the simulation
t = 0:Ts:tmax;
u = 0.1*ones(length(t),1);
x0 = [0 0]';
y_c = lsim(sys_c,u,t,x0);
y_c(:,1) = y_c(:,1) + 0.5;
y_c(:,2) = y_c(:,2) + v1e;

% Discrete response
Ts = 0.1;

sys_d = ss(Ad,Bd,Cd,Dd,Ts);
y_d = lsim(sys_d,u,t,x0);

figure('Name','Pole Zero Map - Linearized and Discretized System','NumberTitle','off');
set(gcf,'defaultTextInterpreter','tex');
pzplot(sys_d,'b');
grid on;
axis equal;
% change marker size
a = findobj(gca,'type','line');
for i = 1:length(a)
    set(a(i),'markersize',12) %change marker size
    set(a(i), 'linewidth',2)  %change linewidth
    set(a(i),'color',sstblue) %change marker size
end
title('Pole Zero Map');

y_d(:,1) = y_d(:,1) + 0.5;
y_d(:,2) = y_d(:,2) + v1e;

sys_d_z = c2d(sys_c, Ts, 'zoh');
y_d_z = lsim(sys_d_z,u,t,x0);
y_d_z(:,1) = y_d_z(:,1) + 0.5;
y_d_z(:,2) = y_d_z(:,2) + v1e;

sys_d_t = c2d(sys_c, Ts, 'tustin');
y_d_t = lsim(sys_d_t,u,t,x0);
y_d_t(:,1) = y_d_t(:,1) + 0.5;
y_d_t(:,2) = y_d_t(:,2) + v1e;

figure('Name','Phase Plot - Linearized and Discretized System','NumberTitle','off');
plot(y_d(:,1),y_d(:,2),'s-','Color',sstblue);
grid on;
hold off;
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend('trajectory');
title('Phase plot');

figure('Name','Simulation of both states for u = 0.1 rad','NumberTitle','off');
subplot(2,1,1);
plot(t,y_c(:,1),'Color',sstgreen,'LineStyle','-');
hold on
stairs(t,y_d(:,1),'Color',sstblue,'LineStyle','-');
hold on
stairs(t,y_d_z(:,1),'Color',sstgray,'LineStyle','-');
hold on
stairs(t,y_d_t(:,1),'Color',sstdarkblue,'LineStyle','-');
grid on;
xlabel('t [s]');
ylabel('$p_{10}$');
title('Position $p_{10}$');
legend('cont.','disc. euler','disc. zoh', 'dic. tustin');

subplot(2,1,2);
plot(t,y_c(:,2),'Color',sstgreen,'LineStyle','-');
hold on
stairs(t,y_d(:,2),'Color',sstblue,'LineStyle','-');
grid on;
xlabel('t [s]');
ylabel('$v_{1}$');
title('Velocity $v_{1}$');
legend('cont.','disc.');
hold off;

% Same but with possibility to simulate some disturbances
nu = size(Bd,2);
nk = tmax/Ts;
TX = 1:nk+1;
x0 = [0 0]';
U = 0.1*ones(nu,nk);
Ud = [0 1]'; % wind and v0 disturbances 
X(:,1) = x0;
X(:,2) = x0;
for k = 2:nk
    X(:,k+1) = Ad*X(:,k) + Bd*U(:,k) + Ed*Ud;
    Y(:,k+1) = Cd*X(:,k+1);
end
Y(1,:) = Y(1,:) + 0.5;
figure('Name','Simulation of both states with disturbances','NumberTitle','off');
subplot(2,1,1);
stairs(t,y_d(:,1),'Color',sstgreen,'LineStyle','-');
hold on
stairs(t,Y(1,:),'Color',sstblue,'LineStyle','-');
grid on;
xlabel('t [s]');
ylabel('$p_{10}$');
title('Position $p_{10}$');
legend('disc.','disc.with E matrix');

v1d = v1e + X(2,:);
subplot(2,1,2);
stairs(t,y_d(:,2),'Color',sstblue,'LineStyle','-');
hold on;
stairs(t,v1d,'-','Color',sstgreen);
grid on;
hold off;
xlabel('t [s]');
ylabel('$$v_{1}$$');
legend('without perturb.','with perturb.');
title('$$v_{1}$$');



%% 1.2 Unconstrained MPC

set(0,'defaultTextInterpreter','latex');

Cl = [1 0];
Cd = [1 0];
clear X Xd Y U dUopt Uopt;

N = 5;
P = 10*eye(1);
Q = 10*eye(1);
R = .1;
xd0 = [0 0]';
nx = size(Bd,1);
nu = size(Bd,2);

% compute matrices
[F,G,Qb,Rb,H,Fd,Gd,Hd,A,B,C] = GetBatchXiMatrices(Ad,Bd,Cd,N,P,Q,R);
Fb = H*F;
Gb = H*G;
Fdb = Hd*Fd;
Gdb = Hd*Gd;

% compute final cost function matrices and control gains
Rt = Gb'*Qb*Gb + Rb,
St = Gb'*Qb,
Ky = Rt^(-1)*St,
K  = Rt^(-1)*St*Fb,

%simulate controlled system:
nk = 200;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
ref = square_dt(nk+N,80); % square reference of amplitude -1 1 on the real system but -1.5 0.5 on simulated linearized system because p10e = 0.5
dist_x1 = 0*0.5*ones(size(ref)).*(Tref>=0);
x0 = [xd0*0 ; Cd*xd0];
U = zeros(nu,nk);
dU = zeros(nu,nk);
Xd(:,1) = xd0;
X(:,1) = x0;
Y(:,1) = Cd*xd0;
Xd(:,2) = xd0;
X(:,2) = x0;
Y(:,2) = Cd*xd0;
for k = 2:nk
    
    % compute initial condition and current reference sequence
    Yb = ref(:,k:k+N)';
    Dxdk = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdk; Cd*Xd(:,k)];
    xk = X(:,k);
    
    % compute optimal incremental control sequence:
    dUopt(:,:,k) = reshape(-(K*xk-Ky*Yb) ,nu,N);

    % set MPC control policy and simulate system:
    Uopt(:,:,k) = U(:,k-1) + dUopt(:,:,k);
    uk = Uopt(:,1,k);
    U(:,k) = uk;
    
    % simulate original system:
    Xd(:,k+1) = Ad*Xd(:,k) + Bd*U(:,k) + [dist_x1(:,k);0];
    Y(:,k+1) = Cd*Xd(:,k+1);
        
end
X(:,k+1) = [ Xd(:,k+1)-Xd(:,k) ; Cd*Xd(:,k+1)];
    
figure('Name','Phase Plot - Unconstrained MPC','NumberTitle','off');
plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
grid on;
hold off;
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend('trajectory');
title('Phase plot');

figure('Name','State Evolution - Unconstrained MPC','NumberTitle','off');
plot(Tref(1:end-5),ref(1:end-5),'gs-');
grid on;
hold on;
plot(TX(1:end-1),Xd(1,1:end-1),'s-','Color',sstblue);
plot(TX(1:end-1),Xd(2,1:end-1),'d-','Color',sstgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('ref','$$x_1$$','$$x_2$$','Location','SouthEast');
title('State evolution');

figure('Name','Control Action - Unconstrained MPC','NumberTitle','off');
plot(1:length(U(1,:,:)),U(1,:,:),'s-','Color',sstblue);
grid on;
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('input');
title('Input');


figure('Name','$$\Delta$$ State evolution - Unconstrained MPC','NumberTitle','off');
plot(TX(1:end-1),X(1,1:end-1),'s-','Color',sstblue);
grid on;
hold on;
plot(TX(1:end-1),X(2,1:end-1),'d-','Color',sstgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$\Delta x(t_k)$$');
legend('$$\Delta x_1$$','$$\Delta x_2$$','Location','SouthEast');
title('$$\Delta$$ State evolution');


figure('Name','Simulation - Unconstrained MPC','NumberTitle','off');
subplot(2,1,1)
plot(TX(1:end-1),Xd(1,1:end-1),'s-','Color',sstblue);
grid on;
hold on;
plot(Tref(1:end-5),ref(1,1:end-5),'s-','Color',sstgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$y(t_k)$$');
legend('output','reference');
title('Output');
subplot(2,1,2)
plot(1:length(U(1,:,:)),U(1,:,:),'s-','Color',sstblue);
grid on;
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('input');
title('Input');

%% 1.3 Constrained MPC

clear X Xd Xd2 Y U U2 dUopt dUopt2 Uopt Uopt2;
Cl = [1 0];
Cd = [1 0];

N = 5; %5
P = 10*eye(1); %0.8
Q = 10*eye(1); %0.8
R = 0.1; %0.01
xd0 = [0 0]';
nx = size(Bd,1);
nu = size(Bd,2);

% compute batch matrices
[F,G,Qb,Rb,H,Fd,Gd,Hd] = GetBatchXiMatrices(Ad,Bd,Cd,N,P,Q,R);
Fb = H*F;
Gb = H*G;
Fdb = Hd*Fd;
Gdb = Hd*Gd;
Rt = Gb'*Qb*Gb + Rb;
St = Gb'*Qb;
Ky = Rt^(-1)*St,
K = Rt^(-1)*St*Fb,

% compute constraints matrices:
u_max = deg2rad(30);
y_max = 5-0.5; %p10max - p10e (to adjust for linearized system)
y_min = 0.2-0.5;
U_max = kron(u_max,ones(N,1));
Y_max = kron(y_max,ones(N+1,1));
Y_min = kron(y_min,ones(N+1,1));
M3 = tril(ones(N*nu));
M4 = ones(N*nu,nu);
Mu = [-M3;M3];
My = [-Gb;Gb];
M = [Mu;My];

%simulate controlled system:
nk = 200;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
ref = square_dt(nk+N,80);% square reference of amplitude -1 1, on the real system this would be -.5 to 1.5 because p10e = 0.5
dist_x1 = 0*.5*ones(size(ref)).*(Tref>=0);
x0 = [xd0*0 ; Cd*xd0];
U = zeros(nu,nk);
U2 = zeros(nu,nk);
Xd(:,1) = xd0;
Xd2(:,1) = xd0;
X(:,1) = x0;
Xd(:,2) = xd0;
Xd2(:,2) = xd0;
X(:,2) = x0;
for k = 2:nk
    
    % compute initial conditions
    Yb = ref(:,k:k+N)';    
    Dxdk = Xd(:,k)-Xd(:,k-1);
    X(:,k) = [ Dxdk; Cd*Xd(:,k)];
    xk = X(:,k);
    u_1 = U(:,k-1);
    wu = [U_max + M4*u_1;U_max - M4*u_1];
    wy = [-Y_min + Fb*xk;Y_max - Fb*xk];
    w = [wu;wy];
    Dxdk2 = Xd2(:,k)-Xd2(:,k-1);
    xk2 = [ Dxdk2; Cd*Xd2(:,k)];
    
    % compute constrained optimal incremental control sequence and MPC policy
    [dUo,Jo,exitflag,output,lambda] = quadprog(2*Rt,2*St*(Fb*xk-Yb),M,w);
    if exitflag~=1
        error('Problems in the Optimization problem.');
    end
    dUopt(:,:,k) = reshape( dUo ,nu,N);
    Uopt(:,:,k) = U(:,k-1) + dUopt(:,:,k);
    U(:,k) = Uopt(:,1,k);
    
    % simulate original system:
    Xd(:,k+1) = Ad*Xd(:,k) + Bd*U(:,k) + Ed*[dist_x1(:,k);0];
    
    % compute unconstrained optimal sequence and MPC policy
    dUopt2(:,:,k) = reshape(-(K*xk2-Ky*Yb) ,nu,N);
    Uopt2(:,:,k) = U2(:,k-1) + dUopt2(:,:,k);
    U2(:,k) = Uopt2(:,1,k);
    
    % simulate system 2 (for comparison):
    Xd2(:,k+1) = Ad*Xd2(:,k) + Bd*U2(:,k);
    
end
X(:,k+1) = [ Xd(:,k+1)-Xd(:,k) ; Cd*Xd(:,k+1)];
    
figure('Name','Phase Plot - Constrained MPC','NumberTitle','off');
plot(Xd(1,:),Xd(2,:),'s-','Color',sstblue);
grid on;
hold on;
plot(Xd2(1,:),Xd2(2,:),'o-','Color',sstgray);
hold off;
xlabel('$$x_1$$');
ylabel('$$x_2$$');
legend('const.','unconst.');
title('Phase plot');

figure('Name','State Evolution - Constrained MPC','NumberTitle','off');
plot(Tref(1,1:end-5),ref(1,1:end-5),'-');
grid on;
hold on;
plot(TX(1:end-1),Xd(1,1:end-1),'s-','Color',sstblue);
plot(TX(1:end-1),Xd(2,1:end-1),'s-','Color',sstgreen);
plot(TX(1:end-1),Xd2(1,1:end-1),'-','Color',sstgray);
plot(TX(1:end-1),Xd2(2,1:end-1),'-','Color',sstlightgray);
hold off;
xlabel('$$t_k$$');
ylabel('$$x(t_k)$$');
legend('ref','$$x_1$$ const.','$$x_2$$ const.','$$x_1$$ unc.','$$x_2$$ unc.','Location','SouthEast');
title('State evolution');

figure('Name','Simulation - Constrained MPC','NumberTitle','off');
subplot(2,1,1)
plot(TU,U,'-','Color',sstblue);
grid on;
hold on;
plot(TU,U2,'-','Color',sstgray);
hold off;
xlabel('$$t_k$$');
ylabel('$$u(t_k)$$');
legend('const.','unc.');
title('Input');

subplot(2,1,2)
plot(TX(1:end-1),Xd(1,1:end-1),'-','Color',sstblue);
grid on;
hold on;
plot(Tref(1,1:end-5),ref(1,1:end-5),'-','Color',sstgreen);
hold off;
xlabel('$$t_k$$');
ylabel('$$y(t_k)$$');
legend('output','reference');
title('Output');

%% 2.1 Code: Unconstrained 3 player centralized, decentralized, distributed


np = 3;
% np = 10;


A22 = -2*beta*v1e/m;%m1 = m
A44 = -2*beta*v1e/m;%v2e = v1e && m2 = m
A66 = -2*beta*v1e/m;%v3e = v1e && m3 = m

Al = [0 -1  0  0  0  0 ;
      0 A22 0  0  0  0 ;
      0  1  0 -1  0  0 ;
      0  0  0 A44 0  0 ;
      0  0  0  1  0 -1 ;
      0  0  0  0  0 A66];

B21 = g*sec(theta1e)^2;
B42 = g*sec(theta1e)^2;%theta2e = theta1e
B63 = g*sec(theta1e)^2;%theta3e = theta1e


Bl = [ 0  0  0 ;
      B21 0  0 ;
       0  0  0 ;
       0 B42 0 ;
       0  0  0 ;
       0  0 B63];

Cl = [1 0 0 0 0 0 ;
      0 0 1 0 0 0 ;
      0 0 0 0 1 0 ];

E21 = 2*beta*we/m;%m1 = m
E41 = 2*beta*we/m;%m2 = m
E61 = 2*beta*we/m;%m3 = m

El = [ 0  1 ;
      E21 0 ;
       0  0 ;
      E41 0 ;
       0  0 ;
      E61 0 ];

Ad = eye(size(Al,1)) + Ts*Al;
Bd = Ts * Bl;
Ed = Ts * El;
Cd = Cl;

Ai = {Ad(1:2,1:2),Ad(2:4,2:4),Ad(4:6,4:6)};
Aidec = {Ad(1:2,1:2),Ad(3:4,3:4),Ad(5:6,5:6)};

Bij = { Bd(1:2,1), Bd(1:2,2), Bd(1:2,3)
        Bd(2:4,1), Bd(2:4,2), Bd(2:4,3)
        Bd(4:6,1), Bd(4:6,2), Bd(4:6,3) };
Bijdec = { Bd(1:2,1), Bd(1:2,2), Bd(1:2,3)
           Bd(3:4,1), Bd(3:4,2), Bd(3:4,3)
           Bd(5:6,1), Bd(5:6,2), Bd(5:6,3) };
Ci = {Cd(1,1:2),Cd(2,2:4),Cd(3,4:6)};
Cidec = {Cd(1,1:2),Cd(2,3:4),Cd(3,5:6)};
xd10 = [0 0]';
xd20 = xd10;
x0 = [xd10; xd10; xd10];

nu = size(Bd,2);
ny = size(Cd,1);
nx = size(Bd,1);
nxi = size(Bij{1,1},1);
nui = size(Bij{1,1},2);
nud = size(Ed,2);
% cost parameters
N = 3;
Pi = 10;
Qi = 10;
Ri = .1;
alphai = 1;
% distributed steps parameters
w1 = 0.4;
w2 = 0.3;
w3 = 1-w1-w2;

% compute centralized tracking controller
P = blkdiag(Pi,Pi,Pi);
Q = blkdiag(Qi,Qi,Qi);
R = blkdiag(Ri,Ri,Ri);


[Fbc,Gbc,Qbc,Rbc] = GetYMats(Ad,Bd,Cd,N,P,Q,R);
[Rtc,Sxc,Syc] = GetYMPC(Fbc,Gbc,Qbc,Rbc);

% compute decentralized tracking controllers
[Fbd1,Gbd1,Qbd1,Rbd1] = GetYMats(Aidec{1},Bijdec{1,1},Cidec{1},N,Pi,Qi,Ri);
[Rtd1,Sxd1,Syd1] = GetYMPC(Fbd1,Gbd1,Qbd1,Rbd1);
[Fbd2,Gbd2,Qbd2,Rbd2] = GetYMats(Aidec{2},Bijdec{2,2},Cidec{2},N,Pi,Qi,Ri);
[Rtd2,Sxd2,Syd2] = GetYMPC(Fbd2,Gbd2,Qbd2,Rbd2);
[Fbd3,Gbd3,Qbd3,Rbd3] = GetYMats(Aidec{3},Bijdec{3,3},Cidec{3},N,Pi,Qi,Ri);
[Rtd3,Sxd3,Syd3] = GetYMPC(Fbd3,Gbd3,Qbd3,Rbd3);

% compute distributed tracking controllers
[Fb,Gb,Qb,Rb] = GetYMats(Ai,Bij,Ci,N,Pi,Qi,Ri,alphai);
[Rt,Sx,Sy,Su] = GetYMPC(Fb,Gb,Qb,Rb);

% compute constraints matrices:
%Ui_max = 0.5*ones(N*nui,1);
Mui = [];
wui = [];
%U_max = 0.5*ones(N*nu,1);
Mu = [];
wu = [];

% analyze stability of distributed law:
% L = [ (1-w1)*eye(size(L1))  , w1*L1
%       w2*L2                 , (1-w2)*eye(size(L2)) ];
% maxL = max(abs(eig(L))),
% Kb = [K11,K12;K21,K22];
% Kby = [K11y,K12y;K21y,K22y];
% K = E*(1-L)^(-1)*Kb;
% lbd_cl = eig(A+B*K),

%simulate controlled system:
nk = 200;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
% ref = [1*(Tref>=0)+1.5*(Tref>=40)-1*(Tref>=80);+0.5*(Tref>=0)+1.5*(Tref>=60)+0.5*(Tref>=90);1*(Tref>=0)+2*(Tref>=30)-0.5*(Tref>=70)]; % step reference of amplitude 1, starting at k=30
ref = [1*(Tref>=0)-0.5*(Tref>=80);1*(Tref>=0)-0.5*(Tref>=80);1*(Tref>=0)-0.5*(Tref>=80)];
uc = zeros(nu,nk);
ud = zeros(nu,nk);
u = zeros(nu,nk);
udd = [0*ones(nud-1,nk);0*ones(nud-1,nk)];
U1 = zeros(nui,N,nk);
U2 = zeros(nui,N,nk);
U3 = zeros(nui,N,nk);
xc = zeros(nx,nk+1);
xd = zeros(nx,nk+1);
x = zeros(nx,nk+1);
xc(:,1) = x0;
xd(:,1) = x0;
x(:,1) = x0;
xc(:,2) = x0;
xd(:,2) = x0;
x(:,2) = x0;
for k = 2:nk
    Yb = reshape(ref(:,k:k+N),[],1);
    Yb1 = ref(1,k:k+N)';
    Yb2 = ref(2,k:k+N)';
    Yb3 = ref(3,k:k+N)';
    
    % centralized
    Stc = Sxc*xc(:,k)-Syc*Yb;
    [Uco,Jco,exitflag] = quadprog(Rtc,Stc,Mu,wu);
    if exitflag~=1, error('Problems in centralized.'); end
    Ucopt(:,:,k) = reshape( Uco ,nu,N);
    uc(:,k) = Ucopt(:,1,k);

    % simulate system for centralized MPC
    xc(:,k+1) = Ad*xc(:,k) + Bd*uc(:,k) + Ed*udd(:,k); %simulate joint system

    % decentralized
    Std1 = Sxd1*xd(1:2,k)-Syd1*Yb1;
    Std2 = Sxd2*xd(3:4,k)-Syd2*Yb1;
    Std3 = Sxd3*xd(5:6,k)-Syd3*Yb1;
    [Ud1o,Jd1o,ed1] = quadprog(Rtd1,Std1,Mui,wui);
    [Ud2o,Jd2o,ed2] = quadprog(Rtd2,Std2,Mui,wui);
    [Ud3o,Jd3o,ed3] = quadprog(Rtd3,Std3,Mui,wui);
    if any([ed1,ed2,ed3]~=1), error('Problems in decentralized MPC.'); end
    Ud1opt = reshape( Ud1o ,nui,N);
    Ud2opt = reshape( Ud2o ,nui,N);
    Ud3opt = reshape( Ud3o ,nui,N);
    ud(:,k) = [ Ud1opt(:,1) ; Ud2opt(:,1) ; Ud3opt(:,1) ];

    % simulate system for decentralized MPC
    xd(:,k+1) = Ad*xd(:,k) + Bd*ud(:,k) + Ed*udd(:,k); %simulate joint system
    
   
    % distributed MPC
    x1k = x(1:2,k);
    x2k = x(2:4,k);
    x3k = x(4:6,k);
    U1p = reshape( U1(:,:,k-1) ,[],1);
    U2p = reshape( U2(:,:,k-1) ,[],1);
    U3p = reshape( U3(:,:,k-1) ,[],1);

    % Get optimal sequence for player 1
    St1 = Sx{1,1}*x1k - Sy{1,1}*Yb1;
    [U1o,J1o,exitflag] = quadprog(Rt{1},St1,Mui,wui);
    if exitflag~=1, error('Problems in player 1.'); end
    % Get optimal sequence for player 2
    St2 = Sx{2,1}*x1k - Sy{2,1}*Yb1 + Sx{2,2}*x2k - Sy{2,2}*Yb1 + Su{2,1}*U1p;
    [U2o,J2o,exitflag] = quadprog(Rt{2},St2,Mui,wui);
    if exitflag~=1, error('Problems in player 2.'); end
    % Get optimal sequence for player 3
    St3 = Sx{3,2}*x2k - Sy{3,2}*Yb1 + Sx{3,3}*x3k - Sy{3,3}*Yb1 + Su{3,2}*U2p;
    [U3o,J3o,exitflag] = quadprog(Rt{3},St3,Mui,wui);
    if exitflag~=1, error('Problems in player 3.'); end
    % convex step iterates
    U1pp = U1p;
    U2pp = U2p;
    U3pp = U3p;
    for p = 1:np
        U1pp = w1*U1o + (1-w1)*U1pp;
        U2pp = w2*U2o + (1-w2)*U2pp;
        U3pp = w3*U3o + (1-w3)*U3pp;
    end
    U1(:,:,k) = reshape( U1pp ,nui,N);
    U2(:,:,k) = reshape( U2pp ,nui,N);
    U3(:,:,k) = reshape( U3pp ,nui,N);
    
    u(:,k) = [ U1(:,1,k) ; U2(:,1,k) ; U3(:,1,k) ]; % apply first input to each player

    % simulate system for distributed MPC
    x(:,k+1) = Ad*x(:,k) + Bd*u(:,k) + Ed*udd(:,k); %simulate joint system
end

figure('Name','Simulation of both states for u = 0.1 rad for three drones (1)','NumberTitle','off');
subplot(2,1,1);
plot(Tref,ref(1,:),'s-','Color',sstlightgray);
hold on
plot(Tref,ref(2,:),'s-','Color',sstlightgray);
plot(Tref,ref(3,:),'s-','Color',sstlightgray);
plot(TX,xc(1,:),'.--','Color',sstblue);
plot(TX,xc(3,:),'Color',sstblue);
plot(TX,xc(5,:),'s-','Color',sstblue);
plot(TX,xd(1,:),'.--','Color','r');
plot(TX,xd(3,:),'Color','r');
plot(TX,xd(5,:),'s-','Color','r');
xlabel('$$t_k$$')
ylabel('$p_{ij}$ [m]');
grid on
title('State $$p_{ij}$$ evolution');
legend('$$r_1$$','$$r_2$$','$$r_3$$','$$p10$$ cen.','$$p21$$ cen.','$$p32$$ cen.','$$p10$$ decen.','$$p21$$ decen.','$$p32$$ decen.')
hold off

subplot(2,1,2);
plot(Tref,0*ref(1,:),'s-','Color',sstlightgray);
hold on
plot(Tref,0*ref(2,:),'s-','Color',sstlightgray);
plot(Tref,0*ref(3,:),'s-','Color',sstlightgray);
plot(TX,xc(2,:),'.--','Color',sstblue);
plot(TX,xc(4,:),'Color',sstblue);
plot(TX,xc(6,:),'s-','Color',sstblue);
plot(TX,xd(2,:),'.--','Color','r');
plot(TX,xd(4,:),'Color','r');
plot(TX,xd(6,:),'s-','Color','r');
xlabel('$$t_k$$')
ylabel('$v_{i}$ [m/s]');
grid on
title('State $$v_{i}$$ evolution');
legend('$$r_1$$','$$r_2$$','$$r_3$$','$$v_1$$ cen.','$$v_2$$ cen.','$$v_3$$ cen.','$$v_1$$ decen.','$$v_2$$ decen.','$$v_3$$ decen.')
hold off

figure('Name','Simulation of both states for u = 0.1 rad for three drones (2)','NumberTitle','off');
subplot(2,1,1);
plot(Tref,ref(1,:),'s-','Color',sstlightgray);
hold on
plot(Tref,ref(2,:),'s-','Color',sstlightgray);
plot(Tref,ref(3,:),'s-','Color',sstlightgray);
plot(TX,xc(1,:),'.--','Color',sstblue);
plot(TX,xc(3,:),'Color',sstblue);
plot(TX,xc(5,:),'s-','Color',sstblue);
plot(TX,xd(1,:),'.--','Color','r');
plot(TX,xd(3,:),'Color','r');
plot(TX,xd(5,:),'s-','Color','r');
plot(TX,x(1,:),'.--','Color','g');
plot(TX,x(3,:),'Color','g');
plot(TX,x(5,:),'s-','Color','g');
title('State $$p_{ij}$$ evolution');
xlabel('$$t_k$$');
ylabel('$p_{ij}$ [m]');
grid on
legend('$$r_1$$','$$r_2$$','$$r_3$$','$$p10$$ cen.','$$p21$$ cen.','$$p32$$ cen.','$$p10$$ decen.','$$p21$$ decen.','$$p32$$ decen.','$$p10$$ dist.','$$p21$$ dist.','$$p32$$ dist.')
hold off
subplot(2,1,2);
plot(Tref,0*ref(1,:),'s-','Color',sstlightgray);
hold on
plot(Tref,0*ref(2,:),'s-','Color',sstlightgray);
plot(Tref,0*ref(3,:),'s-','Color',sstlightgray);
plot(TX,xc(2,:),'.--','Color',sstblue);
plot(TX,xc(4,:),'Color',sstblue);
plot(TX,xc(6,:),'s-','Color',sstblue);
plot(TX,xd(2,:),'.--','Color','r');
plot(TX,xd(4,:),'Color','r');
plot(TX,xd(6,:),'s-','Color','r');
plot(TX,x(2,:),'.--','Color','g');
plot(TX,x(4,:),'Color','g');
plot(TX,x(6,:),'s-','Color','g');
xlabel('$$t_k$$');
ylabel('$v_{i}$ [m/s]');
grid on
title('State $$v_{i}$$ evolution');
legend('$$r_1$$','$$r_2$$','$$r_3$$','$$v_1$$ cen.','$$v_2$$ cen.','$$v_3$$ cen.','$$v_1$$ decen.','$$v_2$$ decen.','$$v_3$$ decen.','$$v_1$$ dist.','$$v_2$$ dist.','$$v_3$$ dist.')
hold off

figure('Name','Input Evolution','NumberTitle','off');
subplot(2,1,1);
plot(Tref,ref(1,:),'s-','Color',sstlightgray);
hold on
plot(Tref,ref(2,:),'s-','Color',sstlightgray);
plot(Tref,ref(3,:),'s-','Color',sstlightgray);
plot(TX,xc(1,:),'.--','Color',sstblue);
plot(TX,xc(3,:),'Color',sstblue);
plot(TX,xc(5,:),'s-','Color',sstblue);
plot(TX,xd(1,:),'.--','Color','r');
plot(TX,xd(3,:),'Color','r');
plot(TX,xd(5,:),'s-','Color','r');
plot(TX,x(1,:),'.--','Color','g');
plot(TX,x(3,:),'Color','g');
plot(TX,x(5,:),'s-','Color','g');
xlabel('$$t_k$$');
ylabel('$p_{ij}$ [m]');
grid on
title('Output $$p_{ij}$$ evolution');
legend('$$r_1$$','$$r_2$$','$$r_3$$','$$p10$$ cen.','$$p21$$ cen.','$$p32$$ cen.','$$p10$$ decen.','$$p21$$ decen.','$$p32$$ decen.','$$p10$$ dist.','$$p21$$ dist.','$$p32$$ dist.')
hold off
subplot(2,1,2);
plot(TU,uc(1,:),'.--','Color',sstblue);
hold on
plot(TU,uc(2,:),'Color',sstblue);
plot(TU,uc(3,:),'s-','Color',sstblue);
plot(TU,ud(1,:),'.--','Color','r');
plot(TU,ud(2,:),'Color','r');
plot(TU,ud(3,:),'s-','Color','r');
plot(TU,u(1,:),'.--','Color','g');
plot(TU,u(2,:),'Color','g');
plot(TU,u(3,:),'s-','Color','g');
xlabel('$$t_k$$');
ylabel('$\theta_{i}$ [rad]');
grid on
title('Input evolution');
legend('$$\theta_1$$ cen.','$$\theta_2$$ cen.','$$\theta_3$$ cen.','$$\theta_1$$ decen.','$$\theta_2$$ decen.','$$\theta_3$$ decen.','$$\theta_1$$ dist.','$$\theta_2$$ dist.','$$\theta_3$$ dist.')
hold off



%% 2.2 Code: Constrained 3 player centralized, decentralized, distributed


np = 3;
% np = 10;


A22 = -2*beta*v1e/m;%m1 = m
A44 = -2*beta*v1e/m;%v2e = v1e && m2 = m
A66 = -2*beta*v1e/m;%v3e = v1e && m3 = m

Al = [0 -1  0  0  0  0 ;
      0 A22 0  0  0  0 ;
      0  1  0 -1  0  0 ;
      0  0  0 A44 0  0 ;
      0  0  0  1  0 -1 ;
      0  0  0  0  0 A66];

B21 = g*sec(theta1e)^2;
B42 = g*sec(theta1e)^2;%theta2e = theta1e
B63 = g*sec(theta1e)^2;%theta3e = theta1e


Bl = [ 0  0  0 ;
      B21 0  0 ;
       0  0  0 ;
       0 B42 0 ;
       0  0  0 ;
       0  0 B63];

Cl = [1 0 0 0 0 0 ;
      0 0 1 0 0 0 ;
      0 0 0 0 1 0 ];

E21 = 2*beta*we/m;%m1 = m
E41 = 2*beta*we/m;%m2 = m
E61 = 2*beta*we/m;%m3 = m

El = [ 0  1 ;
      E21 0 ;
       0  0 ;
      E41 0 ;
       0  0 ;
      E61 0 ];

Ad = eye(size(Al,1)) + Ts*Al;
Bd = Ts * Bl;
Ed = Ts * El;
Cd = Cl;

Ai = {Ad(1:2,1:2),Ad(2:4,2:4),Ad(4:6,4:6)};
Aidec = {Ad(1:2,1:2),Ad(3:4,3:4),Ad(5:6,5:6)};

Bij = { Bd(1:2,1), Bd(1:2,2), Bd(1:2,3)
        Bd(2:4,1), Bd(2:4,2), Bd(2:4,3)
        Bd(4:6,1), Bd(4:6,2), Bd(4:6,3) };
Bijdec = { Bd(1:2,1), Bd(1:2,2), Bd(1:2,3)
           Bd(3:4,1), Bd(3:4,2), Bd(3:4,3)
           Bd(5:6,1), Bd(5:6,2), Bd(5:6,3) };
Ci = {Cd(1,1:2),Cd(2,2:4),Cd(3,4:6)};
Cidec = {Cd(1,1:2),Cd(2,3:4),Cd(3,5:6)};
xd10 = [0 0]';
xd20 = xd10;
x0 = [xd10; xd10; xd10];

nu = size(Bd,2);
ny = size(Cd,1);
nx = size(Bd,1);
nxi = size(Bij{1,1},1);
nxi2 = size(Bij{2,2},1);
nui = size(Bij{1,1},2);
nxidec = size(Bijdec{1,1},1);
nuidec = size(Bijdec{1,1},2);
nud = size(Ed,2);
% cost parameters
N = 3;
Pi = 10;
Qi = 10;
Ri = .1;
alphai = 1;
% distributed steps parameters
w1 = 0.4;
w2 = 0.3;
w3 = 1-w1-w2;


% compute centralized tracking controller
P = blkdiag(Pi,Pi,Pi);
Q = blkdiag(Qi,Qi,Qi);
R = blkdiag(Ri,Ri,Ri);

[Fbc,Gbc,Qbc,Rbc,Fc,Gc,Hc] = GetYMats(Ad,Bd,Cd,N,P,Q,R);
[Rtc,Sxc,Syc] = GetYMPC(Fbc,Gbc,Qbc,Rbc);

% compute decentralized tracking controllers
[Fbd1,Gbd1,Qbd1,Rbd1,Fd1,Gd1,Hd1] = GetYMats(Aidec{1},Bijdec{1,1},Cidec{1},N,Pi,Qi,Ri);
[Rtd1,Sxd1,Syd1] = GetYMPC(Fbd1,Gbd1,Qbd1,Rbd1);
[Fbd2,Gbd2,Qbd2,Rbd2,Fd2,Gd2,Hd2] = GetYMats(Aidec{2},Bijdec{2,2},Cidec{2},N,Pi,Qi,Ri);
[Rtd2,Sxd2,Syd2] = GetYMPC(Fbd2,Gbd2,Qbd2,Rbd2);
[Fbd3,Gbd3,Qbd3,Rbd3,Fd3,Gd3,Hd3] = GetYMats(Aidec{3},Bijdec{3,3},Cidec{3},N,Pi,Qi,Ri);
[Rtd3,Sxd3,Syd3] = GetYMPC(Fbd3,Gbd3,Qbd3,Rbd3);

% compute distributed tracking controllers
[Fb,Gb,Qb,Rb,F,G,H] = GetYMats(Ai,Bij,Ci,N,Pi,Qi,Ri,alphai);
[Rt,Sx,Sy,Su] = GetYMPC(Fb,Gb,Qb,Rb);

% compute constraints matrices:

u_max = deg2rad(30);
x_max = [(5-0.5) ; 50];
x_min = [(0.2-0.5) ; -50];

% constraints centralized
U_max = kron(u_max,ones(N*nu,1));
X_max = repmat(x_max,(N+1)*nx/2,1);
X_min = repmat(x_min,(N+1)*nx/2,1);

mGc = cell2mat(Gc);
Mx = [-mGc;mGc];
Mu = [-eye(N*nu);eye(N*nu)];
wu = [U_max;U_max];
M = [Mu;Mx];

% constraints decentralized
Ui_max_dec = kron(u_max,ones(N*nuidec,1));
Xi_max_dec = repmat(x_max,(N+1)*(nxidec/2),1);
Xi_min_dec = repmat(x_min,(N+1)*(nxidec/2),1);

Mui_dec = [-eye(N*nuidec);eye(N*nuidec)];

Mxi1_dec = [-Gd1{1};Gd1{1}];
Mxi2_dec = [-Gd2{1};Gd2{1}];
Mxi3_dec = [-Gd3{1};Gd3{1}];

Mi1_dec = [Mui_dec];%;Mxi1_dec];
Mi2_dec = [Mui_dec];%;Mxi2_dec];
Mi3_dec = [Mui_dec];%;Mxi3_dec];
wui_dec = [Ui_max_dec;Ui_max_dec];

% constraints distributed
Ui_max = kron(u_max,ones(N*nui,1));
Xi1_max = repmat(x_max,(N+1)*(nxi/2),1);
Xi1_min = repmat(x_min,(N+1)*(nxi/2),1);

x2_3_max = [50 ; (5-0.5) ; 50];
x2_3_min = [-50 ; (0.2-0.5) ; -50];
Xi2_3_max = repmat(x2_3_max,(N+1)*(nxi2/3),1);
Xi2_3_min = repmat(x2_3_min,(N+1)*(nxi2/3),1);

Mui = [-eye(N*nui);eye(N*nui)];

Mxi1 = [- G{1} ; G{1}];
Mxi2 = [- G{2} ; G{2}];
Mxi3 = [- G{3} ; G{3}];

Mi1 = [Mui];% ; Mxi1];
Mi2 = [Mui];% ; Mxi2];
Mi3 = [Mui];% ; Mxi3];
wui = [Ui_max ; Ui_max];

% analyze stability of distributed law:
% L = [ (1-w1)*eye(size(L1))  , w1*L1
%       w2*L2                 , (1-w2)*eye(size(L2)) ];
% maxL = max(abs(eig(L))),
% Kb = [K11,K12;K21,K22];
% Kby = [K11y,K12y;K21y,K22y];
% K = E*(1-L)^(-1)*Kb;
% lbd_cl = eig(A+B*K),

%simulate controlled system:
nk = 200;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
% ref = [1*(Tref>=0)+1.5*(Tref>=40)-1*(Tref>=80);+0.5*(Tref>=0)+1.5*(Tref>=60)+0.5*(Tref>=90);1*(Tref>=0)+2*(Tref>=30)-0.5*(Tref>=70)]; % step reference of amplitude 1, starting at k=30
ref = [1*(Tref>=0)-0.5*(Tref>=100);5*(Tref>=0)-0.5*(Tref>=100);6*(Tref>=0)-0.5*(Tref>=100)];
uc = zeros(nu,nk);
ud = zeros(nu,nk);
u = zeros(nu,nk);
udd = [0*ones(nud-1,nk);0*ones(nud-1,nk)];
U1 = zeros(nui,N,nk);
U2 = zeros(nui,N,nk);
U3 = zeros(nui,N,nk);
xc = zeros(nx,nk+1);
xd = zeros(nx,nk+1);
x = zeros(nx,nk+1);
xc(:,1) = x0;
xd(:,1) = x0;
x(:,1) = x0;
xc(:,2) = x0;
xd(:,2) = x0;
x(:,2) = x0;

for k = 2:nk
    Yb = reshape(ref(:,k:k+N),[],1);
    Yb1 = ref(1,k:k+N)';
    Yb2 = ref(2,k:k+N)';
    Yb3 = ref(3,k:k+N)';
    

    % compute constraints centralized
    xk = xc(:,k);
    mFc = cell2mat(Fc);
    wx = [-X_min + mFc*xk;X_max - mFc*xk];
    w = [wu;wx];

    % centralized
    Stc = Sxc*xc(:,k)-Syc*Yb;
    
    [Uco,Jco,exitflag] = quadprog(Rtc,Stc,M,w);
    if exitflag~=1, error('Problems in centralized.'); end
    Ucopt(:,:,k) = reshape( Uco ,nu,N);
    uc(:,k) = Ucopt(:,1,k);

    % simulate system for centralized MPC
    xc(:,k+1) = Ad*xc(:,k) + Bd*uc(:,k) + Ed*udd(:,k); %simulate joint system


    % compute constraints decentralized
    wxi1_dec = [-Xi_min_dec + Fd1{1}*xd((1:2),k); Xi_max_dec - Fd1{1}*xd((1:2),k)];
    wi1_dec = [wui_dec];%;wxi1_dec];

    wxi2_dec = [-Xi_min_dec + Fd2{1}*xd((3:4),k); Xi_max_dec - Fd2{1}*xd((3:4),k)];
    wi2_dec = [wui_dec];%;wxi2_dec];

    wxi3_dec = [-Xi_min_dec + Fd3{1}*xd((5:6),k); Xi_max_dec - Fd3{1}*xd((5:6),k)];
    wi3_dec = [wui_dec];%;wxi3_dec];

    % decentralized
    Std1 = Sxd1*xd(1:2,k)-Syd1*Yb1;
    Std2 = Sxd2*xd(3:4,k)-Syd2*Yb2;
    Std3 = Sxd3*xd(5:6,k)-Syd3*Yb3;
    [Ud1o,Jd1o,ed1] = quadprog(Rtd1,Std1,Mi1_dec,wi1_dec);
    [Ud2o,Jd2o,ed2] = quadprog(Rtd2,Std2,Mi2_dec,wi2_dec);
    [Ud3o,Jd3o,ed3] = quadprog(Rtd3,Std3,Mi3_dec,wi3_dec);
    if any([ed1,ed2,ed3]~=1), error('Problems in decentralized MPC.'); end
    Ud1opt = reshape( Ud1o ,nui,N);
    Ud2opt = reshape( Ud2o ,nui,N);
    Ud3opt = reshape( Ud3o ,nui,N);
    ud(:,k) = [ Ud1opt(:,1) ; Ud2opt(:,1) ; Ud3opt(:,1) ];

    % simulate system for decentralized MPC
    xd(:,k+1) = Ad*xd(:,k) + Bd*ud(:,k) + Ed*udd(:,k); %simulate joint system
    
   
    % distributed MPC
    x1k = x(1:2,k);
    x2k = x(2:4,k);
    x3k = x(4:6,k);
    U1p = reshape( U1(:,:,k-1) ,[],1);
    U2p = reshape( U2(:,:,k-1) ,[],1);
    U3p = reshape( U3(:,:,k-1) ,[],1);

    % compute constraints distributed
    wxi1 = [-Xi1_min + F{1}*x1k; Xi1_max - F{1}*x1k];
    wi1 = [wui];%;wxi1];
 
    wxi2 = [-Xi2_3_min + F{2}*x2k; Xi2_3_max - F{2}*x2k];
    wi2 = [wui];%;wxi2];

    wxi3 = [-Xi2_3_min + F{3}*x3k; Xi2_3_max - F{3}*x3k];
    wi3 = [wui];%;wxi3];
    
    % Get optimal sequence for player 1
    St1 = Sx{1,1}*x1k - Sy{1,1}*Yb1;
    [U1o,J1o,exitflag] = quadprog(Rt{1},St1,Mi1,wi1);
    if exitflag~=1, error('Problems in player 1.'); end
    % Get optimal sequence for player 2
    St2 = Sx{2,1}*x1k - Sy{2,1}*Yb1 + Sx{2,2}*x2k - Sy{2,2}*Yb2 + Su{2,1}*U1p;
    [U2o,J2o,exitflag] = quadprog(Rt{2},St2,Mi2,wi2);
    if exitflag~=1, error('Problems in player 2.'); end
    % Get optimal sequence for player 3
    St3 = Sx{3,2}*x2k - Sy{3,2}*Yb2 + Sx{3,3}*x3k - Sy{3,3}*Yb3 + Su{3,2}*U2p;
    [U3o,J3o,exitflag] = quadprog(Rt{3},St3,Mi3,wi3);
    if exitflag~=1, error('Problems in player 3.'); end
    % convex step iterates
    U1pp = U1p;
    U2pp = U2p;
    U3pp = U3p;
    for p = 1:np
        U1pp = w1*U1o + (1-w1)*U1pp;
        U2pp = w2*U2o + (1-w2)*U2pp;
        U3pp = w3*U3o + (1-w3)*U3pp;
    end
    U1(:,:,k) = reshape( U1pp ,nui,N);
    U2(:,:,k) = reshape( U2pp ,nui,N);
    U3(:,:,k) = reshape( U3pp ,nui,N);
    
    u(:,k) = [ U1(:,1,k) ; U2(:,1,k) ; U3(:,1,k) ]; % apply first input to each player

    % simulate system for distributed MPC
    x(:,k+1) = Ad*x(:,k) + Bd*u(:,k) + Ed*udd(:,k); %simulate joint system
end

figure('Name','Simulation of both states for u = 0.1 rad for three drones (1)','NumberTitle','off');
subplot(2,1,1);
plot(Tref,ref(1,:),'s-','Color',sstlightgray);
hold on
plot(Tref,ref(2,:),'s-','Color',sstlightgray);
plot(Tref,ref(3,:),'s-','Color',sstlightgray);
plot(TX,xc(1,:),'.--','Color',sstblue);
plot(TX,xc(3,:),'Color',sstblue);
plot(TX,xc(5,:),'s-','Color',sstblue);
plot(TX,xd(1,:),'.--','Color','r');
plot(TX,xd(3,:),'Color','r');
plot(TX,xd(5,:),'s-','Color','r');
xlabel('$$t_k$$')
ylabel('$p_{ij}$ [m]');
grid on
title('State $$p_{ij}$$ evolution');
legend('$$r_1$$','$$r_2$$','$$r_3$$','$$p10$$ cen.','$$p21$$ cen.','$$p32$$ cen.','$$p10$$ decen.','$$p21$$ decen.','$$p32$$ decen.')
hold off

subplot(2,1,2);
plot(Tref,0*ref(1,:),'s-','Color',sstlightgray);
hold on
plot(Tref,0*ref(2,:),'s-','Color',sstlightgray);
plot(Tref,0*ref(3,:),'s-','Color',sstlightgray);
plot(TX,xc(2,:),'.--','Color',sstblue);
plot(TX,xc(4,:),'Color',sstblue);
plot(TX,xc(6,:),'s-','Color',sstblue);
plot(TX,xd(2,:),'.--','Color','r');
plot(TX,xd(4,:),'Color','r');
plot(TX,xd(6,:),'s-','Color','r');
xlabel('$$t_k$$')
ylabel('$v_{i}$ [m/s]');
grid on
title('State $$v_{i}$$ evolution');
legend('$$r_1$$','$$r_2$$','$$r_3$$','$$v_1$$ cen.','$$v_2$$ cen.','$$v_3$$ cen.','$$v_1$$ decen.','$$v_2$$ decen.','$$v_3$$ decen.')
hold off

figure('Name','Simulation of both states for u = 0.1 rad for three drones (2)','NumberTitle','off');
subplot(2,1,1);
plot(Tref,ref(1,:),'s-','Color',sstlightgray);
hold on
plot(Tref,ref(2,:),'s-','Color',sstlightgray);
plot(Tref,ref(3,:),'s-','Color',sstlightgray);
plot(TX,xc(1,:),'.--','Color',sstblue);
plot(TX,xc(3,:),'Color',sstblue);
plot(TX,xc(5,:),'s-','Color',sstblue);
plot(TX,xd(1,:),'.--','Color','r');
plot(TX,xd(3,:),'Color','r');
plot(TX,xd(5,:),'s-','Color','r');
plot(TX,x(1,:),'.--','Color','g');
plot(TX,x(3,:),'Color','g');
plot(TX,x(5,:),'s-','Color','g');
title('State $$p_{ij}$$ evolution');
xlabel('$$t_k$$');
ylabel('$p_{ij}$ [m]');
grid on
legend('$$r_1$$','$$r_2$$','$$r_3$$','$$p10$$ cen.','$$p21$$ cen.','$$p32$$ cen.','$$p10$$ decen.','$$p21$$ decen.','$$p32$$ decen.','$$p10$$ dist.','$$p21$$ dist.','$$p32$$ dist.')
hold off
subplot(2,1,2);
plot(Tref,0*ref(1,:),'s-','Color',sstlightgray);
hold on
plot(Tref,0*ref(2,:),'s-','Color',sstlightgray);
plot(Tref,0*ref(3,:),'s-','Color',sstlightgray);
plot(TX,xc(2,:),'.--','Color',sstblue);
plot(TX,xc(4,:),'Color',sstblue);
plot(TX,xc(6,:),'s-','Color',sstblue);
plot(TX,xd(2,:),'.--','Color','r');
plot(TX,xd(4,:),'Color','r');
plot(TX,xd(6,:),'s-','Color','r');
plot(TX,x(2,:),'.--','Color','g');
plot(TX,x(4,:),'Color','g');
plot(TX,x(6,:),'s-','Color','g');
xlabel('$$t_k$$');
ylabel('$v_{i}$ [m/s]');
grid on
title('State $$v_{i}$$ evolution');
legend('$$r_1$$','$$r_2$$','$$r_3$$','$$v_1$$ cen.','$$v_2$$ cen.','$$v_3$$ cen.','$$v_1$$ decen.','$$v_2$$ decen.','$$v_3$$ decen.','$$v_1$$ dist.','$$v_2$$ dist.','$$v_3$$ dist.')
hold off

figure('Name','Input Evolution','NumberTitle','off');
subplot(2,1,1);
plot(Tref,ref(1,:),'s-','Color',sstlightgray);
hold on
plot(Tref,ref(2,:),'s-','Color',sstlightgray);
plot(Tref,ref(3,:),'s-','Color',sstlightgray);
plot(TX,xc(1,:),'.--','Color',sstblue);
plot(TX,xc(3,:),'Color',sstblue);
plot(TX,xc(5,:),'s-','Color',sstblue);
plot(TX,xd(1,:),'.--','Color','r');
plot(TX,xd(3,:),'Color','r');
plot(TX,xd(5,:),'s-','Color','r');
plot(TX,x(1,:),'.--','Color','g');
plot(TX,x(3,:),'Color','g');
plot(TX,x(5,:),'s-','Color','g');
xlabel('$$t_k$$');
ylabel('$p_{ij}$ [m]');
grid on
title('Output $$p_{ij}$$ evolution');
legend('$$r_1$$','$$r_2$$','$$r_3$$','$$p10$$ cen.','$$p21$$ cen.','$$p32$$ cen.','$$p10$$ decen.','$$p21$$ decen.','$$p32$$ decen.','$$p10$$ dist.','$$p21$$ dist.','$$p32$$ dist.')
hold off
subplot(2,1,2);
plot(TU,uc(1,:),'.--','Color',sstblue);
hold on
plot(TU,uc(2,:),'Color',sstblue);
plot(TU,uc(3,:),'s-','Color',sstblue);
plot(TU,ud(1,:),'.--','Color','r');
plot(TU,ud(2,:),'Color','r');
plot(TU,ud(3,:),'s-','Color','r');
plot(TU,u(1,:),'.--','Color','g');
plot(TU,u(2,:),'Color','g');
plot(TU,u(3,:),'s-','Color','g');
xlabel('$$t_k$$');
ylabel('$\theta_{i}$ [rad]');
grid on
title('Input evolution');
legend('$$\theta_1$$ cen.','$$\theta_2$$ cen.','$$\theta_3$$ cen.','$$\theta_1$$ decen.','$$\theta_2$$ decen.','$$\theta_3$$ decen.','$$\theta_1$$ dist.','$$\theta_2$$ dist.','$$\theta_3$$ dist.')
hold off

%% 2.3 Code: 6 player centralized, decentralized, distributed


%np = 3;
np = 10;
% Centralized

nD = 6; % number of drones

A22 = -2*beta*v1e/m;%m1 = m
A44 = -2*beta*v1e/m;%v2e = v1e && m2 = m
A66 = -2*beta*v1e/m;%v3e = v1e && m3 = m
A88 = -2*beta*v1e/m;%v4e = v1e && m4 = m
Aaa = -2*beta*v1e/m;%v5e = v1e && m5 = m
Abb = -2*beta*v1e/m;%v6e = v1e && m6 = m

Al = [0 -1  0  0  0  0  0  0  0  0  0  0 ;
      0 A22 0  0  0  0  0  0  0  0  0  0 ;
      0  1  0 -1  0  0  0  0  0  0  0  0 ;
      0  0  0 A44 0  0  0  0  0  0  0  0 ;
      0  0  0  1  0 -1  0  0  0  0  0  0 ;
      0  0  0  0  0 A66 0  0  0  0  0  0 ;
      0  0  0  0  0  1  0 -1  0  0  0  0 ;
      0  0  0  0  0  0  0 A88 0  0  0  0 ;
      0  0  0  0  0  0  0  1  0 -1  0  0 ;
      0  0  0  0  0  0  0  0  0 Aaa 0  0 ;
      0  0  0  0  0  0  0  0  0  1  0 -1 ;
      0  0  0  0  0  0  0  0  0  0  0 Abb];

B21 = g*sec(theta1e)^2;
B42 = g*sec(theta1e)^2;%theta2e = theta1e
B63 = g*sec(theta1e)^2;%theta3e = theta1e
B84 = g*sec(theta1e)^2;%theta4e = theta1e
Ba6 = g*sec(theta1e)^2;%theta5e = theta1e
Bb8 = g*sec(theta1e)^2;%theta6e = theta1e

Bl = [ 0  0  0  0  0  0 ;
      B21 0  0  0  0  0 ;
       0  0  0  0  0  0 ;
       0 B42 0  0  0  0 ;
       0  0  0  0  0  0 ;
       0  0 B63 0  0  0 ;
       0  0  0  0  0  0 ;
       0  0  0 B84 0  0 ; 
       0  0  0  0  0  0 ;
       0  0  0  0 Ba6 0 ;
       0  0  0  0  0  0 ;
       0  0  0  0  0 Bb8];

Cl = [1 0 0 0 0 0 0 0 0 0 0 0 ;
      0 0 1 0 0 0 0 0 0 0 0 0 ;
      0 0 0 0 1 0 0 0 0 0 0 0 ;
      0 0 0 0 0 0 1 0 0 0 0 0 ;
      0 0 0 0 0 0 0 0 1 0 0 0 ;
      0 0 0 0 0 0 0 0 0 0 1 0 ];

E21 = 2*beta*we/m;%m1 = m
E41 = 2*beta*we/m;%m2 = m
E61 = 2*beta*we/m;%m3 = m
E81 = 2*beta*we/m;%m4 = m
Ea1 = 2*beta*we/m;%m5 = m
Eb1 = 2*beta*we/m;%m6 = m

El = [ 0  1 ;
      E21 0 ;
       0  0 ;
      E41 0 ;
       0  0 ;
      E61 0 ;
       0  0 ;
      E81 0 ;
       0  0 ;
      Ea1 0 ;
       0  0 ;
      Eb1 0 ];

Ad = eye(size(Al,1)) + Ts*Al;
Bd = Ts * Bl;
Ed = Ts * El;
Cd = Cl;

Aidec = {Ad(1:2,1:2),Ad(3:4,3:4),Ad(5:6,5:6),Ad(7:8,7:8),Ad(9:10,9:10),Ad(11:12,11:12)};

Bijdec = { Bd(1:2,1)  , Bd(1:2,2)  , Bd(1:2,3) , Bd(1:2,4)  , Bd(1:2,5)  , Bd(1:2,6)
        Bd(3:4,1)  , Bd(3:4,2)  , Bd(3:4,3) , Bd(3:4,4)  , Bd(3:4,5)  , Bd(3:4,6)
        Bd(5:6,1)  , Bd(5:6,2)  , Bd(5:6,3) , Bd(5:6,4)  , Bd(5:6,5)  , Bd(5:6,6)
        Bd(7:8,1)  , Bd(7:8,2)  , Bd(7:8,3) , Bd(7:8,4)  , Bd(7:8,5)  , Bd(7:8,6)
        Bd(9:10,1) , Bd(9:10,2) , Bd(9:10,3), Bd(9:10,4) , Bd(9:10,5) , Bd(9:10,6)
        Bd(11:12,1), Bd(11:12,2), Bd(11:12,3), Bd(11:12,4), Bd(11:12,5), Bd(11:12,6)};

Cidec = {Cd(1,1:2),Cd(2,3:4),Cd(3,5:6),Cd(4,7:8),Cd(5,9:10),Cd(6,11:12)};

%distributed matrices

Ai = {Ad(1:2,1:2),Ad(2:4,2:4),Ad(4:6,4:6),Ad(6:8,6:8),Ad(8:10,8:10),Ad(10:12,10:12)};

Bij = { Bd(1:2,1)  , Bd(1:2,2)  , Bd(1:2,3) , Bd(1:2,4)  , Bd(1:2,5)  , Bd(1:2,6)
        Bd(2:4,1)  , Bd(2:4,2)  , Bd(2:4,3) , Bd(2:4,4)  , Bd(2:4,5)  , Bd(2:4,6)
        Bd(4:6,1)  , Bd(4:6,2)  , Bd(4:6,3) , Bd(4:6,4)  , Bd(4:6,5)  , Bd(4:6,6)
        Bd(6:8,1)  , Bd(6:8,2)  , Bd(6:8,3) , Bd(6:8,4)  , Bd(6:8,5)  , Bd(6:8,6)
        Bd(8:10,1) , Bd(8:10,2) , Bd(8:10,3), Bd(8:10,4) , Bd(8:10,5) , Bd(8:10,6)
        Bd(10:12,1), Bd(10:12,2), Bd(10:12,3), Bd(10:12,4), Bd(10:12,5), Bd(10:12,6)};

Ci = {Cd(1,1:2),Cd(2,2:4),Cd(3,4:6),Cd(4,6:8),Cd(5,8:10),Cd(6,10:12)};


xd10 = [0 0]';
xd20 = xd10;
x0 = [xd10; xd10; xd10; xd10; xd10; xd10];

nu = size(Bd,2);
ny = size(Cd,1);
nx = size(Bd,1);
nxi = size(Bij{1,1},1);
nui = size(Bij{1,1},2);
nud = size(Ed,2);

% cost parameters
N = 3;
Pi = 10;
Qi = 10;
Ri = .1;
alphai = 1;
% distributed steps parameters
%wi = 1/nD;

w1 = 0.2;
w2 = 0.2;
w3 = 0.1;
w4 = 0.2;
w5 = 0.2;
w6 = 0.1;
% compute centralized tracking controller
P = blkdiag(Pi,Pi,Pi,Pi,Pi,Pi);
Q = blkdiag(Qi,Qi,Qi,Qi,Qi,Qi);
R = blkdiag(Ri,Ri,Ri,Ri,Ri,Ri);

[Fbc,Gbc,Qbc,Rbc] = GetYMats(Ad,Bd,Cd,N,P,Q,R);
[Rtc,Sxc,Syc] = GetYMPC(Fbc,Gbc,Qbc,Rbc);

% compute decentralized tracking controllers
[Fbd1,Gbd1,Qbd1,Rbd1] = GetYMats(Aidec{1},Bijdec{1,1},Cidec{1},N,Pi,Qi,Ri);
[Rtd1,Sxd1,Syd1,Sud1,K1,Ky1,L1] = GetYMPC(Fbd1,Gbd1,Qbd1,Rbd1);
[Fbd2,Gbd2,Qbd2,Rbd2] = GetYMats(Aidec{2},Bijdec{2,2},Cidec{2},N,Pi,Qi,Ri);
[Rtd2,Sxd2,Syd2,Sud2,K2,Ky2,L2] = GetYMPC(Fbd2,Gbd2,Qbd2,Rbd2);
[Fbd3,Gbd3,Qbd3,Rbd3] = GetYMats(Aidec{3},Bijdec{3,3},Cidec{3},N,Pi,Qi,Ri);
[Rtd3,Sxd3,Syd3,Sud3,K3,Ky3,L3] = GetYMPC(Fbd3,Gbd3,Qbd3,Rbd3);
[Fbd4,Gbd4,Qbd4,Rbd4] = GetYMats(Aidec{4},Bijdec{4,4},Cidec{4},N,Pi,Qi,Ri);
[Rtd4,Sxd4,Syd4,Sud4,K4,Ky4,L4] = GetYMPC(Fbd4,Gbd4,Qbd4,Rbd4);
[Fbd5,Gbd5,Qbd5,Rbd5] = GetYMats(Aidec{5},Bijdec{5,5},Cidec{5},N,Pi,Qi,Ri);
[Rtd5,Sxd5,Syd5,Sud5,K5,Ky5,L5] = GetYMPC(Fbd5,Gbd5,Qbd5,Rbd5);
[Fbd6,Gbd6,Qbd6,Rbd6] = GetYMats(Aidec{6},Bijdec{6,6},Cidec{6},N,Pi,Qi,Ri);
[Rtd6,Sxd6,Syd6,Sud6,K6,Ky6,L6] = GetYMPC(Fbd6,Gbd6,Qbd6,Rbd6);

% compute distributed tracking controllers
[Fb,Gb,Qb,Rb] = GetYMats(Ai,Bij,Ci,N,Pi,Qi,Ri,alphai);
[Rt,Sx,Sy,Su] = GetYMPC(Fb,Gb,Qb,Rb);

% compute constraints matrices:
%Ui_max = 0.5*ones(N*nui,1);
Mui = [];
wui = [];
%U_max = 0.5*ones(N*nu,1);
Mu = [];
wu = [];

% analyze stability of distributed law:
% L = [ (1-w1)*eye(size(L1))  , w1*L1
%       w2*L2                 , (1-w2)*eye(size(L2)) ];
% maxL = max(abs(eig(L))),
% Kb = [K11,K12;K21,K22];
% Kby = [K11y,K12y;K21y,K22y];
% K = E*(1-L)^(-1)*Kb;
% lbd_cl = eig(A+B*K),

%simulate controlled system:
nk = 200;
TU = 1:nk;
TX = 1:nk+1;
Tref = 1:nk+N;
% ref = [1*(Tref>=0)+1.5*(Tref>=40)-1*(Tref>=80);+0.5*(Tref>=0)+1.5*(Tref>=60)+0.5*(Tref>=90);1*(Tref>=0)+2*(Tref>=30)-0.5*(Tref>=70)]; % step reference of amplitude 1, starting at k=30
%ref = [1*(Tref>=0)-0.5*(Tref>=80);+0.5*(Tref>=20)+1.5*(Tref>=60);-0.5*(Tref>=30)+1.5*(Tref>=70);1*(Tref>=0)-0.5*(Tref>=80);+0.5*(Tref>=20)+1.5*(Tref>=60);-0.5*(Tref>=30)+1.5*(Tref>=70)];
%ref = [1*(Tref>=0);+0.5*(Tref>=20)+1.5*(Tref>=60);-0.5*(Tref>=30)+1.5*(Tref>=70);1*(Tref>=0)-0.5*(Tref>=80);+0.5*(Tref>=20)+1.5*(Tref>=60);-0.5*(Tref>=30)+1.5*(Tref>=70)];
ref = [0.5*(Tref>=10);1*(Tref>=20); 1.5*(Tref>=30);2*(Tref>=10);2.5*(Tref>=20); 3*(Tref>=30)];

uc = zeros(nu,nk);
ud = zeros(nu,nk);
u = zeros(nu,nk);
udd = [0*ones(nud-1,nk);0*ones(nud-1,nk)];
%udd = [[0*ones(nud-1,nk/2) 0*ones(nud-1,nk/4) 100*ones(nud-1,nk/4)];[0*ones(nud-1,nk/4) 1*ones(nud-1,nk/4) 0*ones(nud-1,nk/2)]];
U1 = zeros(nui,N,nk);
U2 = zeros(nui,N,nk);
U3 = zeros(nui,N,nk);
U4 = zeros(nui,N,nk);
U5 = zeros(nui,N,nk);
U6 = zeros(nui,N,nk);
xc = zeros(nx,nk+1);
xd = zeros(nx,nk+1);
x = zeros(nx,nk+1);
xc(:,1) = x0;
xd(:,1) = x0;
x(:,1) = x0;
xc(:,2) = x0;
xd(:,2) = x0;
x(:,2) = x0;
for k = 2:nk
    Yb = reshape(ref(:,k:k+N),[],1);
    Yb1 = ref(1,k:k+N)';
    Yb2 = ref(2,k:k+N)';
    Yb3 = ref(3,k:k+N)';
    Yb4 = ref(4,k:k+N)';
    Yb5 = ref(5,k:k+N)';
    Yb6 = ref(6,k:k+N)';
    
    % centralized
    Stc = Sxc*xc(:,k)-Syc*Yb;
    [Uco,Jco,exitflag] = quadprog(Rtc,Stc,Mu,wu);
    if exitflag~=1, error('Problems in centralized.'); end
    Ucopt(:,:,k) = reshape( Uco ,nu,N);
    uc(:,k) = Ucopt(:,1,k);

    % simulate system for centralized MPC
    xc(:,k+1) = Ad*xc(:,k) + Bd*uc(:,k) + Ed*udd(:,k); %simulate joint system

    % decentralized
    Std1 = Sxd1*xd(1:2,k)-Syd1*Yb1;
    Std2 = Sxd2*xd(3:4,k)-Syd2*Yb2;
    Std3 = Sxd3*xd(5:6,k)-Syd3*Yb3;
    Std4 = Sxd4*xd(7:8,k)-Syd4*Yb4;
    Std5 = Sxd5*xd(9:10,k)-Syd5*Yb5;
    Std6 = Sxd6*xd(11:12,k)-Syd6*Yb6;
    [Ud1o,Jd1o,ed1] = quadprog(Rtd1,Std1,Mui,wui);
    [Ud2o,Jd2o,ed2] = quadprog(Rtd2,Std2,Mui,wui);
    [Ud3o,Jd3o,ed3] = quadprog(Rtd3,Std3,Mui,wui);
    [Ud4o,Jd4o,ed4] = quadprog(Rtd4,Std4,Mui,wui);
    [Ud5o,Jd5o,ed5] = quadprog(Rtd5,Std5,Mui,wui);
    [Ud6o,Jd6o,ed6] = quadprog(Rtd6,Std6,Mui,wui);
    if any([ed1,ed2,ed3,ed4,ed5,ed6]~=1), error('Problems in decentralized MPC.'); end
    Ud1opt = reshape( Ud1o ,nui,N);
    Ud2opt = reshape( Ud2o ,nui,N);
    Ud3opt = reshape( Ud3o ,nui,N);
    Ud4opt = reshape( Ud4o ,nui,N);
    Ud5opt = reshape( Ud5o ,nui,N);
    Ud6opt = reshape( Ud6o ,nui,N);
    ud(:,k) = [ Ud1opt(:,1) ; Ud2opt(:,1) ; Ud3opt(:,1) ; Ud4opt(:,1) ; Ud5opt(:,1) ; Ud6opt(:,1) ];

    % simulate system for decentralized MPC
    xd(:,k+1) = Ad*xd(:,k) + Bd*ud(:,k) + Ed*udd(:,k); %simulate joint system
    
   
    % distributed MPC
    x1k = x(1:2,k);
    x2k = x(2:4,k);
    x3k = x(4:6,k);
    x4k = x(6:8,k);
    x5k = x(8:10,k);
    x6k = x(10:12,k);
    U1p = reshape( U1(:,:,k-1) ,[],1);
    U2p = reshape( U2(:,:,k-1) ,[],1);
    U3p = reshape( U3(:,:,k-1) ,[],1);
    U4p = reshape( U4(:,:,k-1) ,[],1);
    U5p = reshape( U5(:,:,k-1) ,[],1);
    U6p = reshape( U6(:,:,k-1) ,[],1);

    % Get optimal sequence for player 1
    St1 = Sx{1,1}*x1k - Sy{1,1}*Yb1;
    [U1o,J1o,exitflag] = quadprog(Rt{1},St1,Mui,wui);
    if exitflag~=1, error('Problems in player 1.'); end
    % Get optimal sequence for player 2
    St2 = Sx{2,1}*x1k - Sy{2,1}*Yb1 + Sx{2,2}*x2k - Sy{2,2}*Yb2 + Su{2,1}*U1p;
    [U2o,J2o,exitflag] = quadprog(Rt{2},St2,Mui,wui);
    if exitflag~=1, error('Problems in player 2.'); end
    % Get optimal sequence for player 3
    St3 = Sx{3,2}*x2k - Sy{3,2}*Yb2 + Sx{3,3}*x3k - Sy{3,3}*Yb3 + Su{3,2}*U2p;
    [U3o,J3o,exitflag] = quadprog(Rt{3},St3,Mui,wui);
    if exitflag~=1, error('Problems in player 3.'); end
     % Get optimal sequence for player 4
    St4 = Sx{4,3}*x3k - Sy{4,3}*Yb3 + Sx{4,4}*x4k - Sy{4,4}*Yb4 + Su{4,3}*U3p;
    [U4o,J4o,exitflag] = quadprog(Rt{4},St4,Mui,wui);
    if exitflag~=1, error('Problems in player 4.'); end
    % Get optimal sequence for player 5
    St5 = Sx{5,4}*x4k - Sy{5,4}*Yb4 + Sx{5,5}*x5k - Sy{5,5}*Yb5 + Su{5,4}*U4p;
    [U5o,J5o,exitflag] = quadprog(Rt{5},St5,Mui,wui);
    if exitflag~=1, error('Problems in player 5.'); end
    % Get optimal sequence for player 6
    St6 = Sx{6,5}*x5k - Sy{6,5}*Yb5 + Sx{6,6}*x6k - Sy{6,6}*Yb6 + Su{6,5}*U5p;
    [U6o,J6o,exitflag] = quadprog(Rt{6},St6,Mui,wui);
    if exitflag~=1, error('Problems in player 6.'); end
    % convex step iterates
    U1pp = U1p;
    U2pp = U2p;
    U3pp = U3p;
    U4pp = U4p;
    U5pp = U5p;
    U6pp = U6p;
    for p = 1:np
        U1pp = w1*U1o + (1-w1)*U1pp;
        U2pp = w2*U2o + (1-w2)*U2pp;
        U3pp = w3*U3o + (1-w3)*U3pp;
        U4pp = w4*U4o + (1-w4)*U4pp;
        U5pp = w5*U5o + (1-w5)*U5pp;
        U6pp = w6*U6o + (1-w6)*U6pp;
    end
    U1(:,:,k) = reshape( U1pp ,nui,N);
    U2(:,:,k) = reshape( U2pp ,nui,N);
    U3(:,:,k) = reshape( U3pp ,nui,N);
    U4(:,:,k) = reshape( U4pp ,nui,N);
    U5(:,:,k) = reshape( U5pp ,nui,N);
    U6(:,:,k) = reshape( U6pp ,nui,N);
    
    u(:,k) = [ U1(:,1,k) ; U2(:,1,k) ; U3(:,1,k); U4(:,1,k) ; U5(:,1,k) ; U6(:,1,k) ]; % apply first input to each player

    % simulate system for distributed MPC
    x(:,k+1) = Ad*x(:,k) + Bd*u(:,k) + Ed*udd(:,k); %simulate joint system
end

figure('Name','Simulation of both states for u = 0.1 rad for three drones (1)','NumberTitle','off');
subplot(2,1,1);
plot(Tref,ref(1,:),'s-','Color',sstlightgray);
hold on
plot(Tref,ref(2,:),'s-','Color',sstlightgray);
plot(Tref,ref(3,:),'s-','Color',sstlightgray);
plot(Tref,ref(4,:),'s-','Color',sstlightgray);
plot(Tref,ref(5,:),'s-','Color',sstlightgray);
plot(Tref,ref(6,:),'s-','Color',sstlightgray);
plot(TX,xc(1,:),'.--','Color',sstblue);
plot(TX,xc(3,:),'Color',sstblue);
plot(TX,xc(5,:),'s-','Color',sstblue);
plot(TX,xc(7,:),'.--','Color',sstdarkblue);
plot(TX,xc(9,:),'Color',sstdarkblue);
plot(TX,xc(11,:),'s-','Color',sstdarkblue);
plot(TX,xd(1,:),'.--','Color','r');
plot(TX,xd(3,:),'Color','r');
plot(TX,xd(5,:),'s-','Color','r');
plot(TX,xd(7,:),'.--','Color','r');
plot(TX,xd(9,:),'Color','r');
plot(TX,xd(11,:),'s-','Color','r');
xlabel('$$t_k$$')
ylabel('$p_{ij}$ [m]');
grid on
title('State $$p_{ij}$$ evolution');
legend('$$r_1$$','$$r_2$$','$$r_3$$','$$r_4$$','$$r_5$$','$$r_6$$','$$p10$$ cen.','$$p21$$ cen.','$$p32$$ cen.','$$p43$$ cen.','$$p54$$ cen.','$$p65$$ cen.','$$p10$$ decen.','$$p21$$ decen.','$$p32$$ decen.','$$p43$$ decen.','$$p54$$ decen.','$$p65$$ decen.')
hold off
subplot(2,1,2);
plot(Tref,0*ref(1,:),'s-','Color',sstlightgray);
hold on
plot(Tref,0*ref(2,:),'s-','Color',sstlightgray);
plot(Tref,0*ref(3,:),'s-','Color',sstlightgray);
plot(TX,xc(2,:),'.--','Color',sstblue);
plot(TX,xc(4,:),'Color',sstblue);
plot(TX,xc(6,:),'s-','Color',sstblue);
plot(TX,xc(8,:),'.--','Color',sstdarkblue);
plot(TX,xc(10,:),'Color',sstdarkblue);
plot(TX,xc(12,:),'s-','Color',sstdarkblue);
plot(TX,xd(2,:),'.--','Color','r');
plot(TX,xd(4,:),'Color','r');
plot(TX,xd(6,:),'s-','Color','r');
plot(TX,xd(8,:),'.--','Color','r');
plot(TX,xd(10,:),'Color','r');
plot(TX,xd(12,:),'s-','Color','r');
xlabel('$$t_k$$')
ylabel('$v_{i}$ [m/s]');
grid on
title('State $$v_{i}$$ evolution');
legend('$$r_1$$','$$r_2$$','$$r_3$$','$$v_1$$ cen.','$$v_2$$ cen.','$$v_3$$ cen.','$$v_4$$ cen.','$$v_5$$ cen.','$$v_6$$ cen.','$$v_1$$ decen.','$$v_2$$ decen.','$$v_3$$ decen.','$$v_4$$ decen.','$$v_5$$ decen.','$$v_6$$ decen.')
hold off

figure('Name','Simulation of both states for u = 0.1 rad for three drones (2)','NumberTitle','off');
subplot(2,1,1);
plot(Tref,ref(1,:),'s-','Color',sstlightgray);
hold on
plot(Tref,ref(2,:),'s-','Color',sstlightgray);
plot(Tref,ref(3,:),'s-','Color',sstlightgray);
plot(Tref,ref(4,:),'s-','Color',sstlightgray);
plot(Tref,ref(5,:),'s-','Color',sstlightgray);
plot(Tref,ref(6,:),'s-','Color',sstlightgray);
plot(TX,xc(1,:),'.--','Color',sstblue);
plot(TX,xc(3,:),'Color',sstblue);
plot(TX,xc(5,:),'s-','Color',sstblue);
plot(TX,xc(7,:),'.--','Color',sstblue);
plot(TX,xc(9,:),'Color',sstblue);
plot(TX,xc(11,:),'s-','Color',sstblue);
plot(TX,xd(1,:),'.--','Color','r');
plot(TX,xd(3,:),'Color','r');
plot(TX,xd(5,:),'s-','Color','r');
plot(TX,xd(7,:),'.--','Color','r');
plot(TX,xd(9,:),'Color','r');
plot(TX,xd(11,:),'s-','Color','r');
plot(TX,x(1,:),'.--','Color','g');
plot(TX,x(3,:),'Color','g');
plot(TX,x(5,:),'s-','Color','g');
plot(TX,x(7,:),'.--','Color','g');
plot(TX,x(9,:),'Color','g');
plot(TX,x(11,:),'s-','Color','g');
title('State $$p_{ij}$$ evolution');
xlabel('$$t_k$$');
ylabel('$p_{ij}$ [m]');
grid on
legend('$$r_1$$','$$r_2$$','$$r_3$$','$$r_4$$','$$r_5$$','$$r_6$$','$$p10$$ cen.','$$p21$$ cen.','$$p32$$ cen.','$$p43$$ cen.','$$p54$$ cen.','$$p65$$ cen.','$$p10$$ decen.','$$p21$$ decen.','$$p32$$ decen.','$$p43$$ decen.','$$p54$$ decen.','$$p65$$ decen.','$$p10$$ dist.','$$p21$$ dist.','$$p32$$ dist.','$$p43$$ dist.','$$p54$$ dist.','$$p65$$ dist.')
hold off
subplot(2,1,2);
plot(Tref,0*ref(1,:),'s-','Color',sstlightgray);
hold on
plot(Tref,0*ref(2,:),'s-','Color',sstlightgray);
plot(Tref,0*ref(3,:),'s-','Color',sstlightgray);
plot(TX,xc(2,:),'.--','Color',sstblue);
plot(TX,xc(4,:),'Color',sstblue);
plot(TX,xc(6,:),'s-','Color',sstblue);
plot(TX,xc(8,:),'.--','Color',sstblue);
plot(TX,xc(10,:),'Color',sstblue);
plot(TX,xc(12,:),'s-','Color',sstblue);
plot(TX,xd(2,:),'.--','Color','r');
plot(TX,xd(4,:),'Color','r');
plot(TX,xd(6,:),'s-','Color','r');
plot(TX,xd(8,:),'.--','Color','r');
plot(TX,xd(10,:),'Color','r');
plot(TX,xd(12,:),'s-','Color','r');
plot(TX,x(2,:),'.--','Color','g');
plot(TX,x(4,:),'Color','g');
plot(TX,x(6,:),'s-','Color','g');
plot(TX,x(8,:),'.--','Color','g');
plot(TX,x(10,:),'Color','g');
plot(TX,x(12,:),'s-','Color','g');
xlabel('$$t_k$$');
ylabel('$v_{i}$ [m/s]');
grid on
title('State $$v_{i}$$ evolution');
legend('$$r_1$$','$$r_2$$','$$r_3$$','$$v_1$$ cen.','$$v_2$$ cen.','$$v_3$$ cen.','$$v_4$$ cen.','$$v_5$$ cen.','$$v_6$$ cen.','$$v_1$$ decen.','$$v_2$$ decen.','$$v_3$$ decen.','$$v_4$$ decen.','$$v_5$$ decen.','$$v_6$$ decen.','$$v_1$$ dist.','$$v_2$$ dist.','$$v_3$$ dist.','$$v_4$$ dist.','$$v_5$$ dist.','$$v_6$$ dist.')
hold off

figure('Name','Input Evolution','NumberTitle','off');
subplot(2,1,1);
plot(Tref,ref(1,:),'s-','Color',sstlightgray);
hold on
plot(Tref,ref(2,:),'s-','Color',sstlightgray);
plot(Tref,ref(3,:),'s-','Color',sstlightgray);
plot(Tref,ref(4,:),'s-','Color',sstlightgray);
plot(Tref,ref(5,:),'s-','Color',sstlightgray);
plot(Tref,ref(6,:),'s-','Color',sstlightgray);
plot(TX,xc(1,:),'.--','Color',sstblue);
plot(TX,xc(3,:),'Color',sstblue);
plot(TX,xc(5,:),'s-','Color',sstblue);
plot(TX,xc(7,:),'.--','Color',sstdarkblue);
plot(TX,xc(9,:),'Color',sstdarkblue);
plot(TX,xc(11,:),'s-','Color',sstdarkblue);
plot(TX,xd(1,:),'.--','Color','r');
plot(TX,xd(3,:),'Color','r');
plot(TX,xd(5,:),'s-','Color','r');
plot(TX,xd(7,:),'.--','Color','r');
plot(TX,xd(9,:),'Color','r');
plot(TX,xd(11,:),'s-','Color','r');
plot(TX,x(1,:),'.--','Color','g');
plot(TX,x(3,:),'Color','g');
plot(TX,x(5,:),'s-','Color','g');
plot(TX,x(7,:),'.--','Color','g');
plot(TX,x(9,:),'Color','g');
plot(TX,x(11,:),'s-','Color','g');
xlabel('$$t_k$$')
ylabel('$p_{ij}$ [m]');
grid on
title('Output $$p_{ij}$$ evolution');
legend('$$r_1$$','$$r_2$$','$$r_3$$','$$r_4$$','$$r_5$$','$$r_6$$','$$p10$$ cen.','$$p21$$ cen.','$$p32$$ cen.','$$p43$$ cen.','$$p54$$ cen.','$$p65$$ cen.','$$p10$$ decen.','$$p21$$ decen.','$$p32$$ decen.','$$p43$$ decen.','$$p54$$ decen.','$$p65$$ decen.','$$p10$$ dist.','$$p21$$ dist.','$$p32$$ dist.','$$p43$$ dist.','$$p54$$ dist.','$$p65$$ dist.')
hold off
subplot(2,1,2);
plot(TU,uc(1,:),'.--','Color',sstblue);
hold on
plot(TU,uc(2,:),'Color',sstblue);
plot(TU,uc(3,:),'s-','Color',sstblue);
plot(TU,uc(4,:),'Color',sstblue);
plot(TU,uc(5,:),'s-','Color',sstblue);
plot(TU,uc(6,:),'Color',sstblue);
plot(TU,ud(1,:),'.--','Color','r');
plot(TU,ud(2,:),'Color','r');
plot(TU,ud(3,:),'s-','Color','r');
plot(TU,ud(4,:),'.--','Color','r');
plot(TU,ud(5,:),'Color','r');
plot(TU,ud(6,:),'s-','Color','r');
plot(TU,u(1,:),'.--','Color','g');
plot(TU,u(2,:),'Color','g');
plot(TU,u(3,:),'s-','Color','g');
plot(TU,u(4,:),'.--','Color','g');
plot(TU,u(5,:),'Color','g');
plot(TU,u(6,:),'s-','Color','g');
xlabel('$$t_k$$');
ylabel('$\theta_{i}$ [rad]');
grid on
title('Input evolution');
legend('$$\theta_1$$ cen.','$$\theta_2$$ cen.','$$\theta_3$$ cen.','$$\theta_4$$ cen.','$$\theta_5$$ cen.','$$\theta_6$$ cen.','$$\theta_1$$ decen.','$$\theta_2$$ decen.','$$\theta_3$$ decen.','$$\theta_4$$ decen.','$$\theta_5$$ decen.','$$\theta_6$$ decen.','$$\theta_1$$ dist.','$$\theta_2$$ dist.','$$\theta_3$$ dist.','$$\theta_4$$ dist.','$$\theta_5$$ dist.','$$\theta_6$$ dist.')
hold off


% simulation for paper
% 
% figure('Name','Simulation of both states for u = 0.1 rad for three drones (2)','NumberTitle','off');
% subplot(2,1,1);
% plot(Tref,ref(1,:),'-','Color',sstlightgray);
% hold on
% plot(Tref,ref(2,:),'-','Color',sstlightgray);
% plot(Tref,ref(3,:),'-','Color',sstlightgray);
% plot(Tref,ref(4,:),'-','Color',sstlightgray);
% plot(Tref,ref(5,:),'-','Color',sstlightgray);
% plot(Tref,ref(6,:),'-','Color',sstlightgray);
% plot(TX,xc(1,:),'-','Color',sstblue);
% plot(TX,xc(3,:),'-','Color',sstblue);
% plot(TX,xc(5,:),'-','Color',sstblue);
% plot(TX,xc(7,:),'-','Color',sstblue);
% plot(TX,xc(9,:),'-','Color',sstblue);
% plot(TX,xc(11,:),'-','Color',sstblue);
% plot(TX,xd(1,:),'.--','Color',sstgreen);
% plot(TX,xd(3,:),'.--','Color',sstgreen);
% plot(TX,xd(5,:),'.--','Color',sstgreen);
% plot(TX,xd(7,:),'.--','Color',sstgreen);
% plot(TX,xd(9,:),'.--','Color',sstgreen);
% plot(TX,xd(11,:),'.--','Color',sstgreen);
% plot(TX,x(1,:),':','Color','r');
% plot(TX,x(3,:),':','Color','r');
% plot(TX,x(5,:),':','Color','r');
% plot(TX,x(7,:),':','Color','r');
% plot(TX,x(9,:),':','Color','r');
% plot(TX,x(11,:),':','Color','r');
% title('State $$p_{ij}$$ evolution');
% xlabel('$$t_k$$');
% ylabel('$p_{ij}$ [m]');
% grid on
% legend('$$r_1$$','$$r_2$$','$$r_3$$','$$r_4$$','$$r_5$$','$$r_6$$','$$p10$$ cen.','$$p21$$ cen.','$$p32$$ cen.','$$p43$$ cen.','$$p54$$ cen.','$$p65$$ cen.','$$p10$$ decen.','$$p21$$ decen.','$$p32$$ decen.','$$p43$$ decen.','$$p54$$ decen.','$$p65$$ decen.','$$p10$$ dist.','$$p21$$ dist.','$$p32$$ dist.','$$p43$$ dist.','$$p54$$ dist.','$$p65$$ dist.')
% hold off
% subplot(2,1,2);
% plot(Tref,0*ref(1,:),'-','Color',sstlightgray);
% hold on
% plot(Tref,0*ref(2,:),'-','Color',sstlightgray);
% plot(Tref,0*ref(3,:),'-','Color',sstlightgray);
% plot(TX,xc(2,:),'.--','Color',sstblue);
% plot(TX,xc(4,:),'.--','Color',sstblue);
% plot(TX,xc(6,:),'.--','Color',sstblue);
% plot(TX,xc(8,:),'.--','Color',sstblue);
% plot(TX,xc(10,:),'.--','Color',sstblue);
% plot(TX,xc(12,:),'.--','Color',sstblue);
% plot(TX,xd(2,:),'s-','Color',sstdarkgreen);
% plot(TX,xd(4,:),'s-','Color',sstdarkgreen);
% plot(TX,xd(6,:),'s-','Color',sstdarkgreen);
% plot(TX,xd(8,:),'s-','Color',sstdarkgreen);
% plot(TX,xd(10,:),'s-','Color',sstdarkgreen);
% plot(TX,xd(12,:),'s-','Color',sstdarkgreen);
% plot(TX,x(2,:),'Color','r');
% plot(TX,x(4,:),'Color','r');
% plot(TX,x(6,:),'Color','r');
% plot(TX,x(8,:),'Color','r');
% plot(TX,x(10,:),'Color','r');
% plot(TX,x(12,:),'Color','r');
% xlabel('$$t_k$$');
% ylabel('$v_{i}$ [m/s]');
% grid on
% title('State $$v_{i}$$ evolution');
% legend('$$r_1$$','$$r_2$$','$$r_3$$','$$v_1$$ cen.','$$v_2$$ cen.','$$v_3$$ cen.','$$v_4$$ cen.','$$v_5$$ cen.','$$v_6$$ cen.','$$v_1$$ decen.','$$v_2$$ decen.','$$v_3$$ decen.','$$v_4$$ decen.','$$v_5$$ decen.','$$v_6$$ decen.','$$v_1$$ dist.','$$v_2$$ dist.','$$v_3$$ dist.','$$v_4$$ dist.','$$v_5$$ dist.','$$v_6$$ dist.')
% hold off

