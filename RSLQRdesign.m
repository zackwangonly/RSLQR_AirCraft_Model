clear all
close all
clc

format short e

disp('****************************** Program Start ****************');

rtd  = 180/pi;
dt = 0.002;
t = 0:dt:2;
w = logspace(-1,3,500);
dd=0.:.001:2*pi;
xx1=cos(dd)-1;yy1=sin(dd);
grav = 9.8

%%%%%%%%%%%% PLANT DATA %%%%%%%%%%%%%%%%
% Here Za = Za/V  Zd = Zd/V
% Must multiply by V to get accel      %
% Az = V(Za alp + Zd dele)             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Za = -1.3046;
Zd = -0.2142;
Ma = 47.7109;
Md = -104.8346;

V = 886.78;%fps

% Plant States(no actuator)
Ap = [Za 1;Ma 0];
Bp = [Za;Md]
Cp=[Za*V 0; 1 0;0 1];
Dp=[Zd*V;0;0]

% design a RSLQR command to track Az using state feeback controller 
% form the wiggle system: 
% state vector: [int(e) (fps), AOA (rad), pitch rate q (rps)]'  
Aw = [0.  Za 0.; 0.  Za*V 1.; 0.   Ma    0.]; 
Bw = [   Zd; Zd*V; Md]; 

% Setup range of penalties for the LQR 
Q=0.*Aw; %sets up a matrix Q that is the size of Aw and is all 0s 
R=1;
%qq is a vector which each element scales the LQR Q matrix
qq=logspace(-6,0,500); %these are the varying values of q11 
% The only penalty is the 1,1 element on the error state 
% Desing the gains
for ii = 1 : numel(qq)
    Q(1,1)=qq(ii);
    [Kx_lqr,~,~]=lqr(Aw,Bw,Q,R);
% populate the controller matrices
    Ac =  0.; 
    Bc1 = [1. 0. 0. ]; 
    Bc2 =  -1; 
    Cc = -Kx_lqr(1); 
    Dc1 = [0. -Kx_lqr(2:3) ];
    Dc2 = 0.; 
% Form the closed loop system 
    Z = inv(eye(size(Dc1*Dp))-Dc1*Dp); 
    Acl = [     (Ap+Bp*Z*Dc1*Cp)  (Bp*Z*Cc);
        (Bc1*(Cp+Dp*Z*Dc1*Cp))  (Ac+Bc1*Dp*Z*Cc)]; 
    Bcl = [       Bp*Z*Dc2; 
        (Bc2+Bc1*Dp*Z*Dc2)]; 
    Ccl = [(Cp+Dp*Z*Dc1*Cp) (Dp*Z*Cc)]; 
    Dcl =(Dp*Z*Dc2); 
    sys_cl = ss(Acl,Bcl,Ccl,Dcl);
% Frequency Domain Analysis 
% SS model of loop gain at the plant input Lu 
    A_Lu = [ Ap 0.*Bp*Cc;  Bc1*Cp Ac]; 
    B_Lu = [ Bp; Bc1*Dp]; 
    C_Lu = -[ Dc1*Cp Cc];%change sign for loop gain 
    D_Lu = -[ Dc1*Dp]; 
    sys_Lu = ss(A_Lu,B_Lu,C_Lu,D_Lu); 
    magdb = 20*log10(abs(squeeze(freqresp(sys_Lu,w)))); 
    wc = w*magdb; % LGCF, assumes Lu is a scalar 
    sr = sigma(sys_Lu,w,3); % Stability Robustness 
    srmin = min(abs(sr)); 
    rd = sigma(sys_Lu,w,2); % Return Difference 
    rdmin = min(abs(rd)); 
    T  = freqresp(sys_cl,w); % Complementary Sensitivity 
    S = 1 - T; % Sensitivity 
    T_st(ii,:) = 20*log10(abs(squeeze(T(1,1,:)))); 
    S_st(ii,:) = 20*log10(abs(squeeze(S(1,1,:)))); 
    Tmax = max(T_st(ii,:)); % Inf Norm of T in dB 
    Smax = max(S_st(ii,:)); % Inf Norm of S in dB
% Time Domain Analysis 
    y = step(sys_cl,t); 
    az = y(:,1); %  acceleration (fps2) 
    aze = abs(ones(size(az))-az);  % error for az 
    taur = 0.; taus= 0.; % rise time and settling time 
    fv = aze(numel(aze)); % final value of the error 
    e_n = aze - fv*ones(size(aze)) - 0.36*ones(size(aze)); 
    e_n1 = abs(e_n) + e_n; 
    taur = t*e_n1; % rise time 
    e_n = aze - fv*ones(size(aze)) - 0.05*ones(size(aze)); 
    e_n1 = abs(e_n) + e_n; 
    taus = t*e_n1; % settling time 
    azmin = abs(min(az))*100; % undershoot 
    azmax = (abs(max(az))-1)*100; % overshoot 
    dmax = max(abs(y(:,2)))*rtd*grav; % compute in per g commanded 
    ddmax = max(abs(y(:,3)))*rtd*grav; 
    metric=[qq(ii) rdmin srmin wc taur taus azmin azmax dmax ddmax Tmax] 
    data(ii,:) = metric; 
    az_st(ii,:)     = az'; 
    del_st(ii,:)    = rtd*y(:,2); 
    deldot_st(ii,:) = rtd*y(:,3);
end
i = 1:1:500;
figure 
plot(2-t,az);
hold on;


