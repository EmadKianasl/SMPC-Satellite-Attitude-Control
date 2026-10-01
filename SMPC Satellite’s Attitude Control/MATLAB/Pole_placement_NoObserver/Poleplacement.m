clc;clear;
close all

%% system parameters
j1=1;
j2=0.1;
k=0.091;
b=0.0036;
tfinal = 40; % final time step
stime = 0.1; % size of the time step
taim = 0:stime:tfinal;
samps = size(taim);
u1 = 1; % instrument package control torque
u1=u1+zeros(1,samps(1,2));
%% system matrices

J=[j1 0;0 j2];
K=[k -k;-k k];
B=[b -b;-b b];

%% state space matrices

A=[0 1 0 0;-k/j2 -b/j2 k/j2 b/j2;0 0 0 01;k/j1 b/j1 -k/j1 -b/j1];
B=[0;0;0;1/j1];
C=[1 0 0 0];
D=zeros(1,1);
sys=ss(A,B,C,D);
[Ad,Bd,Cd,Dd]=c2dm(A,B,C,D,0.1,'zoh'); %Discrete systems
%Dsys = ss(Ad,Bd,Cd,Dd,0.1);
transferfunc=tf([0.036 0.9],[1 0.04 1 0 0]);
numGG=conv([0.036 0.9],[0.036 0.9]);
denGG=conv([1 0.04 1 0 0],[1 0.04 1 0 0]);
sysGG=tf(numGG,denGG);
rlocus(sysGG,'blue');

%% state space matrices (Modal)
Fm=[0 1 0 0;-km2/jm2 -bm2/jm2 0 0;0 0 0 1;0 0 -km1/jm1 -bm1/jm1];
Gm=[0;1/jm2;0;1/jm1];
Hm=[1 0 0 0]*phi2;
Dm=zeros(1,1);
sysm=ss(Fm,Gm,Hm,Dm);
figure(4)
lsim(sysm,ones(size(taim)),taim,'b'); 

%% Pole placement
p=[-0.45+0.34j -0.45-0.34j -0.15+1.05j -0.15-1.05j]; % desired poles location corresponding to wbw=0.5 
%K = acker(A,B,p)
%K_cont=place(A,B,p);
K_cont=[-0.2788 0.0546 0.6814 1.1655];
Nbar = rscale(A,B,C,D,K_cont);
%% Closed loop response
%CLT=K*[0;0; 0; 1]*transferfunc/(1+K*[0;0; 0; 1]*transferfunc);
%zpk(CLT);
%s=tf('s');
%numCLT = 0.041774*(s+25);
%denCLT = (s^2+1.042*s+1.022)*(s^2-1.002*s+1.022);
%CLT=numCLT/denCLT;
%zpk(CLT);

%% Closed loop response
sys_cl = ss(A-B*K_cont,B,C,0);
figure(1)
lsim(sys_cl,u1*Nbar,taim,[0 0 0 0]);

CLT=tf(sys_cl);
figure(2)
stepplot(CLT*Nbar);
transferfunc_CL=tf(sys_cl);

%% The Bode plot of the SRL controller
syms r
SRL_controller=K_cont*inv(r*eye(4,4)-A)*B;
SRL_controller=simplify(SRL_controller);
%pretty(SRL_controller);
SRL_controller=tf([5827500 3626618 5623923 1831830],[5000000 198000 5005000 0 0]);
figure(3)
bode(SRL_controller);
