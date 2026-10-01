clc;clear;
close all

%% system parameters
load('results.mat')
j1=1;
j2=0.1;
k=0.091;
b=0.0036;
tfinal = 40; % final time step
stime = 0.1; % size of the time step
taim = 0:stime:tfinal;
samps = size(taim);
s=tf('s');
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

%% initial condition 
%x0=[0.2;0;0;0];
%phi=s*eye(4,4)-A;
%transf=C*inv(phi)*B+D;
%transf0=C*inv(phi)*x0;
%% transfer function 
%transferfunc=tf(sys);
transferfunc=tf([0.036 0.9],[1 0.04 1 0 0]);
%% open loop
figure(1);
rlocus(transferfunc);
figure(2);
bode(transferfunc)

%% closed loop   pd--> c1=0.25*G
%C1=tf([0.5 0.25],[1]);
C1=0.25;
figure(3)
rlocus(C1*transferfunc);
figure(4)
bode(C1*transferfunc)
%% closed loop   pd--> c1=0.25(2s+1)
DC1=tf([0.5 0.25],[1]);
figure(5)
rlocus(DC1*transferfunc);
figure(6)
bode(DC1*transferfunc)

%% closed loop   pd--> c1=0.25(2s+1)*[((s/0.9)^2+1)/((s/25)+1)^2]
DC1=tf([0.5 0.25],[1]);
%DC2=tf([(1/0.81) 0 1],[1/625 2/25 1]);
DC2=tf([62500 0 50625],[81 4050 26325]);
DC3=DC1*DC2;
figure(6)
rlocus(DC2);
figure(7)
bode(DC2)

%% torque
u1 = 1; % instrument package control torque
u1=u1+zeros(1,samps(1,2));

u2 = 0; % instrument package control torque
u2=u2+zeros(1,samps(1,2));
%% closed loop entire
CLT=(DC3*transferfunc)/(1+DC3*transferfunc);
zpk(CLT);
numCLT=13.89*(s+25)*(s^2+0.81)*(s+0.5);
denCLT=(s+42.16)*(s+6.691)*(s^2+0.9674*s+0.6065)*(s^2+0.2234*s+0.822);
CLT=numCLT/denCLT;
x0 = [0 0 0 0];
simul=lsim(CLT,u1,taim,x0);
[A_CL,B_CL,C_CL,D_CL] = tf2ss([13.89 354.2 184.9 286.9 140.6],[1 50.04 341.9 417.2 509.9 286.9 140.6]);
CLss=ss(A_CL,B_CL,C_CL,D_CL);
[Ad_CL,Bd_CL,Cd_CL,Dd_CL]=c2dm(A_CL,B_CL,C_CL,D_CL,0.1,'zoh');
CLDss=ss(Ad_CL,Bd_CL,Cd_CL,Dd_CL,0.1);
figure(8)
rlocus(DC3*transferfunc);
figure(9)
bode(DC3*transferfunc)
grid on


%% sys 2  wn= 2 rad/s
% state space matrices of sys2

sys2=tf([1/50 1],[1/4 0.02 1 0 0]);
figure(10)
rlocus(DC3*sys2);
figure(11)
bode(DC3*sys2)

%% sys 3  collocated
A3=[0 1 0 0;-k/j2 -b/j2 k/j2 b/j2;0 0 0 01;k/j1 b/j1 -k/j1 -b/j1];
B3=[0;0;0;1/j1];
C3=[0 0 1 0;1 0 0 0];
D3=zeros(2,1);
sys3=ss(A3,B3,C3,D3);
transferfunc3=tf(sys3);

%% closed loop   DC-2
DC_2=tf([0.03 0.001],[1]);
figure(12)
rlocus(DC_2);
figure(13)
bode(DC_2)
