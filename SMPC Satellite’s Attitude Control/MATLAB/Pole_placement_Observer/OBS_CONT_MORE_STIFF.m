clc;clear;
close all

%% system parameters
j1=1;
j2=0.1;
k=0.091;
b=0.0036;
tfinal = 30; % final time step
stime = 0.1; % size of the time step
taim = 0:stime:tfinal;
samps = size(taim);
u1 = 1; % instrument package control torque
u1=u1+zeros(1,samps(1,2));

%% system parameters (Modal)
jm1=1.1;
jm2=11;
km1=0;
km2=11.011;
bm1=0;
bm2=0.4402;


%% system matrices

J=[j1 0;0 j2];
K=[k -k;-k k];
B=[b -b;-b b];
Jm=[jm1 0;0 jm2];
Km=[km1 0;0 km2];
Bm=[bm1 0;0 bm2];
phi=[1 1;1 -10]; %modal matrix
phi2=[-10 0 1 0;0 0 0 0;1 0 1 0;0 0 0 0]; %extended modal matrix
%% state space matrices
[F,G,H,D]=tf2ss([1/50 1],[1/4 0.02 1 0 0]);
[A,B,C,V]=c2dm(F,G,H,D,0.005,'zoh'); %Discrete systems
Dsys=ss(A,B,C,V,0.005);
%F=[0 1 0 0;-k/j2 -b/j2 k/j2 b/j2;0 0 0 01;k/j1 b/j1 -k/j1 -b/j1];
%G=[0;0;0;1/j1];
%H=[1 0 0 0];
%D=zeros(1,1);
%transferfunc=tf([0.036 0.9],[1 0.04 1 0 0]);
sys=ss(F,G,H,D);
%[A,B,C,V]=c2dm(F,G,H,D,0.005,'zoh'); %Discrete systems
%Dsys=ss(A,B,C,V,0.005);



%% state space matrices (Modal)
Fm=[0 1 0 0;-km2/jm2 -bm2/jm2 0 0;0 0 0 1;0 0 -km1/jm1 -bm1/jm1];
Gm=[0;1/jm2;0;1/jm1];
Hm=[1 0 0 0]*phi2;
Dm=zeros(1,1);
sysm=ss(Fm,Gm,Hm,Dm);
figure(4)
lsim(sysm,ones(size(taim)),taim,'b'); 
hold on 
lsim(sys,ones(size(taim)),taim,'r');
%% Pole placement
p_cont=[-0.45+0.34j -0.45-0.34j -0.15+1.05j -0.15-1.05j]; % desired poles location corresponding to wbw=0.5 
%Kc = acker(A,B,p_cont)
Kc=place(F,G,p_cont);
%Kc=[-0.26139 0.041515 0.65538 1.1604];  %manual
%Kc=[-0.2788 0.0546 0.6814 1.1655];  %book
Kr = rscale(F,G,H,D,Kc);
%% observer
syms r l1 l2 l3 l4 
%chr=r*eye(4,4)-(A-[l1;l2;l3;l4]*[1 0 0 0]);
%det(chr);
%eq1=l1+0.0396-22.04;
%eq2=0.0396*l1+l2+1.001-243.9253;
%eq3=0.091*l1+0.0036*l2+0.91*l3+0.036*l4-1577.053476;
%eq4=0.091*l2+0.91*l4-5014.27063556;
%sol=solve(eq1,eq2,eq3,eq4);
%L_obs=vpa([sol.l1 sol.l2 sol.l3 sol.l4]',8); %manual
%Ke = [22.0004 242.0531 1512.84032 5485.9822]'; %manual
%computer algorithm
p_obs=[-7.7-3.12j -7.7+3.12j -3.32-7.85j -3.32+7.85j];
Ke=place(F',H',p_obs);
Ke=Ke';
%Book
%Ke=[22 242.3 1515.4 5503.9 ]'; %book

%% compensator state space (pole placement + estimator)
F_c=F-Ke*H;
G_c=[G Ke];
H_c=-Kc;%[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];%
D_c=[D D];%[0 0;0 0;0 0;0 0];%
Csys = ss(F_c,G_c,H_c,D_c);
Csys_tf=tf(Csys);
ss_C=ss(Csys_tf);
%Csys_tf=zpk(Csys_tf);
%num_C={[-1.16 -26.23 -297.4 -1722] [-7362 -3708 -6648 -1976] };
%den_C=[1 22.04 243.9 1577 5014];
%tf_C = tf(num_C,den_C);
%ss_C=ss(tf_C);
[A_c,B_c,C_c,V_c]=c2dm(ss_C.A,ss_C.B,ss_C.C,ss_C.D,0.005,'zoh'); %Discrete system
%[A_c,B_c,C_c,V_c]=c2dm(F_c,G_c,H_c,D_c,0.005,'zoh');
DCsys=ss(A_c,B_c,C_c,V_c,0.005);
%% compensator tf (pole placement + estimator)
s=tf('s');
%CT=-Kc*inv(s*eye(4,4)-F+G*Kc+Ke*H)*Ke; % 7.8 franklin book
%pretty(simplify(CT));
%CT=tf([-3727242495000 -1941038863492 -3386976266199 -1012659437790],20*[25000000 580127500 6763463090 46987780845 169194202391]); %book
CT = tf([7362 3708 6648 1976],[1 22.04 243.9 1577 5014]);
CT=zpk(CT);
figure(1)
bode(CT)


%% CLT
%with error system
%A_CL=[A-B*K_cont B*K_cont;zeros(size(A)) A-L_obs*C]; 
%B_CL=[B*Nbar;zeros(size(B))];
%C_CL=[C zeros(size(C))];
%with compensator system
F_CL=[F -G*Kc;Ke*H F-Ke*H-G*Kc];  %compensator system
G_CL=[G*Kr;G*Kr];
H_CL=[H zeros(size(H))];
sys_CL = ss(F_CL,G_CL,H_CL,0);  
figure(2)
x0=[0.0 0 0 0];
lsim(sys_CL,ones(size(taim)),taim,[x0 x0]);
CLT=tf(sys_CL);
CLT=minreal(CLT);
figure(3)
%rlocus(transferfunc,'r')
%hold on
%rlocus(CLT,'b')