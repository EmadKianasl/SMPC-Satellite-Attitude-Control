%% main

clc;clear;
close all
%% newmark algorithm params

tf = 40; % final time step
dt = 0.01; % size of the time step
t = 0:dt:tf;
samples = size(t);

%% system parameters
j1=1;
j2=0.1;
k=0.091;
b=0.036;
J=[j1 0;0 j2];
K=[k -k;-k k];
B=[b -b;-b b];
%forces
f1 = 0; % instrument package control torque
f1=f1+zeros(1,samples(1,2));
f2 = 0; % main satellite control torque
f2=f2+zeros(1,samples(1,2));
%% control_newmark

kp = 0.25;
kd = 0.5;
ki=0;
[deplc, velc, acclc,t,f1c,f2c,error] = newmark_PID(2,J,K,B,kp,ki,kd,tf,dt); %controlled system
R=[f1;f2];
[depl, vel, accl,t] = newmark(2,R,J,K,B,tf,dt); %uncontrolled system


%% tf closed loop
%% figures

figure(1)
subplot(2,1,1);
stairs(t,deplc(1,:),'LineWidth',2) ;
hold on
stairs(t,depl(1,:),'LineWidth',2) ;
grid on
legend('theta1c','theta-1')
ylabel('theta1 (rad)') ;
subplot(2,1,2);
stairs(t,deplc(2,:),'LineWidth',2) ;
hold on
stairs(t,depl(2,:),'LineWidth',2) ;
grid on
legend('theta2c','theta-2')
xlabel('time (sec)') ; ylabel('theta2 (rad)') ;

figure(2)
subplot(2,1,1);
stairs(t,velc(1,:),'LineWidth',2) ;
hold on 
stairs(t,vel(1,:),'LineWidth',2) ;
%ylim([-0.4 0.4]);
grid on
legend("d-theta1c","d-theta1");
ylabel("d-theta1 (rad/sec)") ;
subplot(2,1,2);
stairs(t,velc(2,:),'LineWidth',2) ;
hold on 
stairs(t,vel(2,:),'LineWidth',2) ;
%ylim([-0.4 0.4]);
grid on
legend("d-theta2c","d-theta2");
xlabel('time (sec)') ; ylabel("d-theta2 (rad/sec)") ;

figure(3)
subplot(2,1,1);
plot(t,acclc(1,:),'LineWidth',2) ;
grid on
hold on
stairs(t,accl(1,:),'LineWidth',2) ;
legend('dd-theta1c','dd-theta1');
ylabel('dd-theta1 (rad/sec^2)') ;
subplot(2,1,2);
stairs(t,acclc(2,:),'LineWidth',2) ;
grid on
hold on 
stairs(t,accl(2,:),'LineWidth',2) ;
legend('dd-theta2c','dd-theta2');
xlabel('time (sec)') ; ylabel('dd-theta2 (rad/sec^2)') ;

figure(5)
subplot(2,1,1);
stairs(t,f1c,'LineWidth',2) ;
grid on
ylabel('f1 (N)') ;
subplot(2,1,2);
stairs(t,f2c,'LineWidth',2) ;
grid on
xlabel('time (sec)') ; ylabel("f2 (N)") ;