
function [depl, vel, accl,t,f1c,f2c,error] = newmark_PID(sdof,m,k,c,kp,ki,kd,tf,dt)

% INPUT :                                                                %
%           sdof : System degree's of freedom                            %
%           [R]  : External applied load                                 %
%           [m]  : Assembeled Mass Matrix                                %
%           [k]  : Assembeled Stiffness MAtrix                           %
%           [c]  : Damping Matrix                                        %
%                                                                        %
% OUTPUT :                                                               %
%           depl : Displacement Response                                 %
%           vel  : Velocity                                              %
%           accl : Acceleration                                          %
%                                                                        %
%========================================================================%
clc
ti = 0 ;     % initial time step
t = ti:dt:tf;
nt = fix((tf-ti)/dt);      % number of time steps
% initializing the displacement, velocity and acceleration matrices

%initialization
depl = zeros(sdof,nt) ;
vel = zeros(sdof,nt) ;
accl = zeros(sdof,nt);
Reff = zeros(sdof,nt) ;
%load('el.mat');
%sizet=size(t);
%sizeel=size(elcen(2,:));
%ddxg = elcen(2,:);
%ddxg=[ddxg,zeros(1,sizet(1,2)-sizeel(1,2))];

%ddxg(t>6) = 0;
%f1c0 = zeros(1,nt+1) ;
%f2c0 = zeros(1,nt+1) ;
f1c = zeros(1,nt+1) ;
f2c = zeros(1,nt+1) ;
accl(:,1) = [0;0];
depl(:,1) = [0;0.2];
vel(:,1) = [0;0];
%f1c(:,1) = 0;
%f2c(:,1) = 0;

%external foreces
%R=[f1c- m(1,1)*ddxg;f2c- m(2,2)*ddxg];
% Solve for initial accelerations
%accl(:,1) = inv(m)*(R(:,1)-c*vel(:,1)-k*depl(:,1)); 
% Parameters for Newmark time integration
alpha = 0.25 ;delta = 0.5 ;
% Calculating integration constants
a0 = 1/(alpha*dt^2); a1 = delta/(alpha*dt) ; a2 = 1/(alpha*dt) ;
a3 = (1/(2*alpha))-1 ; a4 = (delta/alpha)-1 ;a5 = (dt/2)*(delta/alpha-2) ;
a6 = dt*(1-delta) ; a7 = delta*dt ;
% calculating effectvie stiffness matrix
keff = k+a0*m+a1*c ;

r = 1; %referrence input
e = 0; %
eprev = 0;
integral_prev = 0;
error = [];
% time step starts

for it = 1:nt
    
    e = r-depl(1,it);
    error=[error;e];
    %esum = e+eprev;
   % integeral = integral_prev + e*dt;
 %   f1c_past(it) = kp*e + ki*integeral + kd*(e-eprev)/dt;
 %   enext = (((f1c_past(it)-kp*e - ki*integeral)*2*dt)/kd)+eprev;
 %   f1c(it) = kp*e + ki*integeral + kd*(enext-eprev)/(2*dt);
    if size(error)> 2
        f1c(it) = kp*e + kd*(error(it-2)-4*eprev+3*e)/(2*dt);%+ ki*integeral ;
    else
        f1c(it) = kp*e + kd*(e-eprev)/dt;%+ ki*integeral ;
    end
    
%    if size(error)> 2
 %       f1c(it) = 1.558*f1c(it-1)-0.6065*f1c(it-2)+604.1*f1c0(it)-1208*f1c0(it-1)+604.1*f1c0(it-2); 
 %   else
%        f1c(it) = 604.1*f1c0(it);
%    end
    
    R = [f1c;f2c];
        % Solve for initial accelerations
    if it == 1
        accl(:,1) = inv(m)*(R(:,1)-c*vel(:,1)-k*depl(:,1)); 
    end
    
    Reff(:,it) = R(:,it)+m*(a0*depl(:,it)+a2*vel(:,it)+a3*accl(:,it)).....
                          +c*(a1*depl(:,it)+a4*vel(:,it)+a5*accl(:,it));
% solving for displacements at time (it+dt)
    depl(:,it+1)= keff\Reff(:,it);
             
% calculating velocities and accelerations at time (it+dt)
    accl(:,it+1) = a0*(depl(:,it+1)-depl(:,it))-a2*vel(:,it)....
             -a3*accl(:,it) ;
             vel(:,it+1) = vel(:,it)+a6*accl(:,it)+a7*accl(:,it+1);
    eprev = e;
  %  integral_prev = integeral;
 
end

%f1c(nt+1)= f1c(nt);