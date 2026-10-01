clc;clear;
close all

load ('results')
samples2=0:5999;
%%FIGURES sigma=6
stepinfo_sys=stepinfo(Yc,samples)
stepinfo_syss=stepinfo(Ycs,samples)
stepinfo_sysn=stepinfo(Ycn,samples)
stepinfo_syssn=stepinfo(Ycsn,samples)
%% controlled sys
ax1=figure();
ax1.Color='white';

stairs(samples,Yc,'red','Linewidth',1.5);
hold on
stairs(samples,Ycs,'--','Color','blue','Linewidth',1.5);
xlim([0 30]);
ylabel('Angular deviation (rad)','Fontsize',14,'interpreter', 'Latex')
xlabel('Time (sec)','Fontsize',14,'interpreter', 'Latex')
legend({'CL','CL SMPC-COM'},'Fontsize',16,'FontName','Times New Roman')
grid on

%%savefig('disp1.fig')


%% control signal
ax2=figure();
ax2.Color='white';

subplot(2,1,1);
stairs(samples,U,'red','Linewidth',1.5);
hold on
stairs(samples,Us_rec,'--','Color','blue','Linewidth',1.5);
xlim([0 30]);
ylabel('Control signal (N.m)','Fontsize',14,'interpreter', 'Latex')
xlabel('Iteration','Fontsize',14,'interpreter', 'Latex')
legend({'CL','CL SMPC-COM'},'Fontsize',16,'FontName','Times New Roman')
grid on

subplot(2,1,2);
stairs(samples,Un,'red','Linewidth',1.5);
hold on
stairs(samples,Usn_rec,'--','Color','blue','Linewidth',1.5);
xlim([0 30]);
ylabel('Control signal (N.m)','Fontsize',14,'interpreter', 'Latex')
xlabel('Time (sec)','Fontsize',14,'interpreter', 'Latex')
legend({'CL (noisy signal)','CL SMPC-COM (noisy signal)'},'Fontsize',16,'FontName','Times New Roman')
grid on
%%savefig('disp1.fig')

%% controlled sys (noisy)
ax3=figure();
ax3.Color='white';

stairs(samples,Ycn,'red','Linewidth',1.5);
hold on
stairs(samples,Ycsn,'--','Color','blue','Linewidth',1.5);
xlim([0 30]);
ylabel('Angular deviation (rad)','Fontsize',14,'interpreter', 'Latex')
xlabel('Time (sec)','Fontsize',14,'interpreter', 'Latex')
legend({'CL (noisy signal)','CL SMPC-COM (noisy signal)'},'Fontsize',16,'FontName','Times New Roman')
grid on

%%savefig('disp1.fig')

%% server observation (control signal)
ax6=figure();
ax6.Color='white';
subplot(2,1,1);
for i=1:length(Us(1,:))
    scatter(samples2,Us(:,i),'Linewidth',1.5)
    hold on
end
ylabel('Exchanged ciphertext messages','Fontsize',14,'interpreter', 'Latex')
xlabel('Iteration','Fontsize',14,'interpreter', 'Latex')
    
leg={};
for i=1:length(Us(1,:))
    leg=[leg,strcat('Server',{' '},num2str(i))];
end
legend(leg,'Fontsize',16,'FontName','Times New Roman')
%savefig('Distributed Control Signal.fig')
subplot(2,1,2);
for i=1:length(Us(1,:))
    stairs(samples2,Us(:,i),'Linewidth',1.5)
    hold on
end
xlabel('Iteration','Fontsize',14,'interpreter', 'Latex')