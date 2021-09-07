clc
clear all

IC = [0.01, 0.7, 0];
tspan = 0:0.1:2000;

[t,S] = ode15s(@(t,S) li_rinzel_dynamic_IP3(t,S,0,1.5), tspan, IC);

IC1 = S(end,:);
tspan1 = 0:0.1:200;

%% ATP Application: constant ctot

[t,S0] = ode15s(@(t,S) li_rinzel_dynamic_IP3(t,S,0,1.5), tspan1, IC1);
[t,S1] = ode15s(@(t,S) li_rinzel_dynamic_IP3(t,S,0.6,1.5), tspan1, IC1);
[t,S2] = ode15s(@(t,S) li_rinzel_dynamic_IP3(t,S,0.8,1.5), tspan1, IC1);
[t,S3] = ode15s(@(t,S) li_rinzel_dynamic_IP3(t,S,1,1.5), tspan1, IC1);
[t,S4] = ode15s(@(t,S) li_rinzel_dynamic_IP3(t,S,2,1.5), tspan1, IC1);
[t,S5] = ode15s(@(t,S) li_rinzel_dynamic_IP3(t,S,10,1.5), tspan1, IC1);
[t,S6] = ode15s(@(t,S) li_rinzel_dynamic_IP3(t,S,100,1.5), tspan1, IC1);

figure
plot(t,S6(:,1),'linewidth',2);
xlabel('Time(s)');
ylabel('Cytosolic Calcium(\muM)');
set(gca,'fontsize',14);
title('ATP = 100');

figure(1)
plot(t,S0(:,1),'linewidth',1);
hold on
plot(t,S1(:,1),'linewidth',1);
plot(t,S2(:,1),'linewidth',1);
plot(t,S3(:,1),'linewidth',1);
plot(t,S4(:,1),'linewidth',1);
plot(t,S5(:,1),'linewidth',1);
plot(t,S6(:,1),'linewidth',1);
hold off
legend('0','0.6','0.8','1','2','10','100');
xlabel('Time(s)');
ylabel('Cytosolic Calcium (\muM)');
set(gca,'FontSize',14);