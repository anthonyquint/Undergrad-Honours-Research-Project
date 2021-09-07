clc
clear all

IC = [0.2, 0, 0.7, 0.05, 0, 0.7, 200 1 0 0 0 0 0 0 0 0 0 0 0 0.1 60000 300000 22750 118750 6500]; %adding the P2x7 model initial conditions too 
tspan = 0:0.1:2000;

[t,S] = ode15s(@(t,S) shell_bulk_dynamic_IP3(t,S,0), tspan, IC);


IC1 = S(end,:); 
IC1(25) = 6500; 

tspan1 = 0:0.1:200;


%Cacyt = @(Cas,Cab) Cab + (1/48)*Cas;
%what is this^

Cacyt = @(Cas,Cab, V_cell) Cab + (1./(V_cell - (0.1.*V_cell))./((V_cell - ((4/3).*(pi).*((((V_cell.*(3./4).*(1./pi)).^(1./3)) - (0.08))).^3)) - (((V_cell - ((4./3).*(pi).*((((V_cell.*(3./4).*(1./pi)).^(1./3)) - (0.08))).^3))./V_cell).*((0.1.*V_cell))))).*Cas;

%% ATP Application: constant ctot

[t,S0] = ode15s(@(t,S) shell_bulk_dynamic_IP3(t,S,0.5), tspan1, IC1);
[t,S1] = ode15s(@(t,S) shell_bulk_dynamic_IP3(t,S,1), tspan1, IC1);
[t,S2] = ode15s(@(t,S) shell_bulk_dynamic_IP3(t,S,10), tspan1, IC1);
[t,S3] = ode15s(@(t,S) shell_bulk_dynamic_IP3(t,S,100), tspan1, IC1);
[t,S4] = ode15s(@(t,S) shell_bulk_dynamic_IP3(t,S,1000), tspan1, IC1);
[t,S5] = ode15s(@(t,S) shell_bulk_dynamic_IP3(t,S,10000), tspan1, IC1);
[t,S6] = ode15s(@(t,S) shell_bulk_dynamic_IP3(t,S,0.1), tspan1, IC1);

%%

% 
figure(20) 
plot(t, S4(:,25))
ylabel("Volume (um^3)", 'Fontsize', 14) 
xlabel("Time (s)", 'Fontsize', 14) 
title('Volume vs. time', 'Fontsize', 14) 


%%
JIPR_b = JIPR(S1(:,6),S1(:,4),S1(:,7),S1(:,5));
JL_b = JERLEAK(S1(:,4),S1(:,7));
JD = dif(S1(:,1),S1(:,4),S1(:,2),S1(:,5));
JS_b = JSERCA(S1(:,4));
%What's this for^?

% figure
% subplot(1,2,1)
% plot(t,JIPR_b)
% hold on
% plot(t,JL_b)
% plot(t,JD)
% plot(t,JIPR_b+JL_b+JD);
% plot(t,JS_b);
% hold off
% legend('IPR','Leak','Diff','In','SERCA');
% 
% subplot(1,2,2)
% plot(t,S1(:,4));


% figure(3)
% hold on 
% plot(t,S4(:,2));
% plot(t,S4(:,5));
% hold off
% legend({'Shell IP3','Bulk IP3'});

%%
figure(1)
subplot(1,3,1)
plot(t,S0(:,1),'linewidth',1);
hold on
plot(t,S1(:,1),'linewidth',1);
plot(t,S2(:,1),'linewidth',1);
plot(t,S3(:,1),'linewidth',1);
plot(t,S4(:,1),'linewidth',1);
plot(t,S5(:,1),'linewidth',1);
hold off
legend({'.5','1','10','100','1000','10000'});
xlabel('Time(s)');
ylabel('Calcium (\muM)');
set(gca,'FontSize',14);
title('Shell')

subplot(1,3,2)
plot(t,S0(:,4),'linewidth',1);
hold on
plot(t,S1(:,4),'linewidth',1);
plot(t,S2(:,4),'linewidth',1);
plot(t,S3(:,4),'linewidth',1);
plot(t,S4(:,4),'linewidth',1);
plot(t,S5(:,4),'linewidth',1);
hold off
legend({'.5','1','10','100','1000','10000'});
xlabel('Time(s)');
ylabel('Calcium (\muM)');
set(gca,'FontSize',14);
title('Bulk')

% figure
% plot(t,S0(:,7),'linewidth',1);
% hold on
% plot(t,S1(:,7),'linewidth',1);
% plot(t,S2(:,7),'linewidth',1);
% plot(t,S3(:,7),'linewidth',1);
% plot(t,S4(:,7),'linewidth',1);
% plot(t,S5(:,7),'linewidth',1);
% hold off
% legend({'.5','1','10','100','1000','10000'});
% xlabel('Time(s)');
% ylabel('Cytosolic Calcium (\muM)');
% set(gca,'FontSize',14);
% title('ER')


subplot(1,3,3)
plot(t,Cacyt(S0(:,1),S0(:,4),S0(:,25)),'linewidth',1);
hold on
plot(t,Cacyt(S1(:,1),S1(:,4),S1(:,25)),'linewidth',1);
plot(t,Cacyt(S2(:,1),S2(:,4),S2(:,25)),'linewidth',1);
plot(t,Cacyt(S3(:,1),S3(:,4),S3(:,25)),'linewidth',1);
plot(t,Cacyt(S4(:,1),S4(:,4),S4(:,25)),'linewidth',1);
plot(t,Cacyt(S5(:,1),S5(:,4),S5(:,25)),'linewidth',1);
hold off
set(gca,'FontSize',14);
legend({'.5','1','10','100','1000','10000'});
xlabel('Time(s)');
ylabel('Calcium (\muM)');
title('Cytosol');


figure(2)
subplot(1,2,1)
plot(t,S0(:,2));
hold on
plot(t,S1(:,2));
plot(t,S2(:,2));
plot(t,S3(:,2));
plot(t,S4(:,2));
plot(t,S5(:,2));
hold off
title('IP3 Shell');

subplot(1,2,2)
plot(t,S0(:,5));
hold on
plot(t,S1(:,5));
plot(t,S2(:,5));
plot(t,S3(:,5));
plot(t,S4(:,5));
plot(t,S5(:,5));
hold off
title('IP3 Bulk');

%My work. Plotting delaneys plots for figure 1 on their own individual
%graphs.
% figure(3)
% plot(t,Cacyt(S0(:,1),S0(:,4)),'linewidth',1);
% hold off
% set(gca,'FontSize',14);
% legend({'.5'});
% xlabel('Time(s)');
% ylabel('Calcium (\muM)');
% title('Cytosol');
% ylim([0 1.6])

% figure(100)
% subplot(2, 5, 2)
% plot(t,Cacyt(S1(:,1),S1(:,4)),'linewidth',1);
% hold off
% set(gca,'FontSize',14);
% legend({'1'});
% xlabel('Time(s)');
% ylabel('Calcium (\muM)');
% title('(b)');
% ylim([0 1.6])
% xlim([0 114])
% 
% subplot(2, 5, 3)
% plot(t,Cacyt(S2(:,1),S2(:,4)),'linewidth',1);
% hold off
% set(gca,'FontSize',14);
% legend({'10'});
% xlabel('Time(s)');
% ylabel('Calcium (\muM)');
% title('(c)');
% ylim([0 1.6])
% xlim([0 114])
% 
% subplot(2, 5, 4)
% plot(t,Cacyt(S3(:,1),S3(:,4)),'linewidth',1);
% hold off
% set(gca,'FontSize',14);
% legend({'100'});
% xlabel('Time(s)');
% ylabel('Calcium (\muM)');
% title('(d)');
% ylim([0 1.6])
% xlim([0 114])
% 
% subplot(2, 5, 5)
% plot(t,Cacyt(S4(:,1),S4(:,4)),'linewidth',1);
% hold off
% set(gca,'FontSize',14);
% legend({'1000'});
% xlabel('Time(s)');
% ylabel('Calcium (\muM)');
% title('(e)');
% ylim([0 1.6])
% xlim([0 114])

% figure(8)
% plot(t,Cacyt(S5(:,1),S5(:,4)),'linewidth',1);
% hold off
% set(gca,'FontSize',14);
% legend({'10000'});
% xlabel('Time(s)');
% ylabel('Calcium (\muM)');
% title('Cytosol');
% ylim([0 1.6])

% subplot(2, 5, 1)
% plot(t,Cacyt(S6(:,1),S6(:,4)),'linewidth',1);
% hold off
% set(gca,'FontSize',14);
% legend({'0.1'});
% xlabel('Time(s)');
% ylabel('Calcium (\muM)');
% title('(a)');
% ylim([0 1.6])
% xlim([0 114])

%plotting the cystolic graphs altogether on one plot, in one figure
figure(10)
plot(t,Cacyt(S0(:,1),S0(:,4),S0(:,25)),'linewidth',1);
hold on
plot(t,Cacyt(S1(:,1),S1(:,4),S1(:,25)),'linewidth',1);
plot(t,Cacyt(S2(:,1),S2(:,4),S2(:,25)),'linewidth',1);
plot(t,Cacyt(S3(:,1),S3(:,4),S3(:,25)),'linewidth',1);
plot(t,Cacyt(S4(:,1),S4(:,4),S4(:,25)),'linewidth',1);
plot(t,Cacyt(S5(:,1),S5(:,4),S5(:,25)),'linewidth',1);
hold off
set(gca,'FontSize',14);
legend({'.5','1','10','100','1000','10000'});
xlabel('Time(s)');
ylabel('Calcium (\muM)');
title('Cytosol');
ylim([0 1.6])

% figure(13) 
% plot(t,Cacyt(S4(:,1),S4(:,4)),'linewidth',1);
% legend("[ATP] = 1000um", 'Fontsize', 14) 
% xlabel('Time(s)', 'Fontsize', 14);
% ylabel('Calcium (\muM)', 'Fontsize', 14);
% title('Cytosolic Calcium', 'Fontsize', 14);

% figure(14) 
% 
% subplot(1,3,1) 
% plot(t,S4(:,1),'linewidth',1);
% legend({'1000'});
% xlabel('Time(s)');
% ylabel('Calcium (\muM)');
% set(gca,'FontSize',14);
% title('Shell')
% 
% subplot(1,3,2) 
% plot(t,S4(:,4),'linewidth',1);
% legend({'1000'});
% xlabel('Time(s)');
% ylabel('Calcium (\muM)');
% set(gca,'FontSize',14);
% title('Bulk')
% 
% subplot(1,3,3) 
% plot(t,Cacyt(S4(:,1),S4(:,4)),'linewidth',1);
% legend({'1000'});
% xlabel('Time(s)');
% ylabel('Calcium (\muM)');
% set(gca,'FontSize',14);
% title('Cytosol')

% figure(101)
% subplot(1, 2, 1)
% plot(t,Cacyt(S4(:,1),S4(:,4)),'linewidth',1);
% hold off
% set(gca,'FontSize',14);
% legend({'1000'});
% xlabel('Time(s)');
% ylabel('Calcium (\muM)');
% title('(m)');
% ylim([0 1.6])
% xlim([0 114])

% Figure of e-6 by itself 
figure(11)
plot(t,Cacyt(S2(:,1),S2(:,4),S2(:,25)),'linewidth',1);
xlabel('Time(s)');
ylabel('Calcium (\muM)');
title('Cytosol');
xlim([0 114]);







