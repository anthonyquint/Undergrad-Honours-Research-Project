clc
clear all

IC = [0.5, 0.7, 0.01, 0.7, 200];
tspan = 0:0.01:2000;

[t,S] = ode15s(@(t,S) shell_bulk(t,S,0,0), tspan, IC);
    
Catot_s = @(Cas,Cer) (Cas + 0.02*Cer);
Catot_b = @(Cab,Cer) (Cab + 0.98*Cer);
Catot = @(Catot_s,Catot_b) Catot_s + Catot_b;

figure(1)
subplot(2,1,1)
plot(t,S(:,1));
hold on
plot(t,S(:,3));
hold off
legend('Shell Ca','Bulk Ca');
subplot(2,1,2)
plot(t,S(:,5));
title('ER Ca');

Ca_tot_s_rest = Catot_s(S(:,1),S(:,5));
Ca_tot_b_rest = Catot_b(S(:,3),S(:,5));

figure(2)
plot(t,Ca_tot_s_rest,'linewidth',2);
hold on
plot(t,Ca_tot_b_rest,'linewidth',2);
hold off
legend('Total Shell Ca Rest','Total Bulk Ca Rest');

% S starting at 0.0202 (outside of osc. regime for all IP3)
% B starting at 1 (within osc. regime for some IP3)

%% 

tspan1 = 0:0.001:200;
IC1 = S(end,:);

[t,S0] = ode15s(@(t,S) shell_bulk(t,S,10,10),tspan1,IC1);

Catot_s_stim = Catot_s(S0(:,1),S0(:,5));
Catot_b_stim = Catot_b(S0(:,3),S0(:,5));

figure(3)
subplot(1,4,1)
plot(t,S0(:,1),'linewidth',2);
title('Shell');
subplot(1,4,2)
plot(t,S0(:,3),'linewidth',2);
title('Bulk');
subplot(1,4,3)
plot(t,Catot_s_stim,'linewidth',2);
title('Total Shell');
subplot(1,4,4)
plot(t,Catot_b_stim,'linewidth',2);
title('Total Bulk');









