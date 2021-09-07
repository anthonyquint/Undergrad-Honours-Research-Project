clc
clear all

IC = [0.01, 0.7];
tspan = 0:0.1:2000;

[t,S] = ode15s(@(t,S) li_rinzel(t,S,1.5,0), tspan, IC);
disp(S)

IC1 = S(end,:);
tspan1 = 0:0.1:200;

%% Relationship between Dip and Peak Size for increasing IP3

p = 4:0.5:50;
dip = zeros(size(p));
peak = zeros(size(p));
for i = 1:length(p) 
    ip3 = p(i);
    [t,S6] = ode15s(@(t,S) li_rinzel(t,S,1.5,ip3), tspan1, IC1);
    y = S6(:,1);
    peak(i) = max(y);
    y_ss = y(end);
    index = find(y == max(y));
    y = y(index:end);
    dip(i) = y_ss-min(y);
    i = i+1;
end

figure
scatter(peak,dip,'filled');
hold on
f = polyfit(peak,dip,1);
p_plot = f(1)*peak + f(2);
plot(peak,p_plot,'linewidth',2,'color','r');
hold off

% the larger the peak, the smaller the dip
% approximate linear relationship

figure
subplot(2,1,1)
scatter(p,dip,'filled')
subplot(2,1,2)
scatter(p,peak,'filled')


% dip saturates because peak saturates 
    % at least for c_tot = 1.5

%% What about increases c_tot?


c = 0.1:0.1:10;
dip_c = zeros(size(c));
peak_c = zeros(size(c));
for i = 1:length(c) 
    ca = c(i);
    [t,S6] = ode15s(@(t,S) li_rinzel(t,S,ca,4), tspan1, IC1);
    y = S6(:,1);
    peak_c(i) = max(y);
    y_ss = y(end);
    index = find(y == max(y));
    y = y(index:end);
    dip_c(i) = y_ss-min(y);
    i = i+1;
end

figure
scatter(peak_c,dip_c,'filled');



