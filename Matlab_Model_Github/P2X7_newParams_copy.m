function dS = P2X7_newParams_copy(t,S,ATP)

k = 1/31.625; % conversion constant from BzATP data to ATP data
A = k*ATP*10^-6;

%% Parameters

k1_X7 = 0.3;
k2_X7 = 40000;
k3_X7 = 2.4;
k4_X7 = 50000;
k5_X7 = 1.58;
k6_X7 = 7000;
L1_X7 = 0.0001;
L2_X7 = 0.004;
L3_X7 = 0.3;
H1_X7 = 0.001;
H2_X7 = 0.01; % H2_C2
H2 = 0.01; % H2_Q1
H2_ = 0.7; % H2_Q2
H3_X7 = 0;
k8_X7 = 10000;

%% Model

dS=zeros(size(S));

% Markov Model Middle Row
dS(1) = H1_X7*S(9) + k1_X7*S(2) + L1_X7*S(8) - 3*k2_X7*A*S(1);%C1
dS(2) = H3_X7*S(10)+3*k2_X7*A*S(1) + 2*k3_X7*S(3) - (k1_X7 + 2*k4_X7*A +H2_X7)*S(2); %C2
dS(3) = 2*k4_X7*A*S(2) + 3*k5_X7*S(4) - (2*k3_X7 + k6_X7*A + H2)*S(3); %Q1
dS(4) = k6_X7*A*S(3) + L2_X7*S(5) - (3*k5_X7 + L3_X7 + H2_)*S(4);%Q2

% Markov Model Bottom Row
dS(5) = k8_X7*A*S(6) + L3_X7*S(4) - (3*k1_X7 + L2_X7)*S(5);%Q3
dS(6) = 2*k8_X7*A*S(7) + 3*k1_X7*S(5) - (2*k1_X7 +k8_X7*A)*S(6);%Q4
dS(7) = 3*k8_X7*A*S(8) + 2*k1_X7*S(6) - (k1_X7 + 2*k8_X7*A)*S(7);%C3
dS(8) = k1_X7*S(7) - (L1_X7 + 3*k8_X7*A)*S(8);%C4

% Markov Model Top Row (Desensitization)
dS(9) = k1_X7*S(10) - (3*k2_X7*A + H1_X7)*S(9); %D1
dS(10) = 3*k2_X7*A*S(9) +2*k3_X7*S(11) + H2_X7*S(2) - (k1_X7 + 2*k4_X7*A + H3_X7)*S(10); %D2
dS(11) = 2*k4_X7*A*S(10) +3*k5_X7*S(12) + H2*S(3) - (2*k3_X7 +k6_X7*A)*S(11); %D3
dS(12) = k6_X7*A*S(11) +H2_*S(4) -3*k5_X7*S(12);

end

