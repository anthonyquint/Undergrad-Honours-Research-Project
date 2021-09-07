function dS = shell_bulk_dynamic_IP3_MCMC(t,S,ATP, v_SERCA,k_SERCA,v_PMCA,k_PMCA,calcium_diffusion,IP3_diffusion,v_socc,k_socc,v_ERleak,v_IPR)


f = 0.01; % buffering
%% Diffusion between Compartments
[J_diff_c, J_diff_p] = dif_MCMC(S(1),S(4),S(2),S(5),calcium_diffusion,IP3_diffusion);

%% Volume Ratio Calculations 
% vs = 48; % cyt vol/ shell vol   48
% vb = 1; %cyt vol/ bulk vol
% vr = 9.4; % cyt vol/ ER vol

V_cell = S(25); % 6500um^3
l_shell = 0.08; % 0.08 um length of shell compartment, i.e. from exterior to centre

r_cell = (V_cell*(3/4)*(1/pi))^(1/3);
r_b = r_cell - l_shell; 
V_b = (4/3)*(pi)*(r_b)^3; 
V_Bcell = V_b/V_cell; 

V_s = V_cell - V_b; 
V_Scell = V_s/V_cell; 

V_er = 0.1*V_cell; 
V_cyt = V_cell - V_er; 

V_er_s = V_Scell*(V_er); 
V_cyt_s = V_s - V_er_s; 

V_er_b = V_Bcell*(V_er); 
V_cyt_b = V_b - V_er_b;


vs = V_cyt/V_cyt_s; 
vb = V_cyt/V_cyt_b; 
vr = V_cyt/V_er; 

%% Shell
% Calcium 
J_IPR_s = JIPR_MCMC(S(3),S(1),S(7),S(2),v_IPR);
J_leak_s = JERLEAK_MCMC(S(1),S(7),v_ERleak);
J_SERCA_s = JSERCA_MCMC(S(1),v_SERCA,k_SERCA);


%Adding L type calcium to code 
V = Voltage(S(1), S(4), S(21), S(22)); 
mL_inf_s = mL_inf(V); 
Tm_L_s = Tm_L(V); 
J_CaL_type_s = JCaL_type(S(1), S(4), S(20), V, V_cell); 

% flux through P2X7Rs, I added to Delaney's code
J_P2X7_s = JP2X7_new_copy_Calcium(S(10),S(11),S(12),S(13), V, S(1),V_cell); %Switching her constant V=-0.06 to my "V" to see what happens
%I added PMCA to Delaneys code
J_PMCA_s = JPMCA_MCMC(S(1),v_PMCA,k_PMCA); 
%I added Socc to Delaney's Code
J_SOCC_s = JSOCC_MCMC(S(7),v_socc,k_socc);  

% Adding Sodium Flux through P2X 
J_P2X7_Na_s = JP2X7_new_copy_Sodium(S(10),S(11),S(12),S(13), V, S(21),V_cell);
% Adding Potassium Flux through P2X
J_P2X7_K_s = JP2X7_new_copy_Potassium(S(10),S(11),S(12),S(13), V, S(22), V_cell); 
%Adding Sodium Potassium Pump 
J_NaK_ATPase_s = J_NaK_ATPase((10^(-3))*S(21), (10^(-3))*S(22), 140, 4, V); 
%Adding Sodium diffusion from Shell to Bulk 
J_diff_Na = diff_Na(S(21), S(23)); 
%Adding Potassium Diffusion from Shell to bulk 
J_diff_K = diff_K(S(22), S(24)); 


% IP3
a_ATP = 0.03; % IP3 production rate, 0.133
k_ATP = 1; % affinity of P2YR for ATP 2.9
n = 1; % hill coefficient
k_deg = 0.01; %IP3 degradation 0.033

% h (IP3R inactivation)
h_inf_s = h_inf(S(2),S(1));
tau_h_s = Tau_h(S(2),S(1));
%% Bulk

% Calcium
J_IPR_b = JIPR_MCMC(S(6),S(4),S(7),S(5),v_IPR);
J_leak_b = JERLEAK_MCMC(S(4),S(7),v_ERleak);
J_SERCA_b = JSERCA_MCMC(S(4),v_SERCA,k_SERCA);

% h (IP3R inactivation)
h_inf_b = h_inf(S(5),S(4));
tau_h_b = Tau_h(S(5),S(4));

%% ER
vr = 9.4; % cyt vol/ ER vol
fer = 0.025; %ER calcium 
%% Model

dS = zeros(size(S)); 

% Shell calcium
dS(1) = f*vs*(0.15 + J_leak_s - J_SERCA_s - J_diff_c - J_PMCA_s + (1)*J_P2X7_s + J_SOCC_s + J_CaL_type_s) + (1)*f*(vs)*J_IPR_s; %added JPMCA, JINLEAK = 0.15, JP2X7, Jsocc, can add J_CaL_type_s  

z = f*vs*((1)*J_P2X7_s - J_PMCA_s + 0.15 + J_SOCC_s); 

%disp(z) 

% Shell IP3
dS(2) = vs*(a_ATP*(ATP^n/(ATP^n + k_ATP^n)) - k_deg*S(2) - J_diff_p);
% Shell IP3R inactivation
dS(3) = (h_inf_s - S(3))/tau_h_s;

% Bulk calcium
dS(4) = f*vb*(J_leak_b - J_SERCA_b + J_IPR_b + J_diff_c); 
% Bulk IP3
dS(5) = vb*(J_diff_p - k_deg*S(5));
% Bulk IP3R inactivation
dS(6) = (h_inf_b - S(6))/tau_h_b;

% ER calcium
dS(7) = fer*vr*(J_SERCA_s + J_SERCA_b - J_IPR_s - J_IPR_b - J_leak_s - J_leak_b);

% P2X7 
dS(8:19) = P2X7_newParams_copy(t,S(8:19),ATP);

%%%%%%%%%% Adding mL_gating here, i.e. the "m" gating for L type calcium 

dS(20) = (mL_inf_s - S(20))/Tm_L_s;

% Shell Sodium Flux 
dS(21) = vs*((1)*J_P2X7_Na_s - 3*(1000)*J_NaK_ATPase_s - J_diff_Na); 

%Shell Potassium Flux
dS(22) = vs*(2*(1000)*J_NaK_ATPase_s - (1)*J_P2X7_K_s - J_diff_K);

%Bulk Sodium Flux 
dS(23) = vb*(J_diff_Na); 

%Bulk Potassium Flux 
dS(24) = vb*(J_diff_K); 

%Shrinking Volume 
dS(25) = 0; % -3  




end