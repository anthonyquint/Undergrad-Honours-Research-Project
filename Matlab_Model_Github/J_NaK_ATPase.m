function J = J_NaK_ATPase(Nai, Ki, Nae, Ke, V)

%Reaction forward (+) and reaction backward (_) step parameters 
k1_plus = 1050;
k2_plus = 481; 
k3_plus = 2000;
k4_plus = 320; 
k1_minus = 172.1;
k2_minus = 40; 
k3_minus = 79300;
k4_minus = 40; 

K0_dNae = 15.5;
K_dKe = 0.213;
K_dMgATP = 2.51; 
K0_dNai = 2.49;    %in mM
K_dKi = 0.500;     %in mM

y = -0.031;        %y stands for ∆ which is the symbol used in the paper

MgATP = 4.99; 
Pi = 4.95; 
pH = 7.09;
MgADP = 0.06; 
H = 1000*10^(-pH);

R = 8.314;             %Is this the correct R value?  
T = 300;           %Is this the correct value for temperature? 
F = 96485;              %Fill in this value correctly later, Faradays constant?

K_dNai = K0_dNai*exp((y*F*V)/(3*R*T));
Naii = Nai/K_dNai;

Kii = Ki/K_dKi;

K_dNae = K0_dNae*exp(((1 + y)*F*V)/(3*R*T));
Naee = Nae/K_dNae;

Kee = Ke/K_dKe;

MgATPP = MgATP/K_dMgATP;

%Reaction forward (+) and reaction backward (_) step equations

a1_plus = (k1_plus*Naii^3)/((1+ Naii)^3 + (1 + Kii)^2 - 1);
a2_plus = k2_plus;
a3_plus = (k3_plus*Kee^2)/((1+ Naee)^3 + (1 + Kee)^2 - 1);
a4_plus = (k4_plus*MgATPP)/(1 + MgATPP); 
a1_minus = k1_minus*MgADP; 
a2_minus = (k2_minus*Naee^3)/((1+ Naee)^3 + (1 + Kee)^2 - 1);
a3_minus = (k3_minus*Pi*H)/(1 + MgATPP);
a4_minus = (k4_minus*Kii^2)/((1+ Naii)^3 + (1 + Kii)^2 - 1);


%4 state model and the resutling steady state flux representation as in Smith and Crampin (2004)
x = (a1_minus*a2_minus*a3_minus) + (a1_minus*a2_minus*a4_plus) + (a1_minus*a3_plus*a4_plus) + (a2_plus*a3_plus*a4_plus) + (a2_minus*a3_minus*a4_minus) + (a1_plus*a2_minus*a3_minus) + (a1_plus*a2_minus*a4_plus) + (a1_plus*a3_plus*a4_plus) + (a1_minus*a3_minus*a4_minus) + (a2_plus*a3_minus*a4_minus) + (a1_plus*a2_plus*a3_minus) + (a1_plus*a2_plus*a4_plus) + (a1_minus*a2_minus*a4_minus) + (a1_minus*a3_plus*a4_minus) + (a2_plus*a3_plus*a4_minus) + (a1_plus*a2_plus*a3_plus);        

v_NaK = ((a1_plus*a2_plus*a3_plus*a4_plus) - (a1_minus*a2_minus*a3_minus*a4_minus))/x;
        

% Density of Na-K ATPase (putting 1 for now, replace when find real value)
a_NaK = 2200; %in um^-2. Found this in Crampin paper, couldnt find it in Sneyd's. 

%Flux Equation
J = a_NaK*v_NaK; 








