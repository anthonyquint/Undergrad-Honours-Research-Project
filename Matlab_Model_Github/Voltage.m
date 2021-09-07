function V = Voltage(Cas, Cab, Nai, Ki)

%Commented this out since changing Cai to only be Cai = Cas
%Cai = Cab + (1/48)*Cas; %Originally as Cacyt = @(Cas,Cab) Cab + (1/48)*Cas, but whats the @Cas,Cab useful for? I took it out. 
Cai = Cas; 

%Paramters
F = 96485;      % Faraday's constant
R = 8.314;      %Does this need to be changed? 
T = 300; 
%Ki = 125000;    %used paper found June 30, ICF K is 125mM, therefore 125 000 uM 
Ko = 4000;      %used google 
Nao = 140000;   %used google 
%Nai = 24000;    %used paper found June 30, ICF Na is 24mM, therefore 24 000 uM   
Cao = 2000;     %average ECF calcium is 1-3 mM, but our data is in uM, so multiply by 1000

%Permeability coefficients 
PCa = 35;      %Figure this out later 
yK = 1;       %guessing 
JK = 1;       %guessing  
PDK = 1;      %guessing 
PK = 30;       %for PK and PNa, since chang paper has ratio of PNa/PK ot be 0.01, tried to recreate this 
PNa = 1; 

%Voltage formalism from Chang's Paper
V = (R*T/F)*log((PK*Ko + PNa*Nao + 4*PCa*Cao + (yK*JK*PK/PDK))/(PK*Ki + PNa*Nai + 4*PCa*Cai)); 
%V = -0.052; %Gives basically same result 
%disp(V)
%0.04





