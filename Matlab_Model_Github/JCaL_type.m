function J = JCaL_type(Cas, Cab, m, V, V_osteo)

%Commented this out since changing Cai to only be Cai = Cas
%Calculating Intracellular Calcium 
%Cai = Cab + (1/48)*Cas; %Originally as Cacyt = @(Cas,Cab) Cab + (1/48)*Cas, but whats the @Cas,Cab useful for? I took it out.
Cai = Cas; 



%parameters
z = 2; % calcium valence
F = 96485; % Faraday's constant
R = 8.314; %Does this need to be changed? 
T = 300; 
g = 3.5;   %Data say it is 3.5 pS (picosiemens), can I leave the value as is or do I need to change it? 
V_osteo = V_osteo*10^-15; % osteoblast volume
% V_osteo = 6.5*10^-12; % osteoblast volume
Cao = 2000; %average ECF calcium is 1-3 mM, but our data is in uM, so multiply by 1000 

%Other gating variable, already at steady state 
h = 0.00045/(0.00045 + Cai); 


%Nernst Potential for Calcium, from Chang's Paper
ECa = ((R*T)/(z*F))*(log(Cao/Cai)); 
 
%Current Equation from Zeng's Paper
I = g*m*h*(V - ECa);
%disp(h)



%Conversion of Current to Flux via Chang's Formalism 
J = -(10^(-6))*I/(z*F*V_osteo);
%(10^-4)*


