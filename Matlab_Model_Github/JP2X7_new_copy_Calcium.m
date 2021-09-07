function J = JP2X7_new_copy_Calcium(X7O1,X7O2,X7O3,X7O4,V, Cas, V_osteo)

V_osteo = V_osteo*10^-15; % osteoblast volume Originally 2.4*10^-12. But the volume measurment, based on her notes, should be 6.5*10^-12
%V_osteo = 6.5*10^-12; 
z = 2; % calcium valence
F = 96485; % Faraday's constant

Ca_flux_X7 = 0.046; % proportion of P2X7 flux that is calcium 0.046

g_X7 = 7.5*10^-9; 

G_X7 = g_X7*(X7O1+X7O2+X7O3+X7O4); 
%disp(G_X7)

Vss_X7 = 0; % reversal potential


I_X7 = G_X7.*(V-Vss_X7); % current

J = -(10^6)*Ca_flux_X7*I_X7/(z*F*V_osteo); % flux


end