function J = JP2X7_new_copy_Potassium(X7O1,X7O2,X7O3,X7O4,V, Ki, V_osteo)

V_osteo = V_osteo*10^-15; 
%V_osteo = 6.5*10^-12; % osteoblast volume Originally 2.4*10^-12. But the volume measurment, based on her notes, should be 6.5*10^-12
z = 1; % Potassium valence
F = 96485; % Faraday's constant

K_flux_X7 = 0.477; % proportion of P2X7 flux that is Potassium (0.2 is an initial guess)  0.477

g_X7 = 7.5*10^-9; 

G_X7 = g_X7*(X7O1+X7O2+X7O3+X7O4); 

Vss_X7 = 0; % reversal potential


I_X7 = G_X7.*(V-Vss_X7); % current
%disp(Vss_X7)

J = -(10^6)*K_flux_X7*I_X7/(z*F*V_osteo); % flux


end