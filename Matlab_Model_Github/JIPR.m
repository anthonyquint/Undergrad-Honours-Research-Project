function J = JIPR(h,Ca,Cer,IP3)

% Parameters
v_IPR = 15;

%v_IPR = 25;

% Functions
m = m_inf(IP3);
n = n_inf(Ca);
O = O_IPR(m,n,h);

J = v_IPR.*O.*(Cer - Ca);


end

