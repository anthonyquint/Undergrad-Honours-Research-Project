function dS = li_rinzel_dynamic_IP3(t,S,ATP,ctot)

c1 = 0.185;
d1 = 0.13;
d2 = 1.049;
d3 = 0.9434;
d5 = 0.08234;
a2 = 0.2;
v1 = 6;
v2 = 0.11;
v3 = 0.9;
k3 = 0.1;
alpha_IP3 = 0.133;
k_deg = 0.033;
k_IP3 = 2.9;
n = 1.4;

cer = (ctot - S(1))/c1;
m_inf = (S(3)/(S(3)+d1))*(S(1)/(S(1)+d5));
Q2 = d2*(S(3)+d1)/(S(3)+d3);
Tau_h = 1/(a2*(Q2+S(1)));
h_inf = Q2/(Q2+S(1));

dS = zeros(size(S));

dS(1) = -c1*v1*m_inf^3*S(2)^3*(S(1)-cer) - c1*v2*(S(1)-cer) - v3*S(1)^2/(k3^2+S(1)^2);
dS(2) = (h_inf-S(2))/Tau_h;
dS(3) = alpha_IP3*ATP^n/(ATP^n+k_IP3^n) - k_deg*S(3);

end