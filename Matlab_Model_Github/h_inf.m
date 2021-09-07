function h = h_inf(IP3,Ca)

Q = Q2(IP3);

h = Q./(Q+Ca);
end