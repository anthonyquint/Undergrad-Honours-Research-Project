function Q = Q2(IP3)


d1 = 0.13;
d2 = 1.049;
d3 = 0.9434;

Q = d2*((IP3+d1)./(IP3+d3));

% k = 0.04;
% a = 59;
% n = 7.941;
% 
% Q = a*IP3.^n./(IP3.^n + k^n);

end