function [J_diff_c, J_diff_p] = dif(Cs,Cb,Ps,Pb)

% vc = 2.87;
% vp = 1.52;

% vc = 28.24; % dividing by radius (small cell)
% vp = 14.92; 

vc = 12; %7
vp = 7; %7

J_diff_c = vc*(Cs-Cb);

J_diff_p = vp*(Ps-Pb);

end