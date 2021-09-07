function J = JSERCA(Ca)

v_SERCA = 30; %originally 22.5. (30)
%v_SERCA = 1;

k_SERCA = 0.105; %originally 0.105

J = v_SERCA*(Ca.^2./(Ca.^2+k_SERCA.^2));

end