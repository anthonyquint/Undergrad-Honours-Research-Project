function J = JPMCA_MCMC(Ca,v_PMCA,k_PMCA)

%v_PMCA = 8; %8
%k_PMCA = 0.1; %0.1

J = v_PMCA*(Ca.^2./(Ca.^2+k_PMCA.^2));

end