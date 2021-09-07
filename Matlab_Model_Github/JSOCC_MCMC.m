function J = JSOCC_MCMC(Ca,v_socc,k_socc)

%Getting this formalism, and constants, from Paper 4: A Deterministic Model
%Predicts the Properties of Stoachstic Calcium Oscillattions in Airway
%smooth Muscle Cells

%v_socc = 1.6; 
%k_socc = 100; 

J = v_socc*(k_socc.^4./(k_socc.^4+Ca.^4));

end 