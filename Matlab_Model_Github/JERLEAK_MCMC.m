function J = JERLEAK_MCMC(Ca,Cer,v_ERleak)

%v_ERleak = 0.03; %0.03 originally

J = v_ERleak*(Cer-Ca);


end