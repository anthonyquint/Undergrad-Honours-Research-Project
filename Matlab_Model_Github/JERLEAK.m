function J = JERLEAK(Ca,Cer)

v_ERleak = 0.03; %0.03 originally

J = v_ERleak*(Cer-Ca);


end