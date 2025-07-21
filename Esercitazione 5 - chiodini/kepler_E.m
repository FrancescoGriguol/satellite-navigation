function E = kepler_E(e,M)

% fisso errore di tolleranza 
error = 1.0e-6;
%valore iniziale di E
if M < pi
    E = M+e/2;
else
    E = M-e/2;
end

ratio = 1;
while abs(ratio) > error
    ratio = (E-e*sin(E)-M)/(1-e*cos(E));
    E= E-ratio;
end