function dydt = GLONASS_eq(y,acc)

% mem alloc
dydt = zeros(6,1);

% Parameters
mu = 3.9860044e5;
J2 = -1.08263e-3;
w_e = 0.7292115e-4;
r_e = 6.378136e3;

r = sqrt(y(1)^2+y(3)^2+y(5)^2);

dydt(1) = y(2);
dydt(2) = -mu/r^3*y(1)+3/2*J2*(mu*r_e^2)/r^5*y(1)*(1-5*y(5)^2/r^2)+w_e^2*y(1)+2*w_e*y(4)+acc(1);
dydt(3) = y(4);
dydt(4) = -mu/r^3*y(3)+3/2*J2*(mu*r_e^2)/r^5*y(3)*(1-5*y(5)^2/r^2)+w_e^2*y(3)-2*w_e*y(2)+acc(2);
dydt(5) = y(6);
dydt(6) = -mu*y(5)/r^3+3/2*J2*mu*r_e^2/r^5*y(5)*(3-5*y(5)^2/r^2)+acc(3);


end