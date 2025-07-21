function [Xk,Yk,Zk]=GPS_coordinates(t, crs, delta_n, M0, cuc, e, cus, sqrt_a, toe, cic, nodo, cis, i0, crc, long_perigeo, nodo_dot, i_dot)

% GPS usa modello WGS-84

mu = 398600500000000; 
omega_e = 7.2921151467E-05;

% Computational time
tk = t-toe;

% Mean anomaly
Mk = M0+(sqrt(mu)/sqrt_a^3+delta_n)*tk;

% Eccentric anomaly
Ek0 = Mk+e*sin(Mk);
options = optimoptions("fsolve",'Display','none');
Ek = fsolve(@(x)(x-e*sin(x)-Mk),Ek0,options);

% True anomaly
vk = wrapTo2Pi(atan2((sqrt(1-e^2)*sin(Ek)),(cos(Ek)-e)));

% Latitude argument
uk = long_perigeo+vk+cuc*cos(2*(long_perigeo+vk))+cus*sin(2*(long_perigeo+vk));

% Radial distance
rk = sqrt_a^2*(1-e*cos(Ek))+crc*cos(2*(long_perigeo+vk))+crs*sin(2*(long_perigeo+vk));

% Inclination
ik = i0+i_dot*tk+cic*cos(2*(long_perigeo+vk))+cis*sin(2*(long_perigeo+vk));

% Longitude of ascending node
nodok = nodo+(nodo_dot-omega_e)*tk-omega_e*toe;

R3_nodok = [cos(-nodok) sin(-nodok) 0; -sin(-nodok) cos(-nodok) 0; 0 0 1];
R1_ik = [1 0 0; 0 cos(-ik) sin(-ik); 0 -sin(-ik) cos(-ik)];
R3_uk = [cos(-uk) sin(-uk) 0; -sin(-uk) cos(-uk) 0; 0 0 1];

CTS = R3_nodok*R1_ik*R3_uk*[rk 0 0]';

Xk = CTS(1);
Yk = CTS(2);
Zk = CTS(3);

end
