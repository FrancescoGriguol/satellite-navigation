function [x,y,z] = GLONASS_coordinates(ephem, interp_time)
% ephem is a table with PoistionX, VelocityX, AccelerationX, PoistionY, 
% VelocityY, AccelerationY, PoistionZ, VelocityZ, AccelerationZ

% Initial data
pos_0 = [ephem.PositionX,ephem.PositionY,ephem.PositionZ]';
vel_0 = [ephem.VelocityX,ephem.VelocityY,ephem.VelocityZ]';
y0 = [pos_0(1), vel_0(1), pos_0(2), vel_0(2), pos_0(3), vel_0(3)]';

% perturbation vector (luni solare)
acc_p = [ephem.AccelerationX,ephem.AccelerationY,ephem.AccelerationZ]';

% richiamiamo la ode
[t,y_out] = ode45(@(t,y)GLONASS_eq(y,acc_p),interp_time,y0);

x = y_out(:,1);
y = y_out(:,3);
z = y_out(:,5);

