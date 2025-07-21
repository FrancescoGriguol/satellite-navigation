%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code interpolates satellite positions from an SP3 file data
%
% Inputs:   pos - matrix with satellite coordinates from SP3 file [km]
%           time - time span vector [s]
%           interp_time - vector of time where compute the interpolation [hours]
%
% Outputs:  orbit = interpolated coordinates of satellite orbit [km]
%           t = time of interpolated data in seconds 
%
% Coded by Anese Giovanni
% CISAS "Giuseppe Colombo"
% University of Padua
% March 12, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [orbit] = interp_precise_orbits(pos,time,interp_time)

t_hr = time/3600; %hours

n = length(t_hr) - 1; %% cosi ottengo il grado del polinomio, ora ugaule a 10

% memalloc
orbit = zeros(length(interp_time),3);

for i = 1:3
    coeff = polyfit(t_hr,pos(1:length(t_hr),i),n);
    orbit(:,i) = polyval(coeff,interp_time);
    res = polyval(coeff,t_hr) - pos(1:length(t_hr),i);
end

end