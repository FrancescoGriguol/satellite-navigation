function [RecPos,RecPos_lla,unc,DOP,COV_xyz] = GNSS_fix_ti(time_i,obs_i,navmsg,LeapSeconds,RecPos)

% Dept. of Industrial Engineering - University of Padova
% MS course: Satellite Navigation - GNSS
% author: A. Valmorbida, Ph.D
% date: March 27, 2024
% file name: GNSS_fix_ti.m

global wE c mu g
global GNSS_type
global RecPos_lla_approx RecPos_lla_ref
global Tropo Iono Clock Rel

maskAngle = 15; % per non considerare quelli troppo bassi per il multipath
numVar = 4; % tre dimensioni più lo shift del ricevitore

%% get pseudoranges

if GNSS_type == "GPS"
    psr1 = obs_i.C1C; % pseudorange relativo alla portante L1
    psr2 = obs_i.C2S; % pseudorange relativo alla portante L2
    f1 = 10.23*10^6*154; % Hz
    f2 = 10.23*10^6*120; % Hz
elseif GNSS_type == "Galileo"
    psr1 = obs_i.C1B;
    psr2 = obs_i.C7I;
    f1 = 10.23*10^6*154; % Hz
    f2 = 10.23*10^6*115; % Hz
end

%% IONOSPHERIC CORRECTION through measurements combination
% a1 + a2 = 1               => geometry-preserving
% a1 + f1^2/f2^2*a2 = 0     => ionosphere-free combination
a2 = -(f2^2)/(f1^2-f2^2); % => -1.5457
a1 = 1 - a2; % = f1^2/(f1^2-f2^2) => 2.5457

if Iono == 0
    % no ionospere correction
    psr = psr1;
elseif Iono == 1
    % ionosphere correction
    psr = a1*psr1 + a2*psr2;
end

% Gestione dei NaN nelle pseudorange
idx_nan = find(isnan(psr));
if ~isempty(idx_nan)
    psr(idx_nan) = psr1(idx_nan); % Usa psr1 se psr è NaN
end

% Controlla se ci sono ancora NaN
idx_nan = find(isnan(psr));
if ~isempty(idx_nan)
    % Elimina solo le righe rimanenti con NaN
    psr(idx_nan) = [];
    obs_i(idx_nan, :) = [];
    navmsg(idx_nan, :) = [];
end
numSat = size(obs_i,1); % numero dei satelliti in vista

% Controlla se il numero di satelliti è inferiore a 4
if numSat < 4
    RecPos = NaN;
    RecPos_lla = NaN;
    unc = NaN;
    DOP = NaN;
    return;
end


%% SATELLITE CLOCK CORRECTION
DeltaTime = seconds(obs_i.Time - navmsg.Time);
dpsr_svclock = c.*(navmsg.SVClockBias + navmsg.SVClockDrift.*DeltaTime + ...
               navmsg.SVClockDriftRate.*DeltaTime.^2 );    

%% RELATIVISTIC CLOCK CORRECTION

% Time of flight of the signal
time_flight = psr/c+dpsr_svclock/c;

Delta_rel_sat = zeros(numSat,1);
for i=1:numSat
    t_ems = time_i.Hour*3600+time_i.Minute*60+time_i.Second-LeapSeconds-time_flight(i);
    toe = mod(navmsg.Toe(i),(3600*24));
    Dt = t_ems-toe;
    Mk = navmsg.M0(i)+(sqrt(mu)/navmsg.sqrtA(i)^3+navmsg.Delta_n(i))*Dt;
    Ek0 = Mk+navmsg.Eccentricity(i)*sin(Mk);
    options = optimoptions("fsolve",'Display','none');
    Ek = fsolve(@(x)(x-navmsg.Eccentricity(i)*sin(x)-Mk),Ek0,options);
    Delta_rel_sat(i) = -2*(sqrt(mu)*navmsg.sqrtA(i)/c^2)*navmsg.Eccentricity(i)*sin(Ek);
end


%% COMPUTE SATELLITE POSITIONS USING BROADCAST EPHEMERIS

% memalloc
SatPos = NaN(numSat,3);
satID_GPS = NaN(numSat,1);
for i=1:numSat
    t = time_i;
    t.Second = t.Second-time_flight(i)-LeapSeconds;
    [sat_pos,~,satID_GPS(i,1)]=gnssconstellation(t,RINEXData=navmsg(i,:));
    SatPos(i,:) = sat_pos;
    theta = 1*wE*time_flight(i);
    SatPos(i,1) = sat_pos(1)*cos(theta)+sat_pos(2)*sin(theta);
    SatPos(i,2) = -sat_pos(1)*sin(theta)+sat_pos(2)*cos(theta);
    SatPos(i,3) = sat_pos(3);
end

%% TROPOSPHERIC CORRECTION (Saastamoinen)

% Trop_correction = zeros(numSat,1);
% 
% % Grandezze fisiche misurate al ricevitore
% T = 17; % [°C]
% P = 1.017; % [hPa]
% RH = 0.7; % relative humidity
% 
% % Altezza del ricevitore in metri
% h = RecPos_lla_approx(3);
% % Latitudine del ricevitore in radianti
% latitude = deg2rad(RecPos_lla_ref(1));
% 
% % Saturation vapor pressure using Tetens formula
% Es = 6.112 * exp((17.67 * T) / (T + 243.5)); % hPa
% % Water vapor pressure
% E = RH * Es; % hPa
% % Zenith Hydrostatic Delay
%     Kd = 0.002277 * P ./ (1 - 0.0026 * cos(2 * latitude) - 0.00028 * h / 1000);
% 
%     % Zenith Wet Delay
%     Kw = 0.002277 * (1255 / (T + 273.15) + 0.05) * E;
% 
%     % Calcolo delle elevazioni dei satelliti
%     for i = 1:numSat
%         [~, elevation, ~] = lookangles(RecPos_lla_approx, SatPos(i,:),maskAngle);
%         elevation = max(elevation, 0.1); % Elevazione minima per evitare divisioni per zero
%         el_rad = deg2rad(elevation);
% 
%         % Calcolo del ritardo troposferico utilizzando il modello di Saastamoinen
%         Trop_correction(i) = (Kd + Kw) ./ sin(el_rad);
%     end



    % Calculate tropospheric delay using the Saastamoinen model
    % (simplified version, assumes a standard atmosphere)
    lambda = 0.077; % wavelength of GPS signal in meters
    T = 17+273.15; % temperature in Kelvin
    P = 1017; % atmospheric pressure in mbar
    e = 6.11; % water vapor pressure in mbar
    
    h = ecef2lla(RecPos); % receiver altitude in meters
    h = h(3);
    % Calculate tropospheric delay
    Trop_correction = (2.277 * P / T) * (1 + (1255 / T) + (0.05 / h)) * (1 - e / P);
    Trop_correction = Trop_correction * lambda; % convert to meters
    Trop_correction = Trop_correction*ones(numSat,1);

%% APPLY CORRECTIONS

% apply corrections to pseudorange
psr = psr + Clock*dpsr_svclock + Rel*c*Delta_rel_sat + Tropo*Trop_correction;


%% fix RECEIVER POSITION
X = [RecPos 0];

toll = 1e-8;
k_max = 10;
stop_check = 1;
k = 0;

% weigth matrix
WW = diag(ones(1,numSat));

while ((stop_check>toll)&&(k<k_max))

    % update index k
    k = k+1;

    % compute geometric distance between Rec and each Sat
    rho = sqrt(sum((SatPos-RecPos).^2,2));

    % compute psr_c
    psr_c = rho + X(4);

    % compute epsilon
    eps = psr - psr_c;

    % H matrix
    H = NaN(length(numSat),4);
    for j = 1:numSat
        H(j,1) = -(SatPos(j,1)-RecPos(1,1))/rho(j);
        H(j,2) = -(SatPos(j,2)-RecPos(1,2))/rho(j);
        H(j,3) = -(SatPos(j,3)-RecPos(1,3))/rho(j);
        H(:,4) = 1;
    end

    % least square solution
    %Cpp = sqrt(diag(diag(eps*eps')));
    %WW = Cpp\eye(length(Cpp));
    WW = diag(ones(1,numSat));
    AA = H'*WW*H;
    bb = H'*WW*eps;
    DX = (AA\bb)';

    % Update Receiver Position
    X = X + DX;
    RecPos = X(1,1:3);

    stop_check = norm(DX(1:3));

end


%%

RecPos_lla = ecef2lla(RecPos,'WGS84');

rho = sqrt(sum((SatPos-RecPos).^2,2));
psr_c = rho+X(4);
eps = psr-psr_c;
rms = sqrt(eps'*WW*eps/(numSat-numVar));

I = eye(numVar);
COV = (H'*WW*H)\I; % * rms^2;
COV_xyz = COV(1:3,1:3);

% COVARIANCE MATRIX w.r.t. the geodetic coordinates
f = 1/298.257223563;
erad    = 6378137.0 ;   % equatorial radius (meters)
prad    = erad*(1-f);   % polar radius (meters) ; 6356752.3
ecc2 = 1-(prad/erad)^2; % eccentricity al quadrato

% longitude [deg], latitude [deg], altitude [meters]
lat_rad = deg2rad(RecPos_lla(1)); % phi
lon_rad = deg2rad(RecPos_lla(2)); % lambda
alt = RecPos_lla(3);              % h

% Jacobian matrix for the transformation (fi,lam,h) -> (x,y,z)
W_ = sqrt(1-ecc2*(sin(lat_rad))^2);
N = erad/W_;                % NORMAL RADIUS
M = erad*(1-ecc2)/(W_^3);   % MERIDIAN RADIUS

Jac = zeros(3,3);
Jac(1,1) = -(M+alt)*sin(lat_rad)*cos(lon_rad);
Jac(1,2) = -(N+alt)*cos(lat_rad)*sin(lon_rad);
Jac(1,3) = cos(lat_rad)*cos(lon_rad);
Jac(2,1) = -(M+alt)*sin(lat_rad)*sin(lon_rad);
Jac(2,2) = (N+alt)*cos(lat_rad)*cos(lon_rad);
Jac(2,3) = cos(lat_rad)*sin(lon_rad);
Jac(3,1) = -(M+alt)*cos(lat_rad);
Jac(3,2) = 0;
Jac(3,3) = sin(lon_rad);

Jac_1 = Jac\eye(3); % inv(Jac)
COV_geo = Jac_1*COV_xyz*Jac_1'; % Cg nella teoria

% Varianze lla
unc_long = sqrt(COV_geo(1,1));
unc_lat = sqrt(COV_geo(2,2));

% fattore di copertura
lc = 3.7; 
% Varianze ENU
unc_N = lc*M*unc_lat*rms;          %[m]
unc_E = lc*M*unc_long*rms;         %[m]
unc_h = lc*sqrt(COV_geo(3,3))*rms;
unc_tau = lc*sqrt(COV(4,4))*rms;   %[m]

GDOP = sqrt(unc_N^2+unc_E^2+unc_h^2+unc_tau^2)/(lc*rms);
PDOP = sqrt(unc_N^2+unc_E^2+unc_h^2)/(lc*rms);
HDOP = sqrt(unc_N^2+unc_E^2)/(lc*rms);
VDOP = unc_h/(lc*rms);
TDOP = unc_tau/(lc*rms);

% output
unc = [unc_N, unc_E, unc_h];
DOP = [GDOP, PDOP, HDOP, VDOP, TDOP];

end

