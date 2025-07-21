% Dept. of Industrial Engineering - University of Padova
% MS course: Satellite Navigation - GNSS
% author: S. Chiodini, Ph,D; A. Valmorbida, Ph.D
% date: May 7, 2024
% file name: GNSS_LowDyn_KF.m
% -------------------------------------------------------------------------
%% Settings 
clear 
close all
clc

addpath("Acquisizione_SwiftNav_20240410\")
obsdata = rinexread('swift-gnss-20240410-130031.sbp.obs');
obsdata = obsdata.GPS;
obsfileinfo = rinexinfo('swift-gnss-20240410-130031.sbp.obs');

navdata = rinexread('swift-gnss-20240410-130031.sbp.nav');
navdata = navdata.GPS;

% phisical constants definition
c  = 299792458;         % [m/s]
mu = 3.986005e+14;      % WGS-84 Earth Universal Gravitational parameter (m3/s2)
RE = 6378137;           % mean Earth's radius, [m]
wE = 7.2921151467E-05;  % [rad/s] for GPS and Galileo
LeapSeconds = 18;

[~,idx] = unique(obsdata.Time);
% Tempi corrispondendti ad ogni osservazione, avendo considerato un unico
% valore di tempo quando c'erano più satelliti in vista contemporaneamente
time_all_obs = obsdata.Time(idx,:);
% numObs corrisponde al numero di tempi in cui sono  state fatte le
% osservazioni 
numObs = length(time_all_obs);

%% Initialization

% memory allocation
RecPos = zeros(numObs,3);
RecPos_lla = zeros(numObs,3);
unc = zeros(numObs,3);
DOP = zeros(numObs,5);
ti = (0:numObs-1)';
numSat = zeros(numObs,1);

% Posizione approssimativa del ricevitore presa dai dati di osservazione
RecPos0 = obsfileinfo.ApproxPosition;
% lla = ecef2lla(RecPos0);
% 
% figure()
% geoplot(lla(1),lla(2),'.')
% geobasemap topographic

% Discrete time state transition equation that is used for a low dynamic
% receiver case is:
% x(k+1) = F(k)x(k)+Gamma(k)v(k)

% Stato iniziale
% x = [x dx y dy z dz r_tr dr_tr]
% dove x,y,z sono le posizioni del low dynamic receiver e le dx,dy,dz le
% sue velocità. Mentre r e dr sono c*tr
X0 = [RecPos0(1) 0 RecPos0(2) 0 RecPos0(3) 0 0 0]';

% Matrice di stato inizializzata col vettore di stato iniziale
X(:,1) = X0;
X_(:,1) = [0 0 0 0 0 0 0 0]';
% integration time
T = 1;

% state transition matrix
F = [1 T 0 0 0 0 0 0;...
     0 1 0 0 0 0 0 0;...
     0 0 1 T 0 0 0 0;...
     0 0 0 1 0 0 0 0;...
     0 0 0 0 1 T 0 0;...
     0 0 0 0 0 1 0 0;...
     0 0 0 0 0 0 1 T;...
     0 0 0 0 0 0 0 1];

sigmaPseudo = 5; % [m] Pseudorange variance 
sigmaV = 0.1; % [m/s] Velocity noise
sigmaA = 0.2; % [m/s^2]  Accel
Qd22 = 0.0039;

% Process noise covariance matrix
Qv = [sigmaA^1        0        0      0      0;...
             0 sigmaA^2        0      0      0;...
             0        0 sigmaA^2      0      0;...
             0        0        0 0.0114 0.0019;...
             0        0        0 0.0019 0.0039];

% Noise gain at time k
Gamma = [0.5*T^2       0       0 0 0;...
               T       0       0 0 0;...
               0 0.5*T^2       0 0 0;...
               0       T       0 0 0;...
               0       0 0.5*T^2 0 0;...
               0       0       T 0 0;...
               0       0       0 1 0;...
               0       0       0 0 1];

% resulting process noise covariance matrix
QGamma = Gamma*Qv*Gamma';

%% Ciclo filtro di Kalman
for k = 1:numObs

    % get time t_i (è il timestep corrente)
    time_i = time_all_obs(k);

    % get obs data at time t_i
    obs_i = obsdata(obsdata.Time == time_i,:);

    % sort data by Satillte ID
    [~,satIdx_obs_i] = unique(obs_i.SatelliteID);
    obs_i = obs_i(satIdx_obs_i,:);
    numSat(k) = size(obs_i,1);

    % navigation data => get pos of satellite in view
    [~,satIdx_nav_i] = unique(navdata.SatelliteID);
    navmsg = navdata(satIdx_nav_i,:);

    % fix Pos
    [itsc_satID,iobs,inav] = intersect(obs_i.SatelliteID,navmsg.SatelliteID);
    if length(obs_i.SatelliteID) ~= length(navmsg.SatelliteID)
        obs_i = obs_i(iobs,:); 
        navmsg = navmsg(inav,:);
    end
    
    % Get satellite position and velocity from navigation message in ECEF
    % reference frame
    [satPos,satVel,satID] = gnssconstellation(time_i,navmsg);

    %% Get pseudoranges and apply corrections
    
    psr1 = obs_i.C1C; % pseudorange relativo alla portante L1
    psr2 = obs_i.C2S; % pseudorange relativo alla portante L2
    f1 = 10.23*10^6*154; % [Hz]
    f2 = 10.23*10^6*120; % [Hz]
    
    % Ionospheric Correction
    a2 = -(f2^2)/(f1^2-f2^2); 
    a1 = 1 - a2;
    psr = a1*psr1 + a2*psr2;
    
    % Satellite Clock correction
    DeltaTime = seconds(obs_i.Time - navmsg.Time);
    dpsr_svclock = c.*(navmsg.SVClockBias + navmsg.SVClockDrift.*DeltaTime + ...
                   navmsg.SVClockDriftRate.*DeltaTime.^2 );
    
    % Relativistic Clock Correction
    Delta_rel_sat = zeros(numSat(k),1);
    for i = 1:numSat(k)
        deltaT = seconds(time_i-navmsg.Time(i))+3600;
        orbitalPeriod = 2*pi*sqrt(navmsg.sqrtA(i)^6/mu);
        meanAnomaly = navmsg.M0(i)+2*pi*orbitalPeriod/deltaT;
        E = kepler_E(navmsg.Eccentricity(i),meanAnomaly);
        Delta_rel_sat(i) = -2*mu^0.5*navmsg.sqrtA(i)*navmsg.Eccentricity(i)*sin(E)/c^2;
    end
    
    % Compute satellite positions
    time_flight = psr/c+dpsr_svclock/c+Delta_rel_sat;

    for i = 1:numSat(k)
        theta = 1*wE*time_flight(i);
        satPos(i,1) = satPos(i,1)*cos(theta)+satPos(i,2)*sin(theta);
        satPos(i,2) = -satPos(i,1)*sin(theta)+satPos(i,2)*cos(theta);
        satPos(i,3) = satPos(i,3);
    end

    idx_nan = find(isnan(psr));
    psr(idx_nan) = psr1(idx_nan);
    
    % Tropospheric correction

    % Measured pseudoranges (with applied corrections)
    psr = psr+dpsr_svclock+c*Delta_rel_sat;
    
    % Noise Covariance matrix assuming each satellite with the same measurement
    % noise variance 
    R = diag(sigmaPseudo*ones(1,length(psr)));
    
    % Jacobian matrix for pseudoranges
    J = zeros(numSat(k),4);
    for i = 1:numSat(k)
        r = sqrt((satPos(i,1)-X(1,k))^2+(satPos(i,2)-X(3,k))^2+(satPos(i,3)-X(5,k))^2);
        J(i,1) = -(satPos(i,1)-X(1,k))/r;
        J(i,2) = -(satPos(i,2)-X(3,k))/r;
        J(i,3) = -(satPos(i,3)-X(5,k))/r;
        J(i,4) = 1;
    end
   
    % Covariance matrix P
    PA = (J'*R*J)';

    P = [PA(1,1)     0      PA(1,2)    0       PA(1,3)    0        PA(1,4)      0;
           0     sigmaV^2     0        0         0        0          0          0;
         PA(2,1)     0      PA(2,2)    0       PA(2,3)    0        PA(2,4)      0;
           0         0         0     sigmaV^2    0        0          0          0;
         PA(3,1)     0      PA(3,2)    0       PA(3,3)    0        PA(3,4)      0;
           0         0         0       0         0      sigmaV^2     0          0;
         PA(4,1)     0      PA(4,2)    0       PA(4,3)    0        PA(4,4)      0;
         0           0         0       0         0          0          0        2*Qd22];
    
    %% State prediction of the receiver
    
    % Input u(k)==0
    % Predict state estimate
    X_(:,k+1) = F*X(:,k);
    % Predict covariance estimate
    % P_  = PA + Qv; 

    

       


end





















%% Plot finale
lla = ecef2lla([X(1,1:end)',X(3,1:end)',X(5,1:end)']);

figure()
geoplot(lla(:,1),lla(:,2),'.')
geobasemap topographic

% %% Post-processing: geoplot
% 
% figure()
% geoplot(RecPos_lla(:,1),RecPos_lla(:,2),'x');
% hold on;
% geoplot(mean(RecPos_lla(:,1)),mean(RecPos_lla(:,2)),'pm','MarkerSize',10);
% geoplot(lat_fid,lon_fid,'xg');
% geoplot(RecPos_lla_ref(1),RecPos_lla_ref(2),'pr','MarkerSize',10);
% geobasemap satellite % satellite - streets - topography
% legend('estimated solutions','estimated position','fiducial solutions','fiducial position');