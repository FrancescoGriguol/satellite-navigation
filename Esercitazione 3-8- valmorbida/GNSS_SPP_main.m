% Dept. of Industrial Engineering - University of Padova
% MS course: Satellite Navigation - GNSS
% author: A. Valmorbida, Ph.D.
% date: March 27, 2002
% file name: GNSS_SPP_main.m

% clear workspace
clear; close all; clc;

% global variables
global wE c mu g
global GNSS_type
global RecPos_lla_approx RecPos_lla_ref
global Tropo Iono Clock Rel

GNSS_type = "GPS";

% Set 1 to have the correction, 0 to exclude the correction
Tropo = 1;
Iono  = 0; 
Clock = 1;
Rel   = 1;
% Tgd   = 0;

% phisical constants definition
c  = 299792458;         % [m/s]
mu = 3.986005e+14;      % WGS-84 Earth Universal Gravitational parameter (m3/s2)
RE = 6378137;           % mean Earth's radius, [m]
wE = 7.2921151467E-05;  % [rad/s] for GPS and Galileo
LeapSeconds = 18;       % scelto dal 2017
g = 9.80665;            % gravity constant [m/s^2]

%% Load navigation and observation data

addpath("data/");
addpath("Acquisizione_SwiftNav_20240410/");
addpath("brdc_utils\")
addpath("orbit_computation\")
addpath("GNSSdata_29052024_dynamic\")
addpath("GNSSdata_29052024_static\")
obsfilename = "swift-gnss-20240529-132206.sbp.obs";
navfilename = "swift-gnss-20240529-132206.sbp.nav";

% read observation and navigation RINEX files
obsfileinfo = rinexinfo(obsfilename);
obsdata = rinexread(obsfilename);
navdata = rinexread(navfilename);

% compute initial solution in lat-long-alt coordinates
RecPos_lla_approx = ecef2lla(obsfileinfo.ApproxPosition,'WGS84');

% import fiducial solution
load swiftgnss20240529130001_fiducial_solution.mat;
lat_fid = RecPos_lla_ref(:,1);
lon_fid = RecPos_lla_ref(:,2);
alt_fid = RecPos_lla_ref(:,3);
RecPos_ref_totali = lla2ecef(RecPos_lla_ref);
raggio_ref = (RecPos_ref_totali(:,1).^2+RecPos_ref_totali(:,2).^2+RecPos_ref_totali(:,3).^2).^0.5;


% compute reference solution
RecPos_lla_ref = [mean(lat_fid) mean(lon_fid) mean(alt_fid)];
RecPos_ref = lla2ecef(RecPos_lla_ref);
scartolat = lat_fid-RecPos_lla_ref(1);
scartolon = lon_fid-RecPos_lla_ref(2);

%% Plot reference solution

figure()
geoplot(lat_fid,lon_fid,'x','Color',[0.3010 0.7450 0.9330]);
hold on;
geoplot(RecPos_lla_ref(1),RecPos_lla_ref(2),'p','MarkerSize',15,'LineWidth',1,'MarkerFaceColor',[0.4940 0.1840 0.5560],'MarkerEdgeColor','none');
geobasemap satellite % satellite - streets - topographic
legend('Fiducial solutions','Fiducial position','Location','northeast');
title('FIDUCIAL SOLUTION AND RECEIVER POSITION','FontSize',24);
ax = gca;
ax.FontSize = 15;

%% observation data => pseudoranges

% Select the GNSS type data
obsdata = obsdata.(GNSS_type);
navdata = navdata.(GNSS_type);

[~,idx] = unique(obsdata.Time); % della colonna Time elimino tutti i doppioni per capire gli istanti delle misurazioni e non anche le ripetizioni per gli altri satelliti che erano in vista contemporneamente
time_all_obs = obsdata.Time(idx,:); % tiro fuori gli indici dalla colonna senza doppioni
numObs = length(time_all_obs); % numero di istanti temporali in cui sono state prese le misurazioni che è il numero di osservazioni, per ognuno ci sarà il fix di posizione

% memory allocation
RecPos = zeros(numObs,3); % in ecef
RecPos_lla = zeros(numObs,3); % in geodetiche
unc = zeros(numObs,3); % semiassi magg min e theta rotazione ellisse
DOP = zeros(numObs,5); % diluition of precision
ti = (0:numObs-1)'; % tempi
numSat = zeros(numObs,1); % ad ogni istante quanti satelliti sono in vista
SatPos = zeros(numObs,3);
RecPos0 = obsfileinfo.ApproxPosition; % posizione iniziale ecef dal file che ci è stato dato come prima stima
COV_xyz = zeros(3,3,numObs);
for k = 1:numObs

    % get time t_i
    time_i = time_all_obs(k);

    % get obs data at time t_i
    obs_i = obsdata(obsdata.Time == time_i,:);

    % sort data by Satillte ID
    [~,satIdx_obs_i] = unique(obs_i.SatelliteID);
    obs_i = obs_i(satIdx_obs_i,:); % riodino i dati in ordine crescente secondo gli id dei satelliti
    numSat(k) = size(obs_i,1); % quanti satelliti sono in vista in quel momento

    % navigation data => get pos of satellite in view
    [~,satIdx_nav_i] = unique(navdata.SatelliteID);
    navmsg = navdata(satIdx_nav_i,:);

    % fix Pos
    [itsc_satID,iobs,inav] = intersect(obs_i.SatelliteID,navmsg.SatelliteID); % ritorna i dati comuni senza ripetizioni in modo ordinato tra obs e nav
    if length(obs_i.SatelliteID) ~= length(navmsg.SatelliteID)
        obs_i = obs_i(iobs,:); 
        navmsg = navmsg(inav,:);
    end
    
    [RecPos(k,:),RecPos_lla(k,:),unc(k,:),DOP(k,:),COV_xyz] = ...
        GNSS_fix_ti(time_i,obs_i,navmsg,LeapSeconds,RecPos0);

end

rowsWithNaN = any(isnan(RecPos_lla), 2);
RecPos_lla_true = RecPos_lla;
RecPos_lla_true(rowsWithNaN, :) = [];

% Plot degli unc
figure()
subplot(3,1,1); plot(ti,unc(:,2)./3.7,"LineWidth",1,"LineStyle","-","Color",[0.4660 0.6740 0.1880]); xlabel('time [s]'); ylabel('East [m]'); set(gca, 'FontSize', 15);
hold on
plot(ti,mean(unc(:,2)./3.7).*ones(length(ti)),'LineStyle','--','Color','k','LineWidth',1.5)
yticks([0 5 mean(unc(:,2)./3.7) 10])
hold off
subplot(3,1,2); plot(ti,unc(:,1)./3.7,"LineWidth",1,"LineStyle","-","Color",[0 0.4470 0.7410]); xlabel('time [s]'); ylabel('North [m]'); set(gca, 'FontSize', 15);
hold on
plot(ti,mean(unc(:,1)./3.7).*ones(length(ti)),'LineStyle','--','Color','k','LineWidth',1.5)
yticks([0 2 mean(unc(:,1)./3.7) 4])
hold off
subplot(3,1,3); plot(ti,unc(:,3)./3.7,"LineWidth",1,"LineStyle","-","Color",[0.6350 0.0780 0.1840]); xlabel('time [s]'); ylabel('Up [m]'); set(gca, 'FontSize', 15);
hold on
plot(ti,mean(unc(:,3)./3.7).*ones(length(ti)),'LineStyle','--','Color','k','LineWidth',1.5)
yticks([0 5 mean(unc(:,3)./3.7) 10])
hold off
sgtitle('NORTH EAST UP UNCERTAINTY','Fontsize',24);


%% Analisi di Covarianza

% raggio in metri
raggio = (RecPos(:,1).^2+RecPos(:,2).^2+RecPos(:,3).^2).^0.5;

% PER LE STATICHE CONTROLLARE SE LE POSIZIONI VENGONO SPECCHIATE, IN CASO
% CAMBIARE I SEGNI DI ERR ALL'INIZIO. PER L'ACQUISIZIONE FRONTE PIOVEGO
% DECOMMENTARE LA RIGA IN CUI SI TRASPONGONO LE COLONNE DI ERR
err = RecPos_lla_true(:,1:2) - RecPos_lla_ref(1:2);
% err = [err(:,2) err(:,1)]; ROMPERE IN CASO DI NECESSITA'
err = deg2rad(err);
err = err.*raggio;
pos = mean(err);
scartolat = deg2rad(scartolat).*raggio_ref;
scartolon = deg2rad(scartolon).*raggio_ref;

    covariance = cov(err);
    [eigenvector,eigenval] = eig(covariance);

    smallest_eigenval = min(diag(eigenval));
    largest_eigenval = max(diag(eigenval));

    theta1 = atan2(eigenvector(2,2),eigenvector(1,2));
    theta1_deg = rad2deg(theta1);
    R = [cos(theta1) sin(theta1); -sin(theta1) cos(theta1)];

    % uncertainty ellipse
    conf_val = sqrt(-2*log(1-0.997));
    a1 = conf_val*sqrt(largest_eigenval);
    b1 = conf_val*sqrt(smallest_eigenval);

    tt = linspace(0,2*pi,100);

    uu1 = a1*cos(tt);
    vv1 = b1*sin(tt);

    uu_th1 =  uu1*cos(theta1) + vv1*sin(theta1);
    vv_th1 = -uu1*sin(theta1) + vv1*cos(theta1);

    xx1 = uu_th1 + pos(1);
    yy1 = vv_th1 + pos(2);

    figure()
    plot(err(:,1),err(:,2),'x','Color',[0.9290 0.6940 0.1250])
    hold on
    plot(scartolat,scartolon,'x','Color',[0.3010 0.7450 0.9330])
    plot(0,0,'p','MarkerSize',14,'LineWidth',1,'MarkerEdgeColor','none','MarkerFaceColor',[0.4940 0.1840 0.5560])
    plot(pos(1),pos(2),'p','MarkerSize',14,'LineWidth',1,'MarkerEdgeColor','none','MarkerFaceColor','g')
    plot(xx1,yy1,'-','LineWidth',1,'Color',[0.6350 0.0780 0.1840])
    title('POSITION SOLUTION AND UNCERTAINTY ELLIPSE','FontSize',24)
    subtitle('with respect from fiducial position','FontSize',18)
    legend('estimated solution','fiducial solution','fiducial position','estimated position','uncertainty ellipse','Location','northwest')
    grid on
    ax = gca;
    ax.FontSize = 15;
    xlabel('East [m]')
    ylabel('North [m]')
    axis equal

%% Post-processing: geoplot

figure()
geoplot(lat_fid,lon_fid,'x','Color',[0.3010 0.7450 0.9330]);
hold on;
geoplot(RecPos_lla_true(:,1),RecPos_lla_true(:,2),'x','Color',[0.9290 0.6940 0.1250]);
hold on
geoplot(RecPos_lla_ref(1),RecPos_lla_ref(2),'p','MarkerSize',15,'LineWidth',1,'MarkerFaceColor',[0.4940 0.1840 0.5560],'MarkerEdgeColor','none');
hold on
geoplot(mean(RecPos_lla_true(:,1)),mean(RecPos_lla_true(:,2)),'p','MarkerSize',15,'LineWidth',1,'MarkerFaceColor','g','MarkerEdgeColor','none');
geobasemap satellite % satellite - streets - topography
legend('Fiducial solutions','Estimated solutions','Fiducial position','Estimated position','Location','northeast');
title('ESTIMATED SOLUTIONS AND  RECEIVER POSITION OUTCOME','FontSize',24)
ax = gca;
ax.FontSize = 15;

%% Plot DOP
figure()
subplot(2,3,1); plot(ti,DOP(:,1)); xlabel('time [s]'); ylabel('GDOP'); set(gca, 'FontSize', 15);
subplot(2,3,2); plot(ti,DOP(:,2)); xlabel('time [s]'); ylabel('PDOP'); set(gca, 'FontSize', 15);
subplot(2,3,3); plot(ti,DOP(:,3)); xlabel('time [s]'); ylabel('HDOP'); set(gca, 'FontSize', 15);
subplot(2,3,4); plot(ti,DOP(:,4)); xlabel('time [s]'); ylabel('VDOP'); set(gca, 'FontSize', 15);
subplot(2,3,5); plot(ti,DOP(:,5)); xlabel('time [s]'); ylabel('TDOP'); set(gca, 'FontSize', 15);
subplot(2,3,6); plot(ti,numSat); xlabel('time [s]'); ylabel('Num Sat'); set(gca, 'FontSize', 15);
hold on 
plot(ti,4.*ones(length(ti)),'r--')
sgtitle('DOP PARAMETERS AND NUMBER OF SATELLITES','Fontsize',24);
yticks(0:1:12)
ylim([1 12])

%% Skyplot

% Satellite ID (choose satellite ID to study)
SatID = 1:32;

% Time to interpolate
toe_start = obsfileinfo.FirstObsTime;
toe_start = toe_start.Hour*3600+toe_start.Minute*60+toe_start.Second;
toe_end = obsfileinfo.LastObsTime;
toe_end = toe_end.Hour*3600+toe_end.Minute*60+toe_end.Second; 
vector_time = linspace(0,(toe_end-toe_start)/3600,numObs);  % [hours]
interp_time = vector_time;

% Read broadcast file
brdc = rinexread(navfilename); 
info = rinexinfo(navfilename);
% Ephemeris of the selected constellation
brdcEphem = brdc.(GNSS_type);

[SAT,SatID] = brdc_get_sc_data(brdc,SatID,brdcEphem,GNSS_type);

% define mask angle
maskAngle = 15; % [deg]

orbit = zeros(length(vector_time),3);
fname = fieldnames(SAT);

if GNSS_type == "GPS"
    for i = 1 : length(SatID)
        t = toe_start + vector_time*3600;
        for j = 1 : length(t)
            ephem = SAT.(['SatID_', num2str(SatID(i),'%i')]);
            % get coordinates (1-vuol dire prima riga, si estrae già la colonna giusta)
            [Xk,Yk,Zk] = GPS_coordinates(t(j), ephem.Crs(1), ephem.Delta_n(1), ...
                ephem.M0(1), ephem.Cuc(1), ephem.Eccentricity(1), ephem.Cus(1), ...
                ephem.sqrtA(1), ephem.Toe(1), ephem.Cic(1), ephem.OMEGA0(1), ...
                ephem.Cis(1), ephem.i0(1), ephem.Crc(1), ...
                ephem.omega(1), ephem.OMEGA_DOT(1), ephem.IDOT(1));
            orbit(j,:) = [Xk,Yk,Zk]*1e-3; % per salvare i dati in km
        end
        SAT.(string([char(fname(i)),'_brdc'])) = orbit;
    end

elseif GNSS_type == "GLONASS"
    Leap_Seconds = info.LeapSecondParameters.LeapSeconds; % [s]
    for i = 1 : length(SatID)
        t = interp_time*3600-Leap_Seconds;
        ephem = SAT.(['SatID_', num2str(SatID(i),'%i')]);
        %  get GLONASS coordinates
        [x,y,z] = GLONASS_coordinates(ephem(1,:),t);
        orbit = [x,y,z]; % [km]
        SAT.(string([char(fname(i)),'_brdc'])) = orbit;

    end
end

% mem alloc
Az_matrix = NaN(numel(t),numel(SatID));
El_matrix = NaN(numel(t),numel(SatID));
Vis_matrix = NaN(numel(t),numel(SatID));

for i = 1:length(SatID)
    [az,el,vis] = lookangles(RecPos_lla_ref,SAT.(string([char(fname(i)),'_brdc']))*1e3,maskAngle);
    Az_matrix(:,i) = az;
    El_matrix(:,i) = el;
    Vis_matrix(:,i) = vis;
end

% exclude negative elevation
El_matrix(El_matrix<0) = missing;

% Plot
ns = [];
for i = 1:length(t)
    n = find(Vis_matrix(i,:));
    ns = [ns length(n)];
end
figure
plot(vector_time*60*60,ns,'LineWidth',1.5,'LineStyle','-')
title('SATELLITES IN SIGHT',GNSS_type,'FontSize',24)
xlabel('Time of observation [s]')
ylabel('Number of satellites')
yticks(0:1:11)
grid on
ax = gca;
ax.FontSize = 15;

% Skyplot
figure
sp = skyplot(Az_matrix(1,:),El_matrix(1,:),SatID,MaskElevation=maskAngle);
set(sp,'labelfontsize',20)
for idx = 1:size(Az_matrix,1)
set(sp,AzimuthData=Az_matrix(1:idx,:),ElevationData=El_matrix(1:idx,:));
drawnow limitrate
%pause(0.0002)
end
title('SKYPLOT',GNSS_type,'FontSize',24)