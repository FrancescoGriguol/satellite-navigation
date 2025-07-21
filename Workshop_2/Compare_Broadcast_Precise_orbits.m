% Workspace setting

clear; close all; clc;

addpath(genpath('.\'));

%% Initial settings 

% Time to interpolate
interp_time = linspace(1/6,1.2,1000);  % [hours]

% Precise - time 
[gpsweek, tow, doy, dow] = greg2gps([2024,1,31,0,0,0]);
% Broadcast - time
toe_start = 259200; 
vector_time = interp_time;

% SP3 files available at http://navigation-office.esa.int/products/gnss-products/
filename_precise = 'ESA0MGNFIN_20240310000_01D_05M_ORB.SP3';
% RINEX file available at https://cddis.nasa.gov/archive/gnss/data/daily/
filename_broadcast = 'BRDC00IGS_R_20240310000_01D_MN.rnx';

% Selected constellation and satellite
constID = "GPS";
SatID = 4;

%% Read data from SP3 file and get coordinates from sp3 data format

% Read SP3 file
[sp3, Greg_time] = read_sp3_multiconstellation(filename_precise,constID);
[SAT_pr,SatID] = sp3_get_sc_pos(sp3,SatID);

%% Interpolation of SP3

fname = fieldnames(SAT_pr);
t = (juliandate(Greg_time)-juliandate(Greg_time(1,:)))*86400; % time [sec]
% create function that interpolates sat pos based on time to interp, GS_ID
for i = 1:length(SatID)
    [orbit_pr] = interp_precise_orbits(SAT_pr.(string(fname(i))),t(1:11),interp_time);
    SAT_pr.(string([char(fname(i)),'_interp'])) = orbit_pr;
end

%%  Read data from RINEX OBSERVATION FILE

% Read broadcast file
brdc = rinexread(filename_broadcast);
info = rinexinfo(filename_broadcast);

% Ephemeris of the selected constellation
brdcEphem = brdc.(constID);

%% Study selected satellites

[SAT_br,SatID] = brdc_get_sc_data(brdc,SatID,brdcEphem,constID);

%% Orbit computation 

orbit_br = zeros(length(vector_time),3);
fname = fieldnames(SAT_br);

if constID == "GPS"
    for i = 1 : length(SatID)
        t = toe_start + vector_time*3600;
        for j = 1 : length(t)
            ephem = SAT_br.(['SatID_', num2str(SatID(i),'%i')]);
            % get coordinates (1-vuol dire prima riga, si estrae già la colonna giusta)
            [Xk,Yk,Zk] = GPS_coordinates(t(j), ephem.Crs(1), ephem.Delta_n(1), ...
                ephem.M0(1), ephem.Cuc(1), ephem.Eccentricity(1), ephem.Cus(1), ...
                ephem.sqrtA(1), ephem.Toe(1), ephem.Cic(1), ephem.OMEGA0(1), ...
                ephem.Cis(1), ephem.i0(1), ephem.Crc(1), ...
                ephem.omega(1), ephem.OMEGA_DOT(1), ephem.IDOT(1));
            orbit_br(j,:) = [Xk,Yk,Zk]*1e-3; % per salvare i dati in km
        end
        SAT_br.(string([char(fname(i)),'_brdc'])) = orbit_br;
    end

elseif constID == "GLONASS"
    Leap_Seconds = info.LeapSecondParameters.LeapSeconds; % [s]
    for i = 1 : length(SatID)
        t = interp_time*3600-Leap_Seconds;
        ephem = SAT_br.(['SatID_', num2str(SatID(i),'%i')]);
        %  get GLONASS coordinates
        [x,y,z] = GLONASS_coordinates(ephem(1,:),t);
        orbit_br = [x,y,z]; % [km]
        SAT_br.(string([char(fname(i)),'_brdc'])) = orbit_br;

    end
end

%% Plot

eps = zeros(length(vector_time),3);
for i = 1:length(SatID)
    eps = SAT_pr.(string([char(fname(i)),'_interp']))-SAT_br.(string([char(fname(i)),'_brdc']));
end
figure
if interp_time(end)<1.5
    plot(interp_time*60,eps*1000,'LineWidth',1)
    title('Confronto','FontSize',20)
    subtitle('Broadcast Navigation Message - Precise Products','FontSize',18)
    legend('Direzione X','Direzione Y','Direzione Z','FontSize',13)
    xlabel('Tempo di propagazione [min]','FontSize',16)
    ylabel('Differenza di posizione [m]','FontSize',16)
else
    plot(interp_time,eps,'LineWidth',1)
    title('Confronto','FontSize',20)
    subtitle('Broadcast Navigation Message - Precise Products','FontSize',18)
    legend('Direzione X','Direzione Y','Direzione Z','FontSize',13)
    xlabel('Tempo di propagazione [h]','FontSize',16)
    ylabel('Differenza di posizione [km]','FontSize',16)
end





