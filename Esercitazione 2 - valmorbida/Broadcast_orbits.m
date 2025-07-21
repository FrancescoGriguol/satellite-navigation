% workspace setting

close all; clear; clc;

addpath(genpath('.\'));


%% Initial settings

% file available at https://cddis.nasa.gov/archive/gnss/data/daily/
filename = 'BRDC00IGS_R_20240310000_01D_MN.rnx';

% Selected constellation
constID = "GLONASS";

% Satellite ID (choose satellite ID to study)
SatID = [1,15,30];

% Time to interpolate
vector_time = linspace(0,3,1000);  % [hours]
toe_start = 259200;                % tempo iniziale delle effemeridi
interp_time = vector_time;

%%  Read data from RINEX OBSERVATION FILE

% Read broadcast file
brdc = rinexread(filename); % function solo per rinex da versione 3 in poi
info = rinexinfo(filename);

% Ephemeris of the selected constellation
brdcEphem = brdc.(constID);


%% Study selected satellites

[SAT,SatID] = brdc_get_sc_data(brdc,SatID,brdcEphem,constID);


%% Orbit computation 

orbit = zeros(length(vector_time),3);
fname = fieldnames(SAT);

if constID == "GPS"
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

elseif constID == "GLONASS"
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


%% Visualize orbits

sc = satelliteScenario;
sc.StartTime = SAT.SatID_1.Time(1);
sc.StopTime = SAT.SatID_1.Time(1)+hours(2);
% sat
fname = fieldnames(SAT);
t = interp_time*3600;
for i = 1: length(SatID)
    postimeseries = timeseries(SAT.(string([char(fname(i)),'_brdc']))*1e3,t);
    satellite(sc,postimeseries,"CoordinateFrame","ecef","Name",fname(i));
end
satelliteScenarioViewer(sc);



