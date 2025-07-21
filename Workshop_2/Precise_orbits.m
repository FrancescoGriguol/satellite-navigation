% workspace setting

close all; clear; clc;

addpath(genpath('.\')); 

%% Time

[gpsweek, tow, doy, dow] = greg2gps([2024,1,31,0,0,0]);

%% Initial settings

% SP3 files available at http://navigation-office.esa.int/products/gnss-products/
filename = 'ESA0MGNFIN_20240310000_01D_05M_ORB.SP3';

% Selected constellation (ID della costellazione che vogliamo analizzare)
constID = "GPS";

% Satellite ID (choose satellite ID to study)
% si mettono in un vettore i satelliti che si vogliono studiare
SatID = [1,15,30];

% Time to interpolate
% si setta un tempo per l'interpolazione dei dati, e per la stabilità
% dell'interpolazione si esprime il tempo in ore
interp_time = linspace(0,3,1000);  % [hours]

%% Read data from SP3 file

% Read SP3 file
[sp3, Greg_time] = read_sp3_multiconstellation(filename,constID);

%% Get satellite coordinates from sp3 data format
% estrarre i dati per l'interpolazione
[SAT,SatID] = sp3_get_sc_pos(sp3,SatID);

%% Visualize orbits

% set Satellite Scenario Viewer
sc = satelliteScenario;
sc.StartTime = datetime(Greg_time(1,:));
sc.StopTime = datetime(Greg_time(1,:)) + days(1);

% add satellites to satelliteScenario
fname = fieldnames(SAT);
t = (juliandate(Greg_time)-juliandate(Greg_time(1,:)))*86400; % time [sec]
for i = 1 : length(SatID)
    postimeseries = timeseries(SAT.(string(fname(i)))*1e3,t);
    satellite(sc,postimeseries,"CoordinateFrame","ecef","Name",fname(i));
end

% display Satellite Scenario Viewer
satelliteScenarioViewer(sc);


%% Interpolation of SP3

% create function that interpolates sat pos based on time to interp, GS_ID
for i = 1:length(SatID)
    [orbit] = interp_precise_orbits(SAT.(string(fname(i))),t(1:11),interp_time);
    SAT.(string([char(fname(i)),'_interp'])) = orbit;
end


%% Compare SP3 data and interpolated data 

figure();

plot(postimeseries.Time,SAT.SatID_1(:,1),'o');
hold on
plot(interp_time*3600,SAT.SatID_1_interp(:,1));
legend('SP3 points', 'interp');
xlabel('time [s]');
ylabel('X [km]');

