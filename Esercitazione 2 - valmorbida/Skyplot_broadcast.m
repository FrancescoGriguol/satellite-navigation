%% skyplot broadcast

% workspace setting

close all; clear; clc;

addpath(genpath('.\')); 

%% Initial settings

% file available at https://cddis.nasa.gov/archive/gnss/data/daily/
filename = 'BRDC00IGS_R_20240310000_01D_MN.rnx';

% Selected constellation
constID = "GPS";

% Satellite ID (choose satellite ID to study)
SatID = 1:32;

% Time to interpolate
vector_time = linspace(1/6,1.2,1000);  % [hours]
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

%% Initialize skyplot

% define mask angle
maskAngle = 15; % [deg]

% compute reference solution
load('swift-gnss-20221206-124503-fiducial_solution.mat')
RecPos_lla_ref = [mean(lat_fid) mean(lon_fid) mean(hei_fid)];

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

%% Plot
ns = [];
for i = 1:length(t)
    n = find(Vis_matrix(i,:));
    ns = [ns length(n)];
end
figure
plot(vector_time*60,ns,'LineWidth',1)
title('SATELLITI IN VISTA',constID,'FontSize',17)
xlabel('Tempo di osservazione [min]','FontSize',13)
ylabel('Numero satelliti','FontSize',13)
grid on
xticks(0:15:180)
yticks(5:1:12)

%% Skyplot
figure
sp = skyplot(Az_matrix(1,:),El_matrix(1,:),SatID,MaskElevation=maskAngle);
set(sp,'labelfontsize',20)
for idx = 1:size(Az_matrix,1)
set(sp,AzimuthData=Az_matrix(1:idx,:),ElevationData=El_matrix(1:idx,:));
drawnow limitrate
pause(0.05)
end
title('SKYPLOT',constID,'FontSize',17)