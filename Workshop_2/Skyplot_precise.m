%% skyplot precise

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
SatID = 1:1:32;

%% Read data from SP3 file and get satellite coordinates from sp3 data format
[sp3, Greg_time] = read_sp3_multiconstellation(filename,constID);
[SAT,SatID] = sp3_get_sc_pos(sp3,SatID);

%% Initialize skyplot

% define mask angle
maskAngle = 15; % [deg]

% compute reference solution
load('swift-gnss-20221206-124503-fiducial_solution.mat')
RecPos_lla_ref = [mean(lat_fid) mean(lon_fid) mean(hei_fid)];

t = (juliandate(Greg_time)-juliandate(Greg_time(1,:)))*86400; % time [sec]

% mem alloc
Az_matrix = NaN(numel(t),numel(SatID));
El_matrix = NaN(numel(t),numel(SatID));
Vis_matrix = NaN(numel(t),numel(SatID));

fname = fieldnames(SAT);

for i = 1:length(SatID)
    [az,el,vis] = lookangles(RecPos_lla_ref,SAT.(string(fname(i)))*1e3,maskAngle);
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
plot(t/3600,ns,'LineWidth',1)
title('SATELLITI IN VISTA',constID,'FontSize',17)
xlabel('Tempo di osservazione [h]','FontSize',13)
ylabel('Numero satelliti','FontSize',13)
xlim([0 24])
grid on

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