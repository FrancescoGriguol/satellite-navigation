function almanac = HelperGPSAlmanac2Struct(filename)
%Convert the almanac in SEM format to structure
%
%   Note: This is a helper and its API and/or functionality may change
%   in subsequent releases.
%
%   ALMANAC = HelperGPSAlmanac2Struct(FILENAME) converts the almanac
%   information that is available in FILENAME file to structure ALMANAC.
%   FILENAME should have data in System Effective Model (SEM) format.
%
% Parameters in ALMANAC
%   SatellitesPRNIndices
%   WeekNumModulo1024
%   ReferenceTimeOfAlmanac
%   Data.PRNID
%   Data.SVID
%   Data.URAID
%   Data.Eccentricity
%   Data.InclinationDifference
%   Data.RateOfRightAscension
%   Data.RootOfSemiMajorAxis
%   Data.LongitudeOfAscendingNode
%   Data.ArgumentOfPerigee
%   Data.MeanAnomaly
%   Data.SVClockCorrectionCoefficients % [af0, af1] in that order
%   Data.SVHealth % 6 bits with MSB value to be health status flag
%   Data.NAVDataHealth % 3 bit NAV data health indications. If not specified,
%                      % default value of 3 three zero bits is taken
%   Data.SVConfig

%   Copyright 2020-2021 The MathWorks, Inc.

fgps = fopen(filename, 'r');
formatSpec = '%lf';
numSatelliteRecords = fscanf(fgps,formatSpec, [1 1]); % First value is number of satellite records
notUsed = fscanf(fgps,'%s', [1 1]); %#ok<NASGU>  % Second value is name of 
                                                 % almanac which is not
                                                 % used here. Initializing
                                                 % this parameter for
                                                 % carrying the file
                                                 % pointer forward

almanac.SatellitesPRNIndices = zeros(numSatelliteRecords, 1);
almanac.WeekNumModulo1024 = fscanf(fgps, formatSpec, [1 1]);
almanac.ReferenceTimeOfAlmanac = fscanf(fgps, formatSpec, [1 1]);
for iAlmanac = 1:numSatelliteRecords
    % Read Satellite PRN Index (PRNIdx)
    prnidx = fscanf(fgps,formatSpec, [1 1]);
    almanac.SatellitesPRNIndices(iAlmanac) = prnidx;
    almanac.Data(iAlmanac).PRNID = prnidx;
    
    % Read SV ID (SVID)
    svid = fscanf(fgps,formatSpec, [1 1]);
    almanac.Data(iAlmanac).SVID = svid;
    
    % Read Average URA Number (URAIdx)
    uraidx = fscanf(fgps,formatSpec, [1 1]);
    almanac.Data(iAlmanac).URAID = uraidx;
    
    % Read eccentricity (e)
    e = fscanf(fgps,formatSpec, [1 1]);
    almanac.Data(iAlmanac).Eccentricity = e;
    
    % Read Inclination offset (delta_i)
    delta_i = fscanf(fgps,formatSpec, [1 1]);
    almanac.Data(iAlmanac).InclinationDifference = delta_i;
    
    % Read Rate of Right ascension (OmegaDot)
    OmegaDot = fscanf(fgps,formatSpec, [1 1]);
    almanac.Data(iAlmanac).RateOfRightAscension = OmegaDot;
    
    % Read Square root of Semi-major axis (SqrtA)
    SqrtA = fscanf(fgps,formatSpec, [1 1]);
    almanac.Data(iAlmanac).RootOfSemiMajorAxis = SqrtA;
    
    % Read Geographic longitude of Orbital Plane (Omega_0)
    Omega_0 = fscanf(fgps,formatSpec, [1 1]);
    almanac.Data(iAlmanac).LongitudeOfAscendingNode = Omega_0;
    
    % Read Argument of Perigee (omega)
    omega = fscanf(fgps,formatSpec, [1 1]);
    almanac.Data(iAlmanac).ArgumentOfPerigee = omega;
    
    % Read mean anomaly (M_0)
    M_0 = fscanf(fgps,formatSpec, [1 1]);
    almanac.Data(iAlmanac).MeanAnomaly = M_0;
    
    % Read zeroth order clock correction (a_f0)
    a_f0 = fscanf(fgps,formatSpec, [1 1]);
    almanac.Data(iAlmanac).SVClockCorrectionCoefficients = a_f0;
    
    % Read first order clock correction (a_f1)
    a_f1 = fscanf(fgps,formatSpec, [1 1]);
    almanac.Data(iAlmanac).SVClockCorrectionCoefficients(2) = a_f1;
    
    % Read satellite health (SVHealth)
    SVHealth = fscanf(fgps,formatSpec, [1 1]);
    almanac.Data(iAlmanac).SVHealth = SVHealth;
    
    % Read Satellite configuration (SVConfig)
    SVConfig = fscanf(fgps,formatSpec, [1 1]);
    almanac.Data(iAlmanac).SVConfig = SVConfig;
end
fclose(fgps);
end
