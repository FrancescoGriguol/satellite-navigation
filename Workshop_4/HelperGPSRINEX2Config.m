function cfg = HelperGPSRINEX2Config(almfilename,data)
%HelperGPSRINEX2Config Maps the RINEX file elements to the configuration
%object elements.
%
%   Note: This is a helper function and its API and/or functionality may
%   change in subsequent releases.
%
%   CFG = HelperGPSRINEX2Config(ALMFILENAME,DATA) creates configuration
%   object of time HelperGPSNavigationConfig, CFG from the provided DATA.
%   This DATA is output from rinexread method. This function acts to map
%   the elements of such rinexread output to the configuration file
%   appropriately. ALMFILENAME is the almanac file name that is there on
%   the path.

%   Copyright 2022-2023 The MathWorks, Inc.

cfg = HelperGPSNavigationConfig;
for isat = 1:size(data,1)
    % Each of the satellite has got RINEX information that needs to be
    % mapped into the configuration object so as to provide to encoder to
    % generate data bits that must be transmitted.

    cfg(isat).SignalType = "LNAV";
    cfg(isat).AlmanacFileName = almfilename;
    cfg(isat).PRNID = data.SatelliteID(isat);
    cfg(isat).TLMMessage = 0;
    cfg(isat).HOWTOW = round(data.TransmissionTime(isat)/6)+1;
    if data.L2ChannelCodes(isat) == 1
        cfg(isat).CodesOnL2 = "P-code";
    else % if data.L2ChannelCodes == 2
        cfg(isat).CodesOnL2 = "C/A-code";
    end
    cfg(isat).L2PDataFlag = data.L2PDataFlag(isat);

    % CEI config initialization
    cfg(isat).SVHealth = data.SVHealth(isat);
    cfg(isat).IssueOfDataClock = data.IODC(isat);

    % Calculate URAID based on SV accuracy
    uraArr = [0;2.4;3.4;4.85;6.85;9.65;13.65;24;48;96;192;384;768;1536;3072;6144];
    uraid = nnz(uraArr<data.SVAccuracy(isat))-1;
    cfg(isat).URAID = uraid; 
    cfg(isat).WeekNumber = data.GPSWeek(isat);
    cfg(isat).GroupDelayDifferential = data.TGD(isat);
    cfg(isat).SVClockCorrectionCoefficients = [data.SVClockBias(isat);data.SVClockDrift(isat);data.SVClockDriftRate(isat)];

    % Get the time in week from the RINEX information
    d = data.Time(isat);
    d.TimeZone = "UTC";
    [~,timeInWeek] = HelperGPSConvertTime(d);
    cfg(isat).ReferenceTimeOfClock = timeInWeek;
    cfg(isat).SemiMajorAxisLength = data.sqrtA(isat).^2;
    cfg(isat).MeanMotionDifference = data.Delta_n(isat)/pi; % Division by pi to convert radians to semi-circles
    if data.FitInterval(isat) == 4
        cfg(isat).FitIntervalFlag = 0;
    else
        cfg(isat).FitIntervalFlag = 1;
    end

    cfg(isat).Eccentricity = data.Eccentricity(isat);
    cfg(isat).MeanAnomaly = data.M0(isat)/pi; % Division by pi to convert radians to semi-circles
    cfg(isat).ReferenceTimeOfEphemeris = data.Toe(isat);
    cfg(isat).HarmonicCorrectionTerms = [data.Cis(isat); data.Cic(isat); data.Crs(isat); data.Crc(isat); data.Cus(isat); data.Cuc(isat)];
    cfg(isat).IssueOfDataEphemeris = data.IODE(isat);
    cfg(isat).IntegrityStatusFlag = 0;
    cfg(isat).ArgumentOfPerigee = data.omega(isat)/pi; % Division by pi to convert radians to semi-circles
    cfg(isat).RateOfRightAscension = data.OMEGA_DOT(isat)/pi; % Division by pi to convert radians to semi-circles
    cfg(isat).LongitudeOfAscendingNode = data.OMEGA0(isat)/pi; % Division by pi to convert radians to semi-circles
    cfg(isat).Inclination = data.i0(isat)/pi; % Division by pi to convert radians to semi-circles
    cfg(isat).InclinationRate = data.IDOT(isat)/pi; % Division by pi to convert radians to semi-circles
    cfg(isat).AlertFlag = 0;
end
end