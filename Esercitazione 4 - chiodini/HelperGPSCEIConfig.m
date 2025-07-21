classdef HelperGPSCEIConfig < comm.internal.ConfigBase
    %HelperGPSCEIConfig GPS clock, ephemeris and integrity data
    %configuration object
    %
    %   Note: This is a helper and its API and/or functionality may change
    %   in subsequent releases.
    %
    %   CFG = HelperGPSCEIConfig creates a Global Positioning System (GPS)
    %   clock, ephemeris, integrity (CEI) configuration object. This helper
    %   file contains configuration properties that are given in table
    %   6-I-1 in IS-GPS-200L [1].
    %
    %   CFG = HelperGPSCEIConfig(Name,Value) creates a GPS CEI data object,
    %   CFG, with the specified property Name set to the specified Value.
    %   You can specify additional name-value pair arguments in any order
    %   as (Name1,Value1,...,NameN,ValueN).
    %
    %   HelperGPSCEIConfig properties:
    %
    %   SVHealth                          - Satellite vehicle health
    %   SignalHealth                      - L1/L2/L5 signal health
    %                                       indicators
    %   IssueOfDataClock                  - Issue of data clock
    %   URAID                             - User range accuracy index
    %   WeekNumber                        - GPS week number
    %   GroupDelayDifferential            - Group delay differential
    %   SVClockCorrectionCoefficients     - Clock correction coefficients,
    %                                       af0, af1, af2
    %   ReferenceTimeOfClock              - Reference time of clock
    %   SemiMajorAxisLength               - Semi-major axis length of the
    %                                       satellite orbit
    %   ChangeRateInSemiMajorAxis         - Rate of change of semi-major
    %                                       axis length
    %   MeanMotionDifference              - Mean motion difference
    %   RateOfMeanMotionDifference        - Rate of change of mean motion
    %                                       difference
    %   FitIntervalFlag                   - Fit interval flag
    %   Eccentricity                      - Eccentricity of the satellite
    %                                       orbit
    %   MeanAnomaly                       - Mean anomaly at reference time
    %   ReferenceTimeOfEphemeris          - Reference time of ephemeris
    %   HarmonicCorrectionTerms           - Six harmonic correction terms
    %   IssueOfDataEphemeris              - Issue of data ephemeris
    %   IntegrityStatusFlag               - Integrity status flag
    %   ArgumentOfPerigee                 - Argument of perigee at
    %                                       reference time
    %   RateOfRightAscension              - Rate of right ascension
    %   LongitudeOfAscendingNode          - Longitude of ascending node
    %   Inclination                       - Inclination angle of the
    %                                       satellite orbit with respect to
    %                                       equator of Earth
    %   InclinationRate                   - Rate of change of inclination
    %                                       angle
    %   URAEDID                           - Elevation dependent user range
    %                                       accuracy index
    %   InterSignalCorrection             - Inter signal correction terms
    %   ReferenceTimeCEIPropagation       - Reference time of CEI
    %                                       propagation
    %   ReferenceWeekNumberCEIPropagation - Reference week number of CEI
    %                                       propagation
    %   URANEDID                          - Non-elevation dependent user
    %                                       range accuracy indices
    %   AlertFlag                         - Alert flag
    %
    %   References:
    %    [1] IS-GPS-200L. "NAVSTAR GPS Space Segment/Navigation User
    %        Segment Interfaces." GPS Enterprise Space & Missile Systems
    %        Center (SMC) - LAAFB, May 14, 2020.
    %
    %   See also HelperGPSNavigationConfig, HelperGPSNAVDataEncode.

    %   Copyright 2021 The MathWorks, Inc.

    properties
        %SVHealth Satellite vehicle health
        %  Indicate the satellite health as an integer value. This property
        %  is valid only when DataType is set to "LNAV". The default is 0.
        SVHealth = 0
        %SignalHealth L1/L2/L5 signal health indicators.
        %  Indicate the signal health of L1/L2/L5 as a three element array
        %  of binary values. This property is valid only when DataType is
        %  set to "CNAV". The default is a 3 element column vector of
        %  zeros.
        SignalHealth = [0; 0; 0]
        %IssueOfDataClock Issue of data clock
        %  Indicate the issue of data clock as an integer value. In the
        %  encoded message, this value is going to be of 10 bits. This
        %  property is valid only when DataType is set to "LNAV". The
        %  default is 0.
        IssueOfDataClock = 0
        %URAID User range accuracy index
        %  Indicate the user range accuracy as an integer for LNAV data.
        %  This property is valid only when DataType is set to "LNAV". The
        %  default is 0.
        URAID = 0
        %WeekNumber GPS week number
        %  Indicate the GPS week number as an integer value. The default
        %  is 2149.
        WeekNumber = 2149
        %GroupDelayDifferential Group delay differential
        %  Indicate the group delay differential value in seconds. The
        %  default is 0.
        GroupDelayDifferential = 0 % T_GD
        %SVClockCorrectionCoefficients Clock correction coefficients, af0,
        %af1, af2
        %  Indicate the satellite vehicle (SV) clock bias (af0), clock
        %  drift (af1) and clock drift rate (af2) in that order as an array
        %  of three elements. The default is a 3 element column vector with
        %  all zeros.
        SVClockCorrectionCoefficients = [0; 0; 0] % [af0; af1; af2]
        %ReferenceTimeOfClock Reference time of clock
        %  Indicate the reference time of clock in seconds. The default is
        %  0.
        ReferenceTimeOfClock = 0 % t_oc
        %SemiMajorAxisLength Semi-major axis length of the satellite orbit
        %  Indicate the semi-major axis length of the satellite orbit as a
        %  scalar double value in kilo-meters. The default is 26560.
        SemiMajorAxisLength = 26560
        %ChangeRateInSemiMajorAxis Rate of change of semi-major axis length
        %  Indicate the rate of change of semi-major axis length in meters
        %  per second. This property is valid only when DataType is set to
        %  "CNAV". The default is 0.
        ChangeRateInSemiMajorAxis = 0
        %MeanMotionDifference Mean motion difference
        %  Indicate the mean motion difference value as a scalar double.
        %  The default is 0.
        MeanMotionDifference = 0
        %RateOfMeanMotionDifference Rate of change of mean motion
        %difference
        %  Indicate the rate of change of mean motion difference as a
        %  scalar double value. This property is valid only when DataType
        %  is set to "CNAV". The default is 0.
        RateOfMeanMotionDifference = 0
        %FitIntervalFlag Fit interval flag
        %  Indicate the fit interval flag as a binary value. This property
        %  is valid only when DataType is set to "LNAV". The default is 0.
        FitIntervalFlag = 0
        %Eccentricity Eccentricity of the satellite orbit
        %  Indicate the eccentricity of the ellipse in which satellite
        %  orbits as a scalar double value in the range of 0 to 1. The
        %  default is 0.02.
        Eccentricity = 0.02
        %MeanAnomaly Mean anomaly at reference time
        %  Indicate the mean anomaly value as a scalar double value. The
        %  default is 0.
        MeanAnomaly = 0
        %ReferenceTimeOfEphemeris Reference time of ephemeris
        %  Indicate the reference time of epehemeris as a scalar double
        %  value. This value indicates the time witnin a week when the
        %  ephemeris data is updated in seconds. The default is 0.
        ReferenceTimeOfEphemeris = 0 % t_oe
        %HarmonicCorrectionTerms Six harmonic correction terms
        %  Indicate the six harmonic correction terms as a vector of 6
        %  elements. First element is the amplitude of the sine harmonic
        %  correction term to the angle of inclination (C_is). Second
        %  element is the amplitude of the cosine harmonic correction term
        %  to the angle of inclination (C_ic). The third element is the
        %  amplitude of the sine correction term to the orbit radius
        %  (C_rs). The fourth element is the Amplitude of the cosine
        %  correction term to the orbit radius (C_rc). The fifth element is
        %  the amplitude of the sine harmonic correction term to the
        %  argument of latitude (C_us). The sixth element is the amplitude
        %  of the cosine harmonic correction term to the argument of
        %  latitude (C_uc). The default is a column vector of six zeros.
        HarmonicCorrectionTerms = zeros(6,1) % [Cis; Cic; Crs; Crc; Cus; Cuc]
        %IssueOfDataEphemeris Issue of data ephemeris
        %  Indicate the issue of data ephemeris as a scalar integer value.
        %  This property is valid only when DataType is set to "LNAV". The
        %  default is 0.
        IssueOfDataEphemeris = 0
        %IntegrityStatusFlag Integrity status flag
        %  Indicate the signal integrity status as a binary scalar value.
        %  The default is 0.
        IntegrityStatusFlag = 0
        %ArgumentOfPerigee Argument of perigee at reference time
        %  Indicate the argument of perigee of the satellite orbit as a
        %  scalar double value in the units of semi-circles. Argument of
        %  perigee is defined as the angle subtended by the direction of
        %  longitude of ascending node to the perigee. The default is
        %  -0.52.
        ArgumentOfPerigee = -0.52
        %RateOfRightAscension Rate of right ascension
        %  Indicate the rate of change of right ascension as a scalar
        %  double value. The default is 0.
        RateOfRightAscension = 0
        %LongitudeOfAscendingNode Longitude of ascending node
        %  Indicate the longitude of ascending node as a scalar double
        %  value. The default is -0.84.
        LongitudeOfAscendingNode = -0.84
        %Inclination Inclination angle of the satellite orbit with respect
        %to equator of Earth
        %  Indicate the inclinatin angle as a scalar double value in the
        %  units of semi-circles. The default is 0.3.
        Inclination = 0.3 % In semi-circles
        %InclinationRate Rate of change of inclination angle
        %  Indicate the rate of change of inclination angle as a scalar
        %  double in the units of semi-circles/second. The default is 0.
        InclinationRate = 0
        %URAEDID Elevation dependent (ED) user range accuracy (URA) index
        %   Indicate the elevation dependent user range accuracy index as
        %   an integer. This property is valid only when DataType is set to
        %   "CNAV". The default is 0.
        URAEDID = 0
        %InterSignalCorrection Inter signal correction terms
        %  Indicate the inter signal correction (ISC) terms as a vector of
        %  4 elements. First element represents ISC L1C/A. Second element
        %  represents ISC L2C. Third element represents ISC L5I5. Fourth
        %  element represents ISC L5Q5. This property is valid only when
        %  DataType is set to "CNAV". The default is a column vector of 4
        %  zeros.
        InterSignalCorrection = zeros(4,1) % [L1C/A; L2C; L5I5; L5Q5]
        %ReferenceTimeCEIPropagation Reference time of CEI propagation
        %  Indicate the reference time of CEI propagation as a scalar
        %  double value. This property is valid only when DataType is set
        %  to "CNAV". The default value is 0.
        ReferenceTimeCEIPropagation = 0 % t_op
        %ReferenceWeekNumberCEIPropagation Reference week number of CEI
        %propagation
        %  Indicate the reference week number of clock, ephemeris, and
        %  integrity (CEI) paramters propagation. This property is valid
        %  only when DataType is set to "CNAV". The default is 101.
        ReferenceWeekNumberCEIPropagation = 101 % WN_OP
        %URANEDID Non-elevation dependent user range accuracy indices
        %  Indicate the Non-elevation dependent user range accuracy indices
        %  as a vector of three elements. The default is a column vector of
        %  3 zeros.
        URANEDID = [0; 0 ; 0] % [URA_NED0; URA_NED1; URA_NED2]
        %AlertFlag Alert flag
        %  Indicate the alert flag as a binary scalar value. The default is
        %  0.
        AlertFlag = 0
    end

    properties(Hidden)
        %DataType Data type of the navigation data
        %  Indicate the datatype as one of "CNAV" | "LNAV" | "NAV". The
        %  default is "NAV".
        DataType = "NAV"
    end

    methods
        function obj = HelperGPSCEIConfig(varargin)
            %HelperGPSNavigationParameters Construct an instance of this class
            %   Support name-value pair arguments when constructing object.
            obj@comm.internal.ConfigBase(varargin{:});
        end

        % Set methods
        function obj = set.SVHealth(obj,val)
            prop ='SVHealth';
            validateattributes(val,{'double','single','uint8'},{'scalar','nonnegative','integer'},mfilename,prop)
            obj.(prop) = val;
        end

        function obj = set.SignalHealth(obj,val)
            prop ='SignalHealth';
            validateattributes(val,{'double','single','uint8'},{'vector','nonnegative','integer','numel',3},mfilename,prop)
            obj.(prop) = val;
        end
    end

    methods(Access = protected)
        function flag = isInactiveProperty(obj,prop)
            flag = false;
            if any(strcmp(prop,{'SignalHealth','ChangeRateInSemiMajorAxis','RateOfMeanMotionDifference','URAEDID','InterSignalCorrection','ReferenceTimeCEIPropagation','ReferenceWeekNumberCEIPropagation','URANEDID'}))
                flag = ~any(strcmp(obj.DataType,{'CNAV','NAV'}));
            elseif any(strcmp(prop,{'SVHealth','IssueOfDataClock','URAID','FitIntervalFlag','IssueOfDataEphemeris'}))
                flag = ~any(strcmp(obj.DataType,{'LNAV','NAV'}));
            end
        end
    end
end