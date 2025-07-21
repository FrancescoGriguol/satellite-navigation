function [cfg,parityChecks] = HelperGPSLNAVDataDecode(bits,cfg,varargin)
%HelperGPSLNAVDataDecode Decode the GPS LNAV data
%
%   Note: This is a helper and its API and/or functionality may change
%   in subsequent releases.
%   [CFG,PARITYCHECKS] = HelperGPSLNAVDataDecode(BITS,CFG) decodes the BITS
%   and assigns the parameter values in the structure CFG. BITS is a column
%   vector of 300 bits corresponding to one subframe.
%
%   [CFG,PARITYCHECKS] = ...
%   HelperGPSLNAVDataDecode(BITS,CFG,DECODEEVENIFPARITYFAIL) decodes the
%   BITS and returns the decoded parameters in CFG even if parity check is
%   failed when DECODEEVENIFPARITYFAIL is set to true. The default value of
%   DECODEEVENIFPARITYFAIL is false. That is by default when parity check
%   of a word in a subframe is failed, then those parameters are not
%   updated in the output CFG.
%
%   Inputs: 
%    BITS                   - Input bits. Column vector with length must be
%                             equal to 300 as number of bits in a subframe
%                             is 300.
%    CFG                    - A structure that may contain previously
%                             decoded fields or an empty structure.
%    DECODEEVENIFPARITYFAIL - (Optional) Set to true to decode even if
%                             parity check fails.
%
%   Outputs:
%    CFG                    - A structure containing parameters that are
%                             decoded from BITS. Potential properties are
%                             same as that of HelperGPSNavigationConfig
%    PARITYCHECKS           - A vector of 10 bits telling whether the
%                             parity check of each word in a subframe
%                             passed or not.

%   Copyright 2021-2023 The MathWorks, Inc.

if nargin > 2
    decodeEvenIfParityFail = varargin{1};
else
    decodeEvenIfParityFail = false;
end

% Data type is LNAV always
cfg.SignalType = "LNAV";

% Convert the subframe into words of 30 bits each
words = reshape(bits,30,10);
parityChecks = zeros(1,10);

%% Extract word1
wordNumber = 1;
[word1,parityBits,parityChecks(1)] = HelperGPSLNAVWordDecode(words(:,1),wordNumber,[0;0]);

if parityChecks(1) || decodeEvenIfParityFail
    % Identify the preamble. While this information is useless for the
    % receiver, this will help to compare against transmitted value
    cfg.Preamble = bit2int(word1(1:8,1),8);

    % Extract the TLM message
    cfg.TLMMessage = bit2int(word1(9:22),14);

    % Extract integrity status flag
    cfg.IntegrityStatusFlag = word1(23);
end

%% Extract word2
[word2,~,parityChecks(2)] = HelperGPSLNAVWordDecode(words(:,2),wordNumber,parityBits(end-1:end,1));
subframeID = bit2int(word2(20:22),3);

if parityChecks(2) || decodeEvenIfParityFail
    cfg.HOWTOW = bit2int(word2(1:17),17);
    cfg.AlertFlag = word2(18);
    cfg.AntiSpoofFlag = word2(19);
    cfg.SubframeID = subframeID;
end

% Based on the detected subframe ID, process other 8 words
switch(subframeID)
    case 1
        [cfg,parityChecks(3:10)] = decodeSubframe1(words(:,3:end),cfg,decodeEvenIfParityFail);
    case 2
        [cfg,parityChecks(3:10)] = decodeSubframe2(words(:,3:end),cfg,decodeEvenIfParityFail);
    case 3
        [cfg,parityChecks(3:10)] = decodeSubframe3(words(:,3:end),cfg,decodeEvenIfParityFail);
    case 4
        [cfg,parityChecks(3:10)] = decodeSubframe4(words(:,3:end),cfg,decodeEvenIfParityFail);
    case 5
        [cfg,parityChecks(3:10)] = decodeSubframe5(words(:,3:end),cfg,decodeEvenIfParityFail);
    otherwise
        warning("Detected invalid subframe ID");
end

% Arrange the field names of the structure in an order
% Arrange the field names of CEI data set
CEIDataEleorder = ["SVHealth";"SignalHealth";"IssueOfDataClock";"URAID";"WeekNumber"; ...
    "GroupDelayDifferential";"SVClockCorrectionCoefficients";"ReferenceTimeOfClock"; ...
    "SemiMajorAxisLength";"ChangeRateInSemiMajorAxis";"MeanMotionDifference"; ...
    "RateOfMeanMotionDifference";"FitIntervalFlag";"Eccentricity";"MeanAnomaly"; ...
    "ReferenceTimeOfEphemeris";"HarmonicCorrectionTerms";"IssueOfDataEphemeris"; ...
    "IntegrityStatusFlag";"ArgumentOfPerigee";"RateOfRightAscension"; ...
    "LongitudeOfAscendingNode";"Inclination";"InclinationRate";"URAEDID"; ...
    "InterSignalCorrection";"ReferenceTimeCEIPropagation";"ReferenceWeekNumberCEIPropagation"];
eleorder = ["SignalType";"PRNID";"MessageTypes";"FrameIndices";"Preamble"; ...
    "TLMMessage";"HOWTOW";"AntiSpoofFlag";"CodesOnL2";"L2PDataFlag"; ...
    "L2CPhasing";CEIDataEleorder;"AgeOfDataOffset";"NMCTAvailabilityIndicator"; ...
    "NMCTERD";"AlmanacFileName";"Ionosphere";"EarthOrientation";"UTC"; ...
    "DifferentialCorrection";"TimeOffset";"TextMessage";"TextInMessageType36";"TextInMessageType15"];
cfg = orderStructFields(cfg,eleorder);
end

function [cfg,parityChecks] = decodeSubframe1(words,cfg,decodeEvenIfParityFail)

% Decode the received words
numWords = size(words,2);
decwords = zeros(24,numWords); % 24 is the number of data bits in each word
parityChecks = zeros(numWords,1);
parityBits = zeros(6,1);
for iWord = 1:numWords
    [decwords(:,iWord),parityBits,parityChecks(iWord)] = ...
        HelperGPSLNAVWordDecode(words(:,iWord),iWord+2,parityBits(end-1:end,1));
end

% Decoding Word 3 and first bit of word 4
if parityChecks(1) || decodeEvenIfParityFail % Decode only if parity check passes
    cfg.WeekNumber = bit2int(decwords(1:10,1),10);
    codesOnL2 = bit2int(decwords(11:12,1),2);
    switch(codesOnL2)
        case 1
            cfg.CodesOnL2 = "P-code";
        case 2
            cfg.CodesOnL2 = "C/A-code";
        otherwise
            cfg.CodesOnL2 = "Invalid";
    end
    cfg.URAID = bit2int(decwords(13:16,1),4);
    cfg.SVHealth = bit2int(decwords(17:22,1),6);
    if parityChecks(6) % Decode only if parity check passes
        cfg.IssueOfDataClock = bit2int([decwords(23:24,1);decwords(1:8,6)],10);
    end
end

if parityChecks(2) || decodeEvenIfParityFail % Decode only if parity check passes
    cfg.L2PDataFlag = decwords(1,2);
end
% Word4, Word 5 and Word 6 are all reserved except for 1st bit of word 4.

if parityChecks(5) || decodeEvenIfParityFail
    % Decoding word 7. First 16 bits of word 7 are reserved
    cfg.GroupDelayDifferential = ssbit2int(decwords(17:24,5),8,true,2^(-31));
end

if parityChecks(6) || decodeEvenIfParityFail
    % Decoding word8. First 8 bits contains IODC and that is already decoded as
    % part of word 3
    cfg.ReferenceTimeOfClock = ssbit2int(decwords(9:24,6),16,false,2^(4));
end

if parityChecks(7) || decodeEvenIfParityFail
    % Decoding words 9 and 10 for af0,f1,af2
    af2 = ssbit2int(decwords(1:8,7),8,true,2^(-55));
    af1 = ssbit2int(decwords(9:24,7),16,true,2^(-43));
    if isfield(cfg,"SVClockCorrectionCoefficients")
        af0 = cfg.SVClockCorrectionCoefficients(1);
    else
        af0 = 0;
    end
    cfg.SVClockCorrectionCoefficients = [af0;af1;af2];
    if parityChecks(8) || decodeEvenIfParityFail
        af0 = ssbit2int(decwords(1:22,8),22,true,2^(-31));
        cfg.SVClockCorrectionCoefficients(1) = af0;
    end
end
end

function [cfg,parityChecks] = decodeSubframe2(words,cfg,decodeEvenIfParityFail)

% Decode the received words
numWords = size(words,2);
decwords = zeros(24,numWords); % 24 is the number of data bits in each word
parityChecks = zeros(numWords,1);
parityBits = zeros(6,1);
for iWord = 1:numWords
    [decwords(:,iWord),parityBits,parityChecks(iWord)] = ...
        HelperGPSLNAVWordDecode(words(:,iWord),iWord+2,parityBits(end-1:end,1));
end

dontupdateflag = false;
if parityChecks(1) || decodeEvenIfParityFail
    % Decode IODE
    cfg.IssueOfDataEphemeris = bit2int(decwords(1:8,1),8);

    % Decode C_rs
    Crs = ssbit2int(decwords(9:24,1),16,true,2^(-5));
else
    dontupdateflag = true;
end

if parityChecks(2) || decodeEvenIfParityFail
    % Decode MeanMotionDifference
    cfg.MeanMotionDifference = ssbit2int(decwords(1:16,2),16,true,2^(-43));

    if parityChecks(3) || decodeEvenIfParityFail
        % Decode MeanAnomaly
        cfg.MeanAnomaly = ssbit2int([decwords(17:24,2);decwords(1:24,3)],32,true,2^(-31));
    end
end

if parityChecks(4) || decodeEvenIfParityFail
    % Decode C_uc
    Cuc = ssbit2int(decwords(1:16,4),16,true,2^(-29));

    if parityChecks(5) || decodeEvenIfParityFail
        % Decode eccentricity
        cfg.Eccentricity = ssbit2int([decwords(17:24,4);decwords(1:24,5)],32,false,2^(-33));
    end
else
    dontupdateflag = true;
end

if parityChecks(6) || decodeEvenIfParityFail
    % Decode C_us
    Cus = ssbit2int(decwords(1:16,6),16,true,2^(-29));
else
    dontupdateflag = true;
end

% Harmonic correction terms are in the order of [Cis; Cic; Crs; Crc; Cus; Cuc]
if ~dontupdateflag
    if isprop(cfg,"HarmonicCorrectionTerms") || isfield(cfg,"HarmonicCorrectionTerms") % Condition supports both config object and a structure
        cfg.HarmonicCorrectionTerms(3) = Crs;
        cfg.HarmonicCorrectionTerms(5) = Cus;
        cfg.HarmonicCorrectionTerms(6) = Cuc;
    else
        cfg.HarmonicCorrectionTerms = [0;0;Crs;0;Cus;Cuc];
    end
end

if parityChecks(7) || decodeEvenIfParityFail
    % Decode SemiMajorAxisLength
    sqrtA = ssbit2int([decwords(17:24,6);decwords(1:24,7)],32,false,2^(-19));
    cfg.SemiMajorAxisLength = sqrtA^2;
end

if parityChecks(8) || decodeEvenIfParityFail
    % Decode ReferenceTimeOfEphemeris
    cfg.ReferenceTimeOfEphemeris = ssbit2int(decwords(1:16,8),16,false,2^(4));

    % Decode FitIntervalFlag
    cfg.FitIntervalFlag = decwords(17,8);

    % Decode AgeOfDataOffset (AODO)
    cfg.AgeOfDataOffset = ssbit2int(decwords(18:22,8),5,false,900);
end
end

function [cfg,parityChecks] = decodeSubframe3(words,cfg,decodeEvenIfParityFail)

% Decode the received words
numWords = size(words,2);
decwords = zeros(24,numWords); % 24 is the number of data bits in each word
parityChecks = zeros(numWords,1);
parityBits = zeros(6,1);
for iWord = 1:numWords
    [decwords(:,iWord),parityBits,parityChecks(iWord)] = ...
        HelperGPSLNAVWordDecode(words(:,iWord),iWord+2,parityBits(end-1:end,1));
end

dontupdateflag = false;
if parityChecks(1) || decodeEvenIfParityFail
    % Decode C_ic
    Cic = ssbit2int(decwords(1:16,1),16,true,2^(-29));

    if parityChecks(2) || decodeEvenIfParityFail
        % Decode eccentricity
        cfg.LongitudeOfAscendingNode = ssbit2int([decwords(17:24,1);decwords(1:24,2)],32,true,2^(-31));
    end
else
    dontupdateflag = true;
end

if parityChecks(3) || decodeEvenIfParityFail
    % Decode C_is
    Cis = ssbit2int(decwords(1:16,3),16,true,2^(-29));

    if parityChecks(4) || decodeEvenIfParityFail
        % Decode Inclination
        cfg.Inclination = ssbit2int([decwords(17:24,3);decwords(1:24,4)],32,true,2^(-31));
    end
else
    dontupdateflag = true;
end

if parityChecks(5) || decodeEvenIfParityFail
    % Decode C_rc
    Crc = ssbit2int(decwords(1:16,5),16,true,2^(-5));

    if parityChecks(6) || decodeEvenIfParityFail
        % Decode ArgumentOfPerigee
        cfg.ArgumentOfPerigee = ssbit2int([decwords(17:24,5);decwords(1:24,6)],32,true,2^(-31));
    end
else
    dontupdateflag = true;
end

if ~dontupdateflag
    % Harmonic correction terms are in the order of [Cis; Cic; Crs; Crc; Cus; Cuc]
    if isprop(cfg,"HarmonicCorrectionTerms") || isfield(cfg,"HarmonicCorrectionTerms") % Condition supports both config object and a structure
        cfg.HarmonicCorrectionTerms(1) = Cis;
        cfg.HarmonicCorrectionTerms(2) = Cic;
        cfg.HarmonicCorrectionTerms(4) = Crc;
    else
        cfg.HarmonicCorrectionTerms = [Cis;Cic;0;Crc;0;0];
    end
end


if parityChecks(7) || decodeEvenIfParityFail
    % Decode RateOfRightAscension
    cfg.RateOfRightAscension = ssbit2int(decwords(1:24,7),24,true,2^(-43));
end

if parityChecks(8) || decodeEvenIfParityFail
    % Decode IODE
    cfg.IssueOfDataEphemeris = bit2int(decwords(1:8,8),8);

    % Decode InclinationRate (IDOT)
    cfg.InclinationRate = ssbit2int(decwords(9:22,8),14,true,2^(-43));
end
end

function [cfg,parityChecks] = decodeSubframe4(words,cfg,decodeEvenIfParityFail)

% Decode the received words
numWords = size(words,2);
decwords = zeros(24,numWords); % 24 is the number of data bits in each word
parityChecks = zeros(numWords,1);
parityBits = zeros(6,1);
for iWord = 1:numWords
    [decwords(:,iWord),parityBits,parityChecks(iWord)] = ... 
        HelperGPSLNAVWordDecode(words(:,iWord),iWord+2,parityBits(end-1:end,1));
end

% Decode Data ID and SV ID
dataID = decwords(1:2,1);
svID = bit2int(decwords(3:8,1),6);
if parityChecks(1) || decodeEvenIfParityFail
    cfg.DataID = dataID;
    cfg.SVPageID = svID;
end


switch(svID)
    case {[25:28,29:32]}
        % Almanac data in pages 2 to 5 and 7 to 10 in subframe 4
        cfg.Almanac(svID) = almdata(decwords);
    case 52
        % Page 13 in subframe 4
        if parityChecks(1) || decodeEvenIfParityFail
            cfg.NMCTAvailabilityIndicator = bit2int(decwords(9:10,1),2);
        end
        temp1 = decwords(:);
        nmctBin = reshape(temp1(11:end-2), 6,30);
        erd = zeros(30,1);
        numBits = 6;
        scaleFactor = 0.3;
        for iERD = 1:size(nmctBin,2)
            erd(iERD) = ssbit2int(nmctBin(:,iERD),numBits,true,scaleFactor);
        end
        cfg.NMCTERD = erd;
    case 55
        % Page 17 in subframe 4
        temp1 = decwords(:);
        cfg.TextMessage = char(bit2int(temp1(9:end-8),8).');
    case 56
        % Page 18 in subframe 4
        % Ionosphere parameters
        temp1 = decwords(:);
        alphaparamsbin = reshape(temp1(9:40),8,[]);
        betaparamsbin = reshape(temp1(41:72),8,[]);
        alphascales = 2.^[-30;-27;-24;-24];
        betascales = 2.^[11;14;16;16];
        alphaparams = zeros(4,1);
        betaparams = zeros(4,1);
        for iParam = 1:4
            alphaparams(iParam) = ssbit2int(alphaparamsbin(:,iParam),8,true,alphascales(iParam));
            betaparams(iParam) = ssbit2int(betaparamsbin(:,iParam),8,true,betascales(iParam));
        end
        cfg.Ionosphere.Alpha = alphaparams;
        cfg.Ionosphere.Beta = betaparams;

        % UTC parameters
        cfg.UTC.UTCTimeCoefficients(2) = ssbit2int(decwords(:,4),24,true,2^(-50));
        cfg.UTC.UTCTimeCoefficients(1) = ssbit2int([decwords(:,5);decwords(1:8,6)],32,true,2^(-30));
        cfg.ReferenceTimeUTCData = ssbit2int(decwords(9:16,6),8,false,2^(12));
        cfg.TimeDataReferenceWeekNumber = bit2int(decwords(17:24,6),8);
        cfg.PastLeapSecondCount = ssbit2int(decwords(1:8,7),8,true,1);
        cfg.LeapSecondReferenceWeekNumber = ssbit2int(decwords(9:16,7),8,false,1);
        cfg.LeapSecondReferenceDayNumber = bit2int(decwords(17:24,7),8);
        cfg.FutureLeapSecondCount = ssbit2int(decwords(1:8,8),8,true,1);
    case 63
        % Page 25 in subframe 4
        % Decode SV config
        temp1 = decwords(:);
        svconfigbits = reshape(temp1(8+(1:128)),4,[]);
        for iSV = 1:size(svconfigbits,2)
            cfg.Almanac.Data(iSV).SVConfiguration = bit2int(svconfigbits(:,iSV),4);
        end

        % Decode SV health values from page 25 of subframe 4
        svhealthbits1 = decwords(19:24,6);
        svhealthbits2 = reshape(decwords(:,7),6,[]);
        svhealthbits3 = reshape(decwords(1:18,8),6,[]);
        svhealthbits = [svhealthbits1,svhealthbits2,svhealthbits3];
        for iSV = 1:size(svhealthbits,2)
            cfg.AlmanacSVHealth(iSV+24) = bit2int(svhealthbits(:,iSV),6); % SV health of satellite SV 25 to 32. So, add 24 to the index.
        end
    otherwise

end
end

function [cfg,parityChecks] = decodeSubframe5(words,cfg,decodeEvenIfParityFail)

% Decode the received words
numWords = size(words,2);
decwords = zeros(24,numWords); % 24 is the number of data bits in each word
parityChecks = zeros(numWords,1);
parityBits = zeros(6,1);
for iWord = 1:numWords
    [decwords(:,iWord),parityBits,parityChecks(iWord)] = ... 
        HelperGPSLNAVWordDecode(words(:,iWord),iWord+2,parityBits(end-1:end,1));
end

% Decode Data ID and SV ID
dataID = decwords(1:2,1);
svID = bit2int(decwords(3:8,1),6);
if parityChecks(1) || decodeEvenIfParityFail
    cfg.DataID = dataID;
    cfg.SVPageID = svID;
end

switch(svID)
    case num2cell(0:24)
        % If SVPageID is zero, that indicates, that particular satellite is
        % not in the constellation and decoding is not necessary.
        if cfg.SVPageID % Non-zero condition
            % Almanac data
            [cfg.Almanac.Data(cfg.SVPageID),cfg.Almanac.ReferenceTimeOfAlmanac] = almdata(decwords);
        end
    case 51
        % Decode page 25 in subframe 5
        % Decode the SV health data of all the SVs
        cfg.Almanac.ReferenceTimeOfAlmanac = ssbit2int(decwords(9:16,1),8,false,2^(12)); 
        cfg.Almanac.AlmanacWeekNumber = bit2int(decwords(17:24,1),8);
        svhealthbits = reshape(decwords(:,2:7),6,[]);
        for iSV = 1:size(svhealthbits,2)
            cfg.AlmanacSVHealth(iSV) = bit2int(svhealthbits(:,iSV),6);
        end
    otherwise
        warning("Detected invalid SV (page) ID in subframe 5.");
end
end

function [almstruct,ReferenceTimeOfAlmanac] = almdata(decwords)

parityChecks = ones(10,1);
decodeEvenIfParityFail = 1;
if parityChecks(1) || decodeEvenIfParityFail
    almstruct.PRNID = bit2int(decwords(3:8,1),6);
    almstruct.Eccentricity = ssbit2int(decwords(9:24,1),16,false,2^(-21));
end

ReferenceTimeOfAlmanac = ssbit2int(decwords(1:8,2),8,false,2^(12));
if parityChecks(2) || decodeEvenIfParityFail
    almstruct.InclinationDifference = ssbit2int(decwords(9:24,2),16,true,2^(-19));
end

if parityChecks(3) || decodeEvenIfParityFail
    almstruct.RateOfRightAscension = ssbit2int(decwords(1:16,3),16,true,2^(-38));
    almstruct.SVHealth = bit2int(decwords(17:24,3),8);
end

if parityChecks(4) || decodeEvenIfParityFail
    almstruct.RootOfSemiMajorAxis = ssbit2int(decwords(1:24,4),24,false,2^(-11));
end

if parityChecks(5) || decodeEvenIfParityFail
    almstruct.LongitudeOfAscendingNode = ssbit2int(decwords(1:24,5),24,true,2^(-23));
end

if parityChecks(6) || decodeEvenIfParityFail
    almstruct.ArgumentOfPerigee = ssbit2int(decwords(1:24,6),24,true,2^(-23));
end

if parityChecks(7) || decodeEvenIfParityFail
    almstruct.MeanAnomaly = ssbit2int(decwords(1:24,7),24,true,2^(-23));
end

if parityChecks(8) || decodeEvenIfParityFail
    af0 = ssbit2int([decwords(1:8,8);decwords(20:22,8)],11,true,2^(-20));
    af1 = ssbit2int(decwords(9:19,8),11,true,2^(-38));
    almstruct.SVClockCorrectionCoefficients = [af0,af1];
end
end

function val = ssbit2int(bits,numbits,issigned,scalefactor)
%SSBIT2INT Convert signed bit values to integer and scale it

intval = bit2int(bits,numbits,IsSigned=issigned);
val = intval*scalefactor;
end

function newStruct = orderStructFields(oldStruct,eleorder)
eleorder = [eleorder;string(fieldnames(oldStruct))];
for ifield = 1:length(eleorder)
    if isfield(oldStruct,eleorder(ifield))
        newStruct.(eleorder(ifield)) = oldStruct.(eleorder(ifield));
    end

end
end
