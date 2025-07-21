function varargout = HelperGPSConvertTime(varargin)
%HelperGPSConvertTime Convert GPS week and time of week into datetime
%object and vice-versa
%
%   Note: This is a helper function and its API and/or functionality may
%   change in subsequent releases.
%
%   T = HelperGPSConvertTime(W,TOW) convert GPS week number, W and time of
%   week, TOW into a datetime object, T.
%
%   [W,TOW] = HelperGPSConvertTime(T) convert datetime object, T into GPS
%   week number, W and time of week, TOW.
%
%   This function chooses which way to convert based on the data type of
%   first input argument. If the first input argument is of datetime
%   object, then there must be only one input and the function returns two
%   outputs - GPS week number and time of week. Whereas if the first input
%   argument is not datetime object, then two inputs (GPS week number and
%   time of week) must be provided and only one output will be available.

%   Copyright 2022 The MathWorks, Inc.

narginchk(1,2);
gpsStart = datetime(1980, 1, 6, 0, 0, 0, 'TimeZone', 'UTCLeapSeconds');
secondsPerWeek = 604800;
if isa(varargin{1},'datetime')
    % Then number of input arguments is 1 and number of output arguments is
    % 2. This function converts the datetime object into appropriate week
    % number and time in that week

    t = varargin{1};
    if isempty(t.TimeZone)
        t.TimeZone = 'UTC';
    end
    t.TimeZone = 'UTCLeapSeconds';

    % Time from GPS start time.
    dateDiff = seconds(t - gpsStart);
    gpsWeek = floor(dateDiff / secondsPerWeek);
    tow = dateDiff - (gpsWeek * secondsPerWeek);

    varargout{1} = gpsWeek;
    varargout{2} = tow;
else
    % Then number of input arguments is 2 and number of output arguments is
    % 1. Inputs GPS week number and time of week are being converted into
    % datetime object.

    varargout{1} = matlabshared.internal.gnss.GPSTime.getLocalTime(varargin{1},varargin{2});

    gpsWeek = varargin{1};
    tow = varargin{2};
    t = gpsStart + seconds((gpsWeek * secondsPerWeek) + tow);
    t.TimeZone = 'UTC';
end

end