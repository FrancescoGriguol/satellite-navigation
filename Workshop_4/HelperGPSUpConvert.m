classdef HelperGPSUpConvert < comm.internal.Helper
    %HelperGPSUpConvert Apply sinusoidal doppler shift to input signal
    %
    %   Note: This is a helper function and its API and/or functionality
    %   may change in subsequent releases.
    %
    %   UC = HelperGPSUpConvert creates an up-conversion System object with
    %   sinusoidal varying frequency offset, UC. This object up-converts a
    %   baseband signal to higher frequencies.
    %
    %   Step method syntax:
    %
    %   [Y,PHASE] = step(UC,X) applies a sinusoidally varying frequency
    %   offset and upconverts the input baseband signal, X, and returns Y.
    %   The input X is a double precision column X, of dimensions M. M is
    %   the number of time samples in the input signals. The data type of X
    %   and Y are the same. Control the input signal sampling rate using
    %   UC.BasebandSampleRate and output IF signal sample rate using
    %   UC.IFSampleRate.
    %
    %   HelperGPSUpConvert properties (All the properties are non-tunable):
    %
    %   CenterFrequency       - Center frequency to which signal must be up
    %                           converted. Set to zero if the signal must
    %                           be in baseband (Default: 10e6 Hz)
    %   DopplerRate           - Rate of change of frequency offset in
    %                           Hz/sec (-10e3 to 10e3 Hz. Default:3e3 Hz)
    %   PeakDoppler           - Maximum doppler shift in Hz (2e3 to -2e3 Hz/sec. Default: 1e3 Hz/sec)
    %   BasebandSampleRate    - Sample rate of baseband signal (Default:
    %                           10.23e6 Hz)
    %   IFSampleRate          - Sample rate of IF signal (Default: 38.192e6 Hz)

    %   Copyright 2021-2022 The MathWorks, Inc.

    % Public, non-tunable properties
    properties(Nontunable)
        CenterFrequency = 10e6
        %DopplerRate Doppler rate in Hertz/sec
        % Specify the rate of change of doppler shift in Hertz/sec. The
        % doppler shift is modeled in a sinusoidal variation with a period
        % that results in doppler rate. This property is nontunable.
        DopplerRate = 1e3
        %PeakDoppler Peak doppler in Hertz
        %   Specify the maximum doppler shift in Hertz. The doppler shift
        %   will vary between +/- peak doppler. This property is
        %   nontunable.
        PeakDoppler = 3e3
        %BasebandSampleRate Sample rate in Hertz
        %   Specify the sample rate of the input samples in Hz as a double
        %   precision, real, positive scalar.  This property is nontunable.
        BasebandSampleRate = 10.23e6
        %IFSamplingRate
        IFSampleRate = 38.192e6

    end
 

    % Pre-computed constants
    properties(Access = private)
        pBBSampleIndex = 0
        pIFSampleIndex = 0
        pUpsampleFactor
        pDownSampleFactor
    end

    methods
        % Constructor
        function obj = HelperGPSUpConvert(varargin)
            % Support name-value pair arguments when constructing object
            setProperties(obj,nargin,varargin{:})
        end

    end

    methods(Access = protected)
        %% Common functions
        function setupImpl(obj)
            [obj.pUpsampleFactor,obj.pDownSampleFactor] = rat(obj.IFSampleRate/obj.BasebandSampleRate);
        end

        function [y,phase] = stepImpl(obj,u)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            numSamples = length(u);
            sampleIndices = obj.pBBSampleIndex + (0:(numSamples-1));
            ph = sin(obj.DopplerRate*sampleIndices/(obj.PeakDoppler*obj.BasebandSampleRate));
            phase = 2*pi*(obj.PeakDoppler^2)/obj.DopplerRate*ph;
            obj.pBBSampleIndex = obj.pBBSampleIndex + numSamples;
            % obj.pPhaseIntroducer.PhaseOffset = phases(:);
            bbwave = u(:).*exp(1j*phase(:));

            % Up-convert the signal to IF frequency
            % 1. Rate match the baseband signal to IF sampling frequency
            upwave = repelem(bbwave(:),obj.pUpsampleFactor); % Upsample
            rateMatchedData = upwave(1:obj.pDownSampleFactor:end); % Decimate
            % 2. Generate the IF carrier wave
            numSamples = length(rateMatchedData);
            sampleIndices = obj.pIFSampleIndex + (0:(numSamples-1));
            t = sampleIndices/obj.IFSampleRate;
            ph = 2*pi*obj.CenterFrequency*t;
            IFCarrier = exp(1j*ph);
            obj.pIFSampleIndex = obj.pIFSampleIndex + numSamples;
            % 3. Multiply the complex data and complex sinusoid and just
            %    take the real value to get the IF data
            IFComplexData = IFCarrier(:).*rateMatchedData(:);
            if obj.CenterFrequency ~= 0
                y = real(IFComplexData);
            else
                y = IFComplexData;
            end
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            obj.pBBSampleIndex = 0;
            obj.pIFSampleIndex = 0;
        end

        %% Backup/restore functions
        function s = saveObjectImpl(obj)
            % Set properties in structure s to values in object obj

            % Set public properties and states
            s = saveObjectImpl@matlab.System(obj);
            if isLocked(obj)
              s.pSampleIndex = obj.pBBSampleIndex;
            end
        end

        function loadObjectImpl(obj,s,wasLocked)
            % Set properties in object obj to values in structure s
            if wasLocked(obj)
               obj.pBBSampleIndex = s.pSampleIndex;
            end
            % Set public properties and states
            loadObjectImpl@matlab.System(obj,s,wasLocked);
        end
    end
end
