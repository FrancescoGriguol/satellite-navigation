classdef HelperGPSInitialSynchronization < matlab.System
    % HelperGPSInitialSynchronization Perform initial
    % synchronization on GPS signals
    %
    %   Note: This is a helper function and its API and/or functionality
    %   may change in subsequent releases.
    %
    %   INITSYNC = HelperGPSInitialSynchronization creates a
    %   synchronization object. This object estimates the available
    %   satellites in the incoming signal, coarse values of Doppler and
    %   delay in the signal of corresponding satellite.
    %
    %   Step method syntax:
    %
    %   [Y,CORRVAL,MAXCORRVALS] = step(INITSYNC,X) searches for what all
    %   satellite PRN IDs are available in the signal X. X must be a column
    %   vector of length equal to INITSYNC.SampleRate*1e-3. That will be
    %   equal to number of samples in one millisecond of digital input to
    %   the receiver. The reason to choose one millisecond is that C/A-code
    %   is of one millisecond duration. Y is a matrix with number of rows
    %   equal to INITSYNC.PRNIDsToSearch. Number of columns equal to 4.
    %   First column in Y represents the PRNID in question. 2nd column
    %   represents the coarse estimate of the Doppler in the corresponding
    %   row PRN ID. 3rd column represents the coarse delay in the units of
    %   number of C/A-code chips for the corresponding row PRN ID. 4th
    %   column is binary value indicating corresponding row PRN ID is
    %   detected or not (1 represents signal detected). CORRVAL is a
    %   3-dimensional matrix with each 2D matrix representing the
    %   correlation value of the searched satellites. Row in each layer of
    %   the 3D-matrix represents the correlation value for each delay
    %   value. Column in each layer of the 3D-matrix represents the
    %   correlation value for each frequency offset value. Layers are
    %   organized in decreasing order of maximum correlation value within
    %   each layer just like order of Y. MAXCORRVALS is a vector containing
    %   maximum correlation value in each layer of CORRVAL.
    %
    %   HelperGPSInitialSynchronization properties:
    %
    %   PRNIDsToSearch           - (Tunable) The satellite PRN IDs to
    %                              search while acquisition. Typically this
    %                              information is obtained after almanac
    %                              decoding, and for cold start search must
    %                              be done for all GPS satellites. This
    %                              must be a row vector. (Default: 1:32)
    %   CenterFrequency          - (Non-tunable) Center frequency of the
    %                              input signal. (Default: 0 Hz)
    %   SampleRate               - (Non-tunable) Sample rate of the input
    %                              signal (Default: 10.23e6 Hz)
    %   DetectionThresholdFactor - (Non-tunable) Detection threshold factor
    %                               (Default: 1.2)

    % References:
    %  [1] K. Borre, ed., A Software-Defined GPS and Galileo Receiver: A
    %      Single-Frequency Approach, Applied and Numerical Harmonic
    %      Analysis (Boston, Mass: Birkhäuser, 2007).

    %   Copyright 2021-2022 The MathWorks, Inc.

    % Public, tunable properties
    properties
        PRNIDsToSearch = 1:32 % Must be a row vector
    end

    properties(Nontunable)
        %CenterFrequency Center frequency of the input signal
        %   Specify the center frequency of the input signal in Hertz. If
        %   the signal is at baseband, CenterFrequency is 0. Default is 0.
        CenterFrequency = 0
        %SampleRate Sample rate of the input signal
        %   Specify the sample rate of the input signal in Hertz. Default
        %   is 10.23e6.
        SampleRate = 10.23e6
        %DetectionThresholdFactor Detection threshold factor
        %   Specify the detection threshold factor. Value of one means
        %   detection threshold is same as noise threshold. For proper
        %   detection, typically, this value is greater than one. The
        %   default value of 1.2 means that the threshold value is 1.2
        %   times the maximum noise power. The maximum signal correlation
        %   value or power is compared against this threshold value to
        %   detect the valid satellite signal.
        DetectionThresholdFactor = 1.2
    end

    properties(Constant,Hidden)
        pCAChipRate = 1.023e6 % This is constant for now
    end

    % Pre-computed constants
    properties(Access = private)
        % FFT size used for calculation of correlation
        pFFTSize
        % Frequency domain C/A-code chips used to calculate correlation
        pFqyCAChips
        % Complex exponential to save common computations
        pComplexExpo
        % Discrete frequency steps for which Doppler is searched
        pFrequencySteps
        % Frequency domain G1 code of the C/A-code
        pFqyConjG1Code
        % Number of C/A-code chips per one block of C/A-code after rate
        % matching
        pNumCAChips
    end

    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            fs = obj.SampleRate;
            [upfac,downfac] = rat(fs/obj.pCAChipRate);
            numCAChips = round(1023*upfac/downfac);
            obj.pNumCAChips = numCAChips;
            nPoints = numCAChips;
            obj.pFFTSize = nPoints;
            fqycachips = zeros(nPoints,32);
            numextrasampl = nPoints-numCAChips;

            % Generate G1 code of the C/A-code for generating noise
            % variance which is useful to calculate detection threshold for
            % detecting satellites

            % Generator G1 Polynomial g(X) = 1 + X^3 + X^10
            % pnseqG1 = comm.PNSequence("InitialConditions",ones(10,1),...
            %     "Polynomial","X^10 + X^7 + 1",...
            %     "SamplesPerFrame", 1023); % Generator polynomial in comm.PNSequence should be taken in reverse order than what is given in GPS standard
            % g1 = pnseqG1();
            g1 = ([1;1;1;1;1;1;1;1;1;1;0;0;0;1;1;1;0;0;0;1;0;0;1;1;1;0;1;1;0;0;1;0;1;...
                0;1;1;1;0;1;1;1;1;0;1;0;1;0;0;0;1;1;1;1;0;1;0;0;1;0;1;0;1;0;0;0;0;0;...
                1;0;1;1;1;1;1;1;1;1;0;1;0;1;0;1;0;1;0;1;1;1;1;0;1;0;0;0;0;1;1;1;0;1;...
                0;0;1;0;0;0;1;1;0;0;1;0;1;1;0;1;0;1;1;0;0;1;1;1;1;0;1;0;1;1;0;0;0;1;...
                1;0;0;1;1;1;1;1;1;0;0;1;0;1;0;1;0;1;0;0;1;1;0;0;1;1;0;0;1;0;1;0;0;1;...
                1;1;1;1;0;1;0;0;1;1;1;0;0;0;0;1;0;0;0;1;1;0;1;1;0;0;1;0;0;0;1;0;1;0;...
                0;1;1;0;1;1;1;1;0;1;1;1;0;1;0;1;0;1;1;1;0;0;1;1;0;0;1;1;1;0;1;1;1;0;...
                1;1;1;0;0;1;1;1;0;1;0;1;0;0;1;1;1;0;1;0;0;0;0;0;1;1;1;1;0;1;1;0;1;1;...
                1;0;0;0;0;1;1;0;0;0;1;0;0;1;0;1;0;0;1;0;1;1;0;0;1;1;0;1;0;0;0;1;0;0;...
                0;1;0;1;1;0;1;0;0;1;0;1;1;1;0;1;0;0;1;1;0;0;0;1;0;1;1;0;0;0;0;0;0;1;...
                0;1;0;0;1;0;0;1;0;1;1;1;1;1;0;1;1;1;1;0;0;0;1;1;0;0;0;1;1;0;1;1;1;0;...
                1;1;0;0;0;0;1;1;1;1;0;0;1;0;0;1;1;1;0;0;1;0;1;1;0;0;0;1;0;0;0;0;1;1;...
                0;1;1;1;1;1;1;1;0;0;1;1;1;0;0;0;1;1;0;1;0;1;0;0;1;0;1;0;0;0;0;1;0;0;...
                0;0;1;0;0;1;0;1;1;0;1;1;1;1;1;0;1;0;1;1;1;0;0;0;1;0;1;1;1;0;0;1;0;0;...
                0;0;1;1;1;1;1;0;1;1;0;1;0;1;0;1;0;0;0;1;0;1;1;1;1;0;1;1;0;0;1;1;1;0;...
                0;1;1;1;1;1;0;0;0;0;0;1;1;1;0;0;1;0;0;1;0;1;0;1;1;0;0;1;0;1;1;1;1;0;...
                0;1;0;1;1;1;0;0;0;0;0;1;0;1;0;1;1;0;1;1;0;0;1;1;0;0;0;0;1;1;0;1;0;1;...
                1;0;1;1;1;0;1;0;0;0;1;0;1;0;1;1;1;1;1;1;0;1;0;0;0;1;1;1;0;0;1;1;0;1;...
                1;1;0;0;1;0;1;0;0;0;1;1;0;1;0;0;0;0;0;0;1;1;0;0;1;0;0;1;0;0;0;1;0;0;...
                0;0;0;1;0;0;1;1;0;1;1;0;1;0;0;1;1;1;1;0;0;1;1;0;1;0;1;0;1;1;0;0;0;0;...
                1;0;1;1;1;0;1;1;0;1;0;0;0;1;1;0;0;0;0;1;0;0;1;1;1;1;1;1;1;0;1;1;1;0;...
                0;0;1;1;1;1;0;0;0;0;0;0;1;1;1;0;1;1;0;1;1;0;0;0;1;0;1;0;0;0;1;0;0;1;...
                1;0;0;1;0;0;0;0;0;1;1;0;1;0;0;1;0;0;1;1;1;1;0;1;1;1;1;1;0;0;0;1;0;1;...
                0;1;0;1;1;0;1;0;0;0;0;1;0;1;0;0;0;0;0;0;0;1;0;1;1;0;1;1;0;1;1;1;1;0;...
                0;1;1;1;1;0;0;0;1;0;0;0;1;1;1;1;1;1;0;1;1;0;0;0;1;1;1;0;1;0;1;1;0;1;...
                0;1;0;0;0;0;1;1;0;0;1;1;0;1;1;0;0;0;0;0;1;1;0;0;0;0;0;0;0;0;1;1;0;1;...
                1;0;1;1;0;1;0;1;1;1;0;1;0;1;1;1;1;0;0;0;0;1;0;1;0;1;0;0;1;0;0;0;0;1;...
                0;1;1;0;0;1;0;0;1;1;0;0;0;0;0;1;0;0;0;1;0;0;1;0;0;0;0;0;0;1;0;0;0;0;...
                0;0;0;0;0;1;0;0;1;0;0;1;0;0;1;1;0;1;0;0;1;1;0;1;0;1;1;1;1;1;0;0;1;1;...
                0;0;0;1;1;1;1;1;0;0;1;0;0;0;1;1;1;0;1;1;1;1;1;1;0;0;0;0;1;1;1;0;0;0;0;0;0;0]);
            upgcode = repelem(g1,upfac,1);
            rateMatchedgCode = upgcode(1:downfac:end);
            postfixGSamples = [rateMatchedgCode;rateMatchedgCode(1:numextrasampl)];
            obj.pFqyConjG1Code = conj(fft(postfixGSamples));

            for iCode = 1:32
                cacode = gnssCACode(iCode,"GPS");
                upcacode = repelem(cacode(:),upfac,1);
                rateMatchedCACode = upcacode(1:downfac:end);
                postfixsamples = [rateMatchedCACode;rateMatchedCACode(1:numextrasampl)];
                fqycachips(:,iCode) = fft(postfixsamples,nPoints);
            end
            obj.pFqyCAChips = conj(fqycachips);

            % Pre-initialize the complex exponential so that it is easy to
            % work in stepImpl
            maxFreqSearch = 10e3; % In Hz
            freqStep = 500; % In Hz
            allfqysteps = obj.CenterFrequency + (-1*maxFreqSearch:freqStep:maxFreqSearch);
            obj.pFrequencySteps = allfqysteps;
            numSamples = fs*1e-3; % Assuming only one millisecond of data per step call.
            tempPhases = (-2*pi*(0:numSamples-1)/fs).';
            phases = tempPhases.*allfqysteps; % By implicit expansion, this becomes a matrix with each row corresponding to phases of that particular frequency offset.
            obj.pComplexExpo = cos(phases)+1j*sin(phases);

        end

        function [y,corrval,maxcorrvals] = stepImpl(obj,u)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.

            % Output is a matrix with number of rows equal to
            % obj.PRNIDsToSearch amd number of columns equal to 4. First
            % column indicates the PRN ID, seconds column indicates the
            % Doppler frequency for that PRNID,3rd column indicates the
            % code phase offset and 4th column indicates whether signal is
            % detected or not in the corresponding PRNID.
            PRNIDs = obj.PRNIDsToSearch;
            numPRNIDsToSearch = length(PRNIDs);            

            cmplxexpo = obj.pComplexExpo;
            nPoints = obj.pFFTSize;
            fqycachips = obj.pFqyCAChips;
            allfqysteps = obj.pFrequencySteps;
            fc = obj.CenterFrequency;
            numCAChips = obj.pNumCAChips;
            th = obj.DetectionThresholdFactor;

            % Search for each satellite
            y = zeros(numPRNIDsToSearch,4);
            maxcorrvals = zeros(numPRNIDsToSearch,1);
            y(:,1) = PRNIDs(:);
            corrval = zeros(nPoints,length(allfqysteps),numPRNIDsToSearch);

            if numPRNIDsToSearch == 0 % Nothing to process if there are no satellites to search
                return
            end

            CarrierWipe = cmplxexpo.*u; % By implicit expansion, CarrierWipe becomes a matrix

            fsig = fft(CarrierWipe,nPoints); % Frequency domain signal after carrier wipe-off

            % Correlation peak with G1 to get maximum noise power
            noiseTime = ifft(fsig.*obj.pFqyConjG1Code);
            absNoise = (real(noiseTime).^2)+(imag(noiseTime).^2);
            maxNoise = max(absNoise(:));

            for iPRN = 1:numPRNIDsToSearch
                fcorval = fsig.*fqycachips(:,PRNIDs(iPRN)); % Frequency domain correlation value. Implicit expansion
                tcorval = ifft(fcorval,nPoints); % See section 6.4 of [1] for relevant algorithms
                absval = (real(tcorval).^2)+(imag(tcorval).^2);
                [accumaxcorr,maxcorridx] = max(absval(1:numCAChips,:).');
                [maxaccumaxcorr,id] = max(accumaxcorr);
                y(iPRN,2) = allfqysteps(maxcorridx(id)) - fc; % Estimate Doppler offset
                corrval(:,:,iPRN) = absval;
                maxcorrvals(iPRN) = max(absval(:));

                % Translate the code phase offset into a fractional number
                % within 1023 chips.
                % nPoints correspond to 1023 chips offset. So, 1
                % corresponds to 1023/nPoints of chips offset.
                y(iPRN,3) = (id-1)*1023/nPoints; % Estimate code-phase offset

                % Logic to check if signal is detected
                y(iPRN,4) = 0;
                if(maxaccumaxcorr > th*maxNoise)
                    y(iPRN,4) = 1;
                end
            end

            % Sort the correlation values and satellites PRN IDS in the
            % descending order of maximum correlation values. This helps in
            % deciding which satellites to choose for tracking, from the
            % pool of detected satellites.
            [maxcorrvals,sortorder] = sort(maxcorrvals,'descend');

            % Sort the values as per maximum correlation so that tracking
            % modules can be assigned to the highly correlated values
            y = y(sortorder,:);
            corrval = corrval(:,:,sortorder);
        end
    end
end
