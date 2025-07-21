classdef HelperGPSCACodeCarrierTracker < matlab.System
    % HelperGPSCACodeCarrierTracker Track the code phase and carrier
    % frequency of the GPS signal
    %
    %   Note: This is a helper function and its API and/or functionality
    %   may change in subsequent releases.
    %
    %   CCTRACK = HelperGPSCACodeCarrierTracker creates a carrier frequency
    %   and code phase tracking object. Frequency locked loop (FLL), phase
    %   locked loop (PLL) and delay locked loop (DLL) are used for tracking
    %   frequency, phase and code phase respectively in the incoming
    %   signal.
    %
    %   Step method syntax:
    %
    %   [Y,FERR,FQY,PERR,PH,DERR,DELAY] = step(CCTRACK,X) tracks the
    %   frequency, phase and delay in the signal X and returns the
    %   integrated version of signal, Y after tracking. Length of X depends
    %   sampling frequency. X must contain data corresponding to
    %   CCTRACK.PLLIntegrationTime time. For example if
    %   CCTRACK.PLLIntegrationTime is 1, then X must be of length
    %   corresponding to 1 millisecond of data. If sampling frequency is
    %   38.192 MHz, then X contains 38192 samples. X can be complex
    %   baseband signal or real IF signal. Y is a complex vector with each
    %   element representing the integrated value of one millisecond of
    %   data in which both carrier and C/A-code are removed. Therefore,
    %   length of Y will be equal to CCTRACK.PLLIntegrationTime. FERR is
    %   the output of frequency discriminator. FQY is the estimated
    %   frequency offset value. PERR is the phase discriminator output. PH
    %   is the detected phase offset value. DERR is the output of delay
    %   discriminator. DELAY is the estimated delay in the units of number
    %   of C/A-code chips.
    %
    %   HelperGPSCACodeCarrierTracker properties:
    %
    %   CenterFrequency        - (Non-tunable) Center frequency of the
    %                            input signal. (Default: 0 Hz)
    %   SampleRate             - (Non-tunable) Sample rate of the input
    %                            signal (Default: 38.192e6 Hz)
    %   PRNID                  - (Non-tunable) PRN ID of the GPS satellite
    %                            that must be tracked. (Range: integers
    %                            from 1 to 32. Default: 1)
    %   InitialCodePhaseOffset - (Tunable) Initial value of code phase
    %                            offset in the incoming signal that is
    %                            detected by initial synchronization
    %                            algorithms. (Range: 0 to 1023. Can include
    %                            fractional values too. Default: 0)
    %   InitialDopplerShift    - (Tunable) Initial value of doppler shift
    %                            in the signal that is detected by initial
    %                            synchronization algorithms. (Default:
    %                            0Hz).
    %   FLLOrder               - (Non-tunable) Order of the frequency
    %                            tracking loop. (Range: 1 or 2. Default:
    %                            1).
    %   FLLNoiseBandwidth      - (Non-tunable) Noise bandwidth of the FLL.
    %                            (Default: 4 Hz)
    %   PLLOrder               - (Non-tunable) Order of the phase tracking
    %                            loop. (Range: 2 or 3. Default: 2)
    %   PLLNoiseBandwidth      - (Tunable) Noise bandwidth of the PLL.
    %                            (Default: 18 Hz)
    %   PLLIntegrationTime     - (Tunable) Integration time for the PLL in
    %                            millisecond units. FLL integration time
    %                            will be half of this value. Typically once
    %                            the bit synchronization is complete,
    %                            integration time is increased for better
    %                            performance. Hence, this is a tunable
    %                            property. (Default: 1 milliseconds)
    %   DisablePLL             - (Tunable) Option to disable PLL. When
    %                            true, PLL is turned off. (Default: false)
    %   DLLOrder               - (Non-tunable) Order of the delay tracking
    %                            loop. (Range: 1 or 2. Default: 1)
    %   DLLNoiseBandwidth      - (Non-tunable) Noise bandwidth of the DLL.
    %                            (Default: 1 Hz)
    %
    %   References:
    %   [1] Kaplan and Hegarty, Understanding GPS/GNSS, chap. 8.

    %   Copyright 2021-2023 The MathWorks, Inc.

    % Public, tunable properties
    properties
        InitialCodePhaseOffset = 0
        InitialDopplerShift = 0
        DisablePLL = false
        PLLIntegrationTime = 1 % In milliseconds
        PLLNoiseBandwidth = 18
    end

    properties(Nontunable)
        % Signal properties
        PRNID = 1
        CenterFrequency = 0
        SampleRate = 38.192e6 % In Hz

        % Properties of carrier tracking loops
        FLLOrder = 1
        PLLOrder = 2
        FLLNoiseBandwidth = 4

        % Properties of code tracking loop
        DLLOrder = 1
        DLLNoiseBandwidth = 1
    end

    properties(Constant,Hidden)
        ChipRate = 1.023e6 % Chip rate of C/A-code
    end

    % Pre-computed constants
    properties(Access = private)
        % FLL properties
        pFLLNaturalFrequency
        pFLLGain1
        pFLLGain2
        pFLLGain3
        pFLLWPrevious1 = 0
        pFLLWPrevious2 = 0
        pFLLNCOOut = 0

        % PLL properties
        pPLLNaturalFrequency
        pPLLGain1
        pPLLGain2
        pPLLGain3
        pPLLWPrevious1 = 0
        pPLLWPrevious2 = 0
        pPLLNCOOut = 0
        pPreviousPhase = 0;

        % DLL properties
        pDLLGain1
        pDLLGain2
        pDLLGain3
        pDLLWPrevious1 = 0
        pDLLNCOOut = 0
        pDLLNaturalFrequency
        pPromptCode

        % General properties
        pNumIntegSamples
        pSamplesPerChip
        pReferenceCode
        pNumSamplesToAppend = 0
        pBuffer
    end

    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            obj.pNumIntegSamples = obj.SampleRate*obj.PLLIntegrationTime*1e-3; % PLLIntegrationTime is in milliseconds. Hence multiply by 1e-3 to get it into sec

            % Calculate loop parameters for DLL
            if obj.DLLOrder == 1 % Table 8.23 in Kaplan 3rd edition [1]
                obj.pDLLNaturalFrequency = obj.DLLNoiseBandwidth/0.25;
                obj.pDLLGain1 = obj.pDLLNaturalFrequency;
            elseif obj.DLLOrder == 2
                obj.pDLLNaturalFrequency = obj.DLLNoiseBandwidth/0.53;
                obj.pDLLGain1 = obj.pDLLNaturalFrequency.^2;
                obj.pDLLGain2 = obj.pDLLNaturalFrequency*1.414;
            else % obj.DLLOrder == 3
                obj.pDLLNaturalFrequency = obj.DLLNoiseBandwidth/0.7845;
                obj.pDLLGain1 = obj.pDLLNaturalFrequency.^3;
                obj.pDLLGain2 = obj.pDLLNaturalFrequency*1.1;
                obj.pDLLGain3 = obj.pDLLNaturalFrequency*2.4;
            end

            % Calculate loop parameters for FLL
            if obj.FLLOrder == 1 % Table 8.23 in Kaplan 3rd edition [1]
                obj.pFLLNaturalFrequency = obj.FLLNoiseBandwidth/0.25;
                obj.pFLLGain1 = obj.pFLLNaturalFrequency;
            elseif obj.FLLOrder == 2
                obj.pFLLNaturalFrequency = obj.FLLNoiseBandwidth/0.53;
                obj.pFLLGain1 = obj.pFLLNaturalFrequency.^2;
                obj.pFLLGain2 = obj.pFLLNaturalFrequency*1.414;
            else % obj.FLLOrder == 3
                obj.pFLLNaturalFrequency = obj.FLLNoiseBandwidth/0.7845;
                obj.pFLLGain1 = obj.pFLLNaturalFrequency.^3;
                obj.pFLLGain2 = obj.pFLLNaturalFrequency*1.1;
                obj.pFLLGain3 = obj.pFLLNaturalFrequency*2.4;
            end

            % Calculate loop parameters for PLL
            if obj.PLLOrder == 1 % Table 8.23 in Kaplan 3rd edition [1]
                obj.pPLLNaturalFrequency = obj.PLLNoiseBandwidth/0.25;
                obj.pPLLGain1 = obj.pPLLNaturalFrequency;
            elseif obj.PLLOrder == 2
                obj.pPLLNaturalFrequency = obj.PLLNoiseBandwidth/0.53;
                obj.pPLLGain1 = obj.pPLLNaturalFrequency.^2;
                obj.pPLLGain2 = obj.pPLLNaturalFrequency*1.414;
            else % obj.PLLOrder == 3
                obj.pPLLNaturalFrequency = obj.PLLNoiseBandwidth/0.7845;
                obj.pPLLGain1 = obj.pPLLNaturalFrequency.^3;
                obj.pPLLGain2 = obj.pPLLNaturalFrequency*1.1;
                obj.pPLLGain3 = obj.pPLLNaturalFrequency*2.4;
            end

            % Initialize the code
            numCACodeBlocks = obj.PLLIntegrationTime; % Each C/A-code block is of 1 milliseconds.
            code = 1-2*double(gnssCACode(obj.PRNID,"GPS"));
            obj.pSamplesPerChip = obj.SampleRate/obj.ChipRate;
            [upSampleFactor,downSampleFactor] = rat(obj.pSamplesPerChip);
            numSamplesPerCodeBlock = obj.SampleRate*1e-3; % As each code block is of 1e-3 seconds
            upwave1 = repelem(code(:),upSampleFactor,1);
            obj.pPromptCode = repmat(upwave1(1:downSampleFactor:end),numCACodeBlocks,1);

            % Calculate number of samples in delay
            numsamprot = round(obj.InitialCodePhaseOffset*obj.pSamplesPerChip); % Number of samples to rotate
            obj.pNumSamplesToAppend = numSamplesPerCodeBlock-mod(numsamprot,numSamplesPerCodeBlock);
            
        end

        function [y,fqyerr,fqynco,pherr,phnco,delayerr,delaynco] = stepImpl(obj,u)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.

            coarsedelay = obj.pNumSamplesToAppend;
            numSamplesPerCodeBlock = obj.SampleRate*1e-3; % As each code block is of 1e-3 seconds
            finedelay = round(obj.pDLLNCOOut*obj.pSamplesPerChip);

            if length(obj.pBuffer) ~= coarsedelay + finedelay
                numextradelay = coarsedelay + finedelay - length(obj.pBuffer);
                if numextradelay > 0
                    obj.pBuffer = [zeros(numextradelay,1);obj.pBuffer];
                else % numextradelay < 0. Equal to zero is not possible because of the first if condition
                    if abs(numextradelay) < length(obj.pBuffer)
                        % Remove samples from pBuffer itself
                        obj.pBuffer(1:abs(numextradelay)) = [];
                    else
                        n = numSamplesPerCodeBlock + numextradelay;
                        obj.pBuffer = [zeros(n,1);obj.pBuffer];
                    end
                end
            end
            
            % Buffer the input
            integtime = obj.PLLIntegrationTime*1e-3; % PLLIntegrationTime is in milliseconds. Hence multiply by 1e-3 to get it into sec
            [u,obj.pBuffer] = buffer([obj.pBuffer;u(:)],round(obj.SampleRate*integtime));

            % Carrier wipe-off
            fc = obj.CenterFrequency + obj.InitialDopplerShift - obj.pFLLNCOOut;
            t = (0:obj.pNumIntegSamples)/obj.SampleRate;
            phases = 2*pi*fc*t.' + obj.pPreviousPhase - obj.pPLLNCOOut;
            iqsig = u.*exp(-1j*phases(1:end-1));
            obj.pPreviousPhase = phases(end) + obj.pPLLNCOOut;

            % Code wipe-off
            % Update the prompt code appropriately
            
            numSamplesPerHalfChip = round(obj.pSamplesPerChip/2);
            iq_e = iqsig.*circshift(obj.pPromptCode,-1*numSamplesPerHalfChip);
            iq_p = iqsig.*obj.pPromptCode;
            iq_l = iqsig.*circshift(obj.pPromptCode,numSamplesPerHalfChip);
            integeval = sum(iq_e);
            integlval = sum(iq_l);

            millisecdata = reshape(iq_p,[],obj.PLLIntegrationTime); % Each column contains one millisecond of data
            y = sum(millisecdata); % Each element contains integrated value of one millisecond of data
            integpval = sum(y);
            if rem(length(iq_p),2)~=0 % Odd number of samples
                fllin = sum(reshape([iq_p;0],[],2)); % Append a zero
            else
                fllin = sum(reshape(iq_p,[],2));
            end

            % DLL discriminator
            E = abs(integeval);
            L = abs(integlval);
            delayerr = (E-L)/(2*(E+L)); % Non-coherent early minus late normalized detector

            % DLL loop filter
            if obj.DLLOrder == 2
                % 1st integrator
                wcurrent = delayerr*obj.pDLLGain1*integtime + obj.pDLLWPrevious1;
                loopfilterout = (wcurrent + obj.pDLLWPrevious1)/2 + delayerr*obj.pDLLGain2;
                obj.pDLLWPrevious1 = wcurrent; % Acceleration accumulator
            elseif obj.DLLOrder == 1
                loopfilterout = delayerr*obj.pDLLGain1;
            end

            % DLL NCO
            delaynco = obj.pDLLNCOOut + integtime*loopfilterout;
            obj.pDLLNCOOut = delaynco;
            
            % FLL discriminator
            phasor = conj(fllin(1))*fllin(2);
            % phasor = conj(obj.pPreviousIntegPVal)*integpval;
            fqyerr = -1*angle(phasor)/(pi*integtime); % Multiplication by 2 is removed because integtime of FLL is half of that of PLL

            % FLL loop filter
            if obj.FLLOrder == 2
                % 1st integrator
                wcurrent = fqyerr*obj.pFLLGain1*integtime + obj.pFLLWPrevious1;
                loopfilterout = (wcurrent + obj.pFLLWPrevious1)/2 + fqyerr*obj.pFLLGain2;
                obj.pFLLWPrevious1 = wcurrent; % Acceleration accumulator
            elseif obj.FLLOrder == 1
                loopfilterout = fqyerr*obj.pFLLGain1;
            end

            % FLL NCO
            fqynco = obj.pFLLNCOOut + integtime*loopfilterout;
            obj.pFLLNCOOut = fqynco;

            % PLL discriminator
            if obj.DisablePLL
                pherr = 0;
            else
                pherr = atan(real(integpval)/imag(integpval));
            end

            % PLL loop filter
            if obj.PLLOrder == 3
                % 1st integrator
                wcurrent = pherr*obj.pPLLGain1*integtime + obj.pPLLWPrevious1;
                integ1out = (wcurrent + obj.pPLLWPrevious1)/2 + pherr*obj.pPLLGain2;
                obj.pPLLWPrevious1 = wcurrent; % Acceleration accumulator

                % 2nd integrator
                wcurrent = integ1out*integtime + obj.pPLLWPrevious2;
                loopfilterout = (wcurrent + obj.pPLLWPrevious2)/2 + pherr*obj.pPLLGain3;
                obj.pPLLWPrevious2 = wcurrent; % Velocity accumulator
            elseif obj.PLLOrder == 2
                wcurrent = pherr*obj.pPLLGain1*integtime + obj.pPLLWPrevious1;
                loopfilterout = (wcurrent + obj.pPLLWPrevious1)/2 + pherr*obj.pPLLGain2;
                obj.pPLLWPrevious1 = wcurrent; % Velocity accumulator
            end

            % PLL NCO
            phnco = obj.pPLLNCOOut + integtime*loopfilterout;
            obj.pPLLNCOOut = phnco;
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            obj.pBuffer = zeros(obj.pNumSamplesToAppend,1);
            obj.pFLLWPrevious1 = 0;
            obj.pFLLWPrevious2 = 0;
            obj.pFLLNCOOut = 0;
            obj.pPLLWPrevious1 = 0;
            obj.pPLLWPrevious2 = 0;
            obj.pPLLNCOOut = 0;
            obj.pDLLWPrevious1 = 0;
            obj.pDLLNCOOut = 0;
        end

        function processTunedPropertiesImpl(obj)
            % Perform actions when tunable properties change
            % between calls to the System object
            obj.pNumIntegSamples = obj.SampleRate*obj.PLLIntegrationTime;

            [upSampleFactor,downSampleFactor] = rat(obj.SampleRate/obj.ChipRate);
            numSamplesPerCodeBlock = obj.SampleRate*1e-3;
            numsamprot = round(obj.InitialCodePhaseOffset*upSampleFactor/downSampleFactor); % Number of samples to rotate
            obj.pNumSamplesToAppend = numSamplesPerCodeBlock-mod(numsamprot,numSamplesPerCodeBlock);

            % Calculate loop parameters for PLL
            if obj.PLLOrder == 1 % Table 8.23 in Kaplan 3rd edition [1]
                obj.pPLLNaturalFrequency = obj.PLLNoiseBandwidth/0.25;
                obj.pPLLGain1 = obj.pPLLNaturalFrequency;
            elseif obj.PLLOrder == 2
                obj.pPLLNaturalFrequency = obj.PLLNoiseBandwidth/0.53;
                obj.pPLLGain1 = obj.pPLLNaturalFrequency.^2;
                obj.pPLLGain2 = obj.pPLLNaturalFrequency*1.414;
            else % obj.PLLOrder == 3
                obj.pPLLNaturalFrequency = obj.PLLNoiseBandwidth/0.7845;
                obj.pPLLGain1 = obj.pPLLNaturalFrequency.^3;
                obj.pPLLGain2 = obj.pPLLNaturalFrequency*1.1;
                obj.pPLLGain3 = obj.pPLLNaturalFrequency*2.4;
            end
        end
    end

    methods(Access = protected, Static)
        function group = getPropertyGroupsImpl

            gengroup = matlab.system.display.SectionGroup('PropertyList', ...
                {'CenterFrequency', 'SampleRate'});
            initprops = {'PRNID','InitialCodePhaseOffset','InitialDopplerShift'};
            % Define property section(s) for System block dialog
            initGroupTitle = 'Outputs from initial synchronization';
            initgroup = matlab.system.display.SectionGroup('Title', ...
                initGroupTitle, 'PropertyList', initprops);

            fllprops = {'FLLOrder','FLLNoiseBandwidth'};
            % Define property section(s) for System block dialog
            fllGroupTitle = 'Frequency locked loop';
            fllgroup = matlab.system.display.SectionGroup('Title', ...
                fllGroupTitle, 'PropertyList', fllprops);

            pllprops = {'PLLOrder','PLLNoiseBandwidth','PLLIntegrationTime','DisablePLL'};
            pllGroupTitle = 'Phase locked loop';
            pllgroup = matlab.system.display.SectionGroup('Title', ...
                pllGroupTitle, 'PropertyList', pllprops);

            dllprops = {'DLLOrder','DLLNoiseBandwidth'};
            dllGroupTitle = 'Delay locked loop';
            dllgroup = matlab.system.display.SectionGroup('Title', ...
                dllGroupTitle, 'PropertyList', dllprops);
            group = [gengroup,initgroup,fllgroup,pllgroup,dllgroup];
        end
    end
end
