classdef HelperGPSLNAVFrameSynchronizer < matlab.System
    % HelperGPSLNAVFrameSynchronizer GPS LNAV data frame synchronize
    %
    %   Note: This is a helper and its API and/or functionality may change
    %   in subsequent releases.
    %
    %   GPSFS = HelperGPSLNAVFrameSynchronizer creates a Global Positioning
    %   System (GPS) frame synchronization object, GPSFS that calculates the frame
    %   boundary in the incoming GPS legacy navigation (LNAV) data bits.
    %
    %   Step method syntax:
    %
    %   [SYNCIDX,Y,SUBFRAMEIDS] = step(GPSFS,X) Calculates the frame
    %   boundary using the 8 bits preamble that is constant for all
    %   subframes. Once preamble is detected, further checks such as
    %   decoding first 60 bits of subframe, calculating the time of week
    %   so as to confirm the successful frame synchronization. Because the
    %   costas loop used for tracking a GPS signal has a 180 degrees phase
    %   ambiguity, frame synchronization works on both 8 bits preamble and
    %   the inverted version of preamble and corrects the data of phase
    %   ambiguity appropriately. SYNCIDX is the index of subframe beginning
    %   in the input data bits, X. Y is the synchronized subframes with
    %   number of rows equal to 300 and number of columns equal to the
    %   number of subframes that are successfully detected in the input
    %   data X. SUBFRAMEIDS is a row vector indicating the subframe IDs of
    %   each row in Y.
    %
    %   There are no properties for this System object.

    %   Copyright 2021-2022 The MathWorks, Inc.

    properties(Constant, Hidden)
        pTLMPreamble = [1; 0; 0; 0; 1; 0; 1; 1]
        pInvTLMPreamble = [0; 1; 1; 1; 0; 1; 0; 0]
    end

    % Pre-computed constants
    properties(Access = private)
        pPreviousBits
        pIsFrameSyncSuccess
        pHOWTOW = 0
        pSyncIdx
        pCalculateSyncIdx = true
        pAlreadyInverted = false
    end

    methods(Access = protected)

        function [syncidx,y,subframeIDs] = stepImpl(obj,u)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.

            bitsPerSubframe = 300;
            [bits,obj.pPreviousBits] = buffer([obj.pPreviousBits;u],bitsPerSubframe); % Input size must be a multiple of 300
            
            y = zeros(size(bits));
            subframeIDs = zeros(1,size(y,2));

            isubframe = 1;
            cntr = 1;
            while(1)
                % Perform frame synchronization only if frame sync is not
                % successful
                syncbits = bits(:,isubframe);
                if obj.pIsFrameSyncSuccess
                    % Just calculate next HOWTOW value and check if frame sync is behaving
                    % properly
                    word1 = syncbits(1:30);
                    [~,parityBits,parityCheck1] = HelperGPSLNAVWordDecode(word1,1,[0;0]);

                    word2 = syncbits(31:60);
                    [word2dec,~,parityCheck2] = HelperGPSLNAVWordDecode(word2,2,parityBits(end-1:end,1));
                    HOWTOW = bit2int(word2dec(1:17),17);
                    sID = bit2int(word2dec(20:22),3); % Subframe ID
                    if ((~parityCheck2) || (~parityCheck1))
                        obj.pIsFrameSyncSuccess = false;
                    else % Frame sync successful
                        y(:,cntr) = syncbits;
                        subframeIDs(cntr) = sID;
                        cntr = cntr + 1;
                        isubframe = isubframe + 1;
                        % obj.pCalculateSyncIdx = true;
                        obj.pHOWTOW = HOWTOW;
                    end
                    
                end

                if ~obj.pIsFrameSyncSuccess
                    tempbits = [reshape(bits((isubframe-1)*bitsPerSubframe+1:end),[],1);obj.pPreviousBits];
                    [syncidxtemp,sID,IsInvertedMatched] = framesynccore(obj,tempbits(:));

                    if IsInvertedMatched
                        tempbits = double(xor(tempbits,1));
                        obj.pAlreadyInverted = true;
                    end

                    if obj.pIsFrameSyncSuccess
                        tempbits = tempbits(syncidxtemp:end).';
                        [bits,obj.pPreviousBits] = buffer(tempbits(:),bitsPerSubframe);
                        if ~isempty(bits)
                            y(:,cntr) = bits(:,1);
                            subframeIDs(cntr) = sID;
                            cntr = cntr + 1;
                            isubframe = 2; % First one is already in sync and accounted for
                        end
                    else
                        obj.pPreviousBits = tempbits(end-59:end); % Store last two words in buffer as they are not searched
                        break;
                    end
                end

                if isubframe > size(bits,2)
                    break;
                end
            end
            if size(y,2)>=cntr
                y = y(:,1:cntr-1);
                subframeIDs = subframeIDs(:,1:cntr-1);
            end
            syncidx = obj.pSyncIdx;
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            obj.pPreviousBits = zeros(0,1);
            obj.pIsFrameSyncSuccess = false;
            obj.pHOWTOW = 0;
            obj.pSyncIdx = 0;
            obj.pCalculateSyncIdx = true;
            obj.pAlreadyInverted = false;
        end

        function [syncidx, sid, IsInvertedMatched] = framesynccore(obj,u)
            % Match the incoming bits with either 0x8B or 0x74. Length of u
            % must be at least 300. Or this function can behave
            % unexpectedly.
            sid = 0; % Pre-initialize
            numMatchBits = 60;
            syncidx = 0;
            for ibit = 1:length(u)-numMatchBits+1
                reg = u(ibit+(0:7));
                IsRegularMatched = isequal(reg,obj.pTLMPreamble);
                IsInvertedMatched = isequal(reg,obj.pInvTLMPreamble);
                if IsRegularMatched || IsInvertedMatched
                    % Perform further checks to verify if data decoding is
                    % truly successful
                    % Based on which sequence matched, either invert all
                    % bits or keep the bits intact
                    if IsInvertedMatched
                        bits = double(xor(u,1));
                    else
                        bits = u;
                    end


                    % Perform parity check
                    word1 = bits(ibit+(0:29));
                    [~,parityBits,parityCheck1] = HelperGPSLNAVWordDecode(word1,1,[0;0]);

                    word2 =bits(ibit+30+(0:29));
                    [word2dec,~,parityCheck2] = HelperGPSLNAVWordDecode(word2(:),2,parityBits(end-1:end,1));

                    % If either of parity check fails, continue the frame
                    % search process
                    if parityCheck1+parityCheck2 == 2 % Condition that parity check passed
                        
                        HOWTOW = bit2int(word2dec(1:17),17);
                        subframeID = bit2int(word2dec(20:22),3);

                        if mod(HOWTOW-1,5)+1 == subframeID % Check if HOWTOW and subframe ID agree
                            % With confidence, frame sync is successful
                            obj.pIsFrameSyncSuccess = true;
                            sid = subframeID; % Subframe ID
                            syncidx = ibit;
                            if obj.pCalculateSyncIdx == true && ~obj.pAlreadyInverted
                                obj.pSyncIdx = syncidx;
                                obj.pCalculateSyncIdx = false;
                            end
                            
                            obj.pHOWTOW = HOWTOW;
                            break;
                        end
                    end
                end
            end
        end
    end
end
