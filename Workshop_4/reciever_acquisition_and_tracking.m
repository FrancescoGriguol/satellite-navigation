%% Inizializzazione
% parametri da cambiare: frequenza di campionamento, bande passanti dei 
% tracking loops, aumentare il tempo della simulazione da tre a più secondi 
% initialsync.DetectionThresholdFactor = 1.9;

clear
clc
close all

% Inizializzazione e lettura dei dati di navigazione RINEX
addpath("Acquisizione_SwiftNav_20240410\")
rinexdata = rinexread('swift-gnss-20240410-130631.sbp.nav');
% Si rimuovono le costellazioni che non sono di nostro interesse
rinexdata = rmfield(rinexdata,'Galileo');
rinexdata = rmfield(rinexdata,'BeiDou');
rinexdata = rmfield(rinexdata,'SBAS');

%% Configurazione della simulazione
% Durata della simulazione
simulatedDataDuration = 3; % In seconds
% Frequenza di campionamento
samplingRate = 10.23e6;    % In Hz

% Numero di passi della simulazione
stepTime = 1e-3;   % In seconds
numSteps = (simulatedDataDuration/stepTime) + 1;

%% Parametri di posizione e di tracking
% Posizione di Padova in LLA dove è stata fatta la misurazione
rxlat = 45.410885220515276; % Latitude in degrees north
rxlon = 11.891726751826376; % Longitude in degrees east; negative longitude is degrees west
rxel  = 70.041;             % Elevation in meters

% Numero di blocchi del codice C/A per bit
numCACodeBlocksPerBit = 20; % del navigation data, per ogni bit ci sono 20 repliche del C/A code
% Tempo di attesa
rxWaitTime = 100;       % Milliseconds
% Flag per sincro iniziale
performInitSync = true; % Initially this must be set to true

% Initialize maximum number of tracking channels. 
% (Minimum of 4 tracking channels are needed for proper functioning of the
% GPS receiver)
maxNumTrackingChannels = 8; 


% Noise bandwidth of each of the tracking loops
PLLNoiseBW = 90; % In Hz
FLLNoiseBW = 4;  % In Hz
DLLNoiseBW = 3;  % In Hz
% (Scelti in base al compromesso che aumentando la banda passante il rumore
% viene seguito molto velocemente, ma contemporanemente lo si riduce poco.
% Al contrario diminuendo la banda passande si riesce a filtrare meglio il 
% rumore, ma viene seguito più lentamente)

% Bit synchronization parameters
isBitSyncComplete = zeros(maxNumTrackingChannels,1);
numBitsForBitSync = 100; % bit necessari per la sincronizzazione
numWaitingStepsForBitSync = numCACodeBlocksPerBit*numBitsForBitSync;
rxcntr = 1;

minTimeForPosEst = 48.5;                         % In seconds
minStepsForPosEst = minTimeForPosEst/stepTime;
subframeDuration = 6;                            % In seconds
numStepsPerSubframe = subframeDuration/stepTime;

%% Parametri fisici per la stima del rumore termico
% parametri fisici per determinare il thermal noise a livello del
% ricevitore, rumore che viene aggiunto alla waveform del segnale una volta
% raggiunto all'antenna ricevente
c = physconst("LightSpeed"); % Speed of light in m/sec
fe = 1575.42e6;              % GPS L1 frequency in Hz
Dt = 12;                     % Directivity of the transmit antenna in dBi
DtLin = db2pow(Dt);
Dr = 4;                      % Directivity of the receive antenna in dBi
DrLin = db2pow(Dr);
Pt = 44.8;                   % Typical transmission power of a GPS satellite in watts
k = physconst("boltzmann");  % Boltzmann constant in Joules/Kelvin
T = 300;                     % Room temperature in Kelvin
rxBW = 24e6;                 % Bandwidth in Hz
Nr = k*T*rxBW;               % Thermal noise power in watts
rng default;                 % Initializing to default random number generation

%% Satellite scenario e dati di navigazione
% Initialize satellite scenario
sc = satelliteScenario;

% Set up the satellites based on the RINEX data and ground station
sat = satellite(sc,rinexdata,"OrbitPropagator","gps");
rx = groundStation(sc,rxlat,rxlon); % Set up the receiver
ac = access(sat,rx);                % Calculate access between satellites and the receiver


% Get the list of satellites that are considered in satellite scenario from
% the RINEX data
indices = ones(length(sat),1);
for isat = 1:length(sat)
    ele = orbitalElements(sat(isat));

    % Check for match of time of applicability and pseudo-random noise
    % (PRN) IDs so that data from RINEX file is considered for waveform
    % generation
    indices(isat) = find(rinexdata.GPS(:,:).Toe == ele.GPSTimeOfApplicability & ...
        rinexdata.GPS(:,:).SatelliteID == ele.PRN);
end

% play satellite scenario ( si visualizzano i satelliti e le loro orbite
% assieme alla ground station e ai vettori che collegano la ground station
% con i satelliti )

%play(sc)

%% Download e configurazione dell'almanacco GPS
d = datetime(2024,04,10);

baseURL = "https://www.navcen.uscg.gov/sites/default/files/gps/almanac/";
almanacType = "/sem/";
almanacExtension = ".al3";
url = baseURL + d.Year + almanacType + day(d,"dayofyear") + ...
      almanacExtension;

filename = "semalmanac" + "_" + d.Year + "-" + ...
           d.Month + "-" + d.Day + ".al3";

websave(filename,url);

% Generate navigation configuration object
% (Funzione che crea un oggetto di configurazione di navigazione GPS usando
% il file dell'almanacco scaricato e i dati RINEX filtrati (solo GPS)
navcfg = HelperGPSRINEX2Config(filename,rinexdata.GPS(indices,:));

[mintow,locmintow] = min([navcfg(:).HOWTOW]);

% The time of week (TOW) value is such that it must be 1 plus a multiple of
% 5 so that data of subframe 1 is always generated first.
mintow = ceil((mintow-1)/5)*5 + 1;

% HOWTOW is a counter that contains integer values. Incrementing by a value
% of 1 represents 6 seconds of data. The counter increases by a value of 1
% for each new subframe.
% Si regola il time of the week per iniziare sempre con il subframe 1
[navcfg(:).HOWTOW] = deal(mintow);
% Set the starting of frames based on mintow
firstsubframeID = mod(mintow-1,125) + 1;
frameID = ceil(firstsubframeID/5);
allFrameIDs = [frameID:25,1:(frameID-1)];
[navcfg(:).FrameIndices] = deal(allFrameIDs);

% Generate GPS navigation data
% Numero totale di bit nei dati di navigazione GPS
numNavBits = 37500;   
% matrice per i dati di navigazione, una colonna per ogni satellite
navdata = zeros(numNavBits,length(navcfg)); 
for isat = 1:length(navcfg)
    % funzione che genera i dati di navigazione per ciascun satellite
    % basandosi sulla configurazione navcfg
    navdata(:,isat) = HelperGPSNAVDataEncode(navcfg(isat));
end

% impostazione sul tempo di inizio e di fine della simulazione
sc.StartTime = HelperGPSConvertTime(navcfg(locmintow).WeekNumber, ...
   mintow*subframeDuration);

sc.StopTime = sc.StartTime + seconds(simulatedDataDuration);
sc.SampleTime = stepTime;

%play(sc)

acstats = accessStatus(ac);
% This example runs for a maximum of 2 minutes of data and in that
% duration, access does not change. Hence, consider only the access status
% of the first sample time.
satindices = find(acstats(:,1));
% numero di satelliti con accesso
numsat = length(satindices);

%% Doppler shift

% Calculate Doppler shift over time for all the visible satellites
[fShift,timeOut,dopplerInfo] = dopplershift(sat(satindices),rx,Frequency=fe);

figure(1)
plot(timeOut,fShift(:,:)/1e3,'LineWidth',1)   % Doppler shift in kilohertz (kHz)
xlim([timeOut(1) timeOut(end)])
title("Doppler shift over time",'FontSize',20)
xlabel("Simulation Time",'FontSize',16)
ylabel("Doppler Shift [kHz]",'FontSize',16)
grid on
satName = [];
for isat = 1:length(navcfg)
    satName = [satName; ['GPS ID ' sprintf('%2i',indices(isat))]];
end
legend(satName)



%% Distance between satellites and receiver, pathloss and delay
% Get the states of satellites in the sky over time. Also, shuffle the
% dimensions such that each row corresponds to one satellite data
% ( states restituisce le posizioni dei satelliti in coordinate ECEF mentre
% permute sposata le dimensioni dell'array mettendo il tempo alla prima
% posizione)
initsatpos = permute(states(sat,"CoordinateFrame","ecef"),[3 1 2]);

% Compute the distance between satellites and receiver over time
satdist = zeros(numsat,numSteps);
for istep = 1:numSteps
    satdist(:,istep) = pseudoranges([rxlat,rxlon,rxel], ...
        initsatpos(satindices,:,istep),"RangeAccuracy",0);
end

% Plot delle distanze del satellite in funzione del tempo
figure(2)
plot(timeOut,satdist/1e3,'LineWidth',1) 
xlim([timeOut(1) timeOut(end)])
title("Distance between satellites and receiver over time",'FontSize',20)
xlabel("Simulation Time",'FontSize',16)
ylabel("Distance [km]",'FontSize',16)
grid on
legend(satName)

% Ritardo del segnale
delays = satdist/c;

%% C/A Code
% Power at receiver from free space pathloss equation
SqrtPr = sqrt(Pt*DtLin*DrLin)*(1./(4*pi*(fe+fShift).*delays));
timeinweek = mintow*6;
PRNIDs = [navcfg(:).PRNID];
disp("Available satellites - " + num2str(PRNIDs(satindices)))

% Generazione dei C/A codes
cacodes = gnssCACode(PRNIDs,"GPS");
caCodeBlockDuration = 1e-3;                     % Constant value
numcacodeblocks = stepTime/caCodeBlockDuration; % Each C/A-code block is of 1 millisecond duration

% Plot del C/A code per un dato satellite
% index = randsample(satindices,1); 
selected_cacodes = cacodes(:,1);%index
figure(3)
plot(selected_cacodes, 'LineWidth',0.2);
title(['C/A Code for Satellite PRN ', num2str(PRNIDs(1))], 'FontSize', 20);%index
xlabel('Chip Index', 'FontSize', 16);
xlabel('Chip Index', 'FontSize', 16);
ylabel('Code Value', 'FontSize', 16);
yticks([0 1]);
ylim([-1 2])
xlim([0 1023])
xticks([0 200 400 600 800 1023])
grid on;

%% Rate match

% To rate match the C/A-code with P-code, repeat each element 10 times
% (si adatta la frequenza del nostro codice C/A con quella del P code)
cacodesig1 = double(repelem(1-2*repmat(cacodes,numcacodeblocks,1),10,1));

% Rate match the C/A-code with the sampling rate
% (si fa la stessa cosa ma per match con la frequenza di campionamento
% specificata)
[upfac,downfac] = rat(samplingRate/10.23e6);
idx = 1:downfac:upfac*size(cacodesig1,1);
cacodesig = cacodesig1(ceil(idx/upfac),:);

% P CODE IS NOT USED IN THIS EXAMPLE! A FAKE P CODE IS GENERATE TO COMPLETE
% THE WAVEFORM GENERATION
numSamples = stepTime*samplingRate;
pcode = (1-2*repmat([1;0],numSamples/2,1))/sqrt(2);
% Serve per compensare l'offset di fase e frequenza
pfo = comm.PhaseFrequencyOffset(FrequencyOffsetSource = "Input port", ...
    SampleRate = samplingRate);

% Initialize static delay object
dynamicDelayRange = 20000;
staticdelay = round(delays(:,1)*samplingRate - dynamicDelayRange);
if nnz(staticdelay<0)~=0
    staticdelay = zeros(numsat,1);
end
staticDelayObj = dsp.Delay("Length",staticdelay);

% Initialize Variable Fractional Delay object for modeling dynamic delay
vfd = dsp.VariableFractionalDelay("InterpolationMethod","Farrow", ...
    "MaximumDelay",65535);
% Diagramma di costellazione del segnale tracciato
rxconstellation = comm.ConstellationDiagram(1,ShowReferenceConstellation=false, ...
    Title="Constellation diagram of signal at the output of tracking");

initialsync = gnssSignalAcquirer;
initialsync.SampleRate = samplingRate;
% imposta la frequenza intermedia a zero
initialsync.IntermediateFrequency = 0;        % Baseband signal
% initialsync.DetectionThresholdFactor = 1.9;
% Un valore più alto di questo fattore implica una soglia più rigorosa, 
% rendendo il rilevamento del segnale più selettivo. Al contrario, un valore
% più basso rende il rilevamento più sensibile, permettendo di captare 
% segnali più deboli, ma aumentando anche il rischio di falsi positivi.

% Properties that store outputs of initial synchronization
[doppleroffsets,codephoffsets] = deal(zeros(1,maxNumTrackingChannels));

% Properties required for storing outputs from tracking module
[accuPh,accuFqy,accuFqyErr,accuPhErr,accuIntegWave,accuDelay,accuDelayErr] = ...
                deal(zeros(numSteps,maxNumTrackingChannels));

% Properties to store outputs of bit synchronization and frame
% synchronization
[maxTransitionLocation,sampleCounter] = deal(zeros(maxNumTrackingChannels,1));
syncidx = zeros(maxNumTrackingChannels,1);

% Property to store output of data decoder
deccfg = cell(maxNumTrackingChannels,1);

% Initialize maximum number of steps for which the simulation chain runs
maxSimSteps = 50/stepTime;

tic % Start of simulation
for istep = 1:numSteps
    %% Generate waveform
    % calcolo indice bit di navigazione corrente
    bitidx = floor((istep-1)/numCACodeBlocksPerBit)+1;
    % crea p code per tutti i satelliti
    allpcode = repmat(pcode,1,numsat);
    % Get navigation bit of each satellite at the corresponding step time
    % estrae e converte i bit di navigazione
    navbitsig = 1-2*navdata(bitidx,satindices);
    % combina p code c/a code e i bit di navigazione per generare il
    % segnale I/Q 
    iqsig = (allpcode + 1j*cacodesig(:,satindices)).*navbitsig; % Implicit expansion !!!!!!!

    %% Propagation channel
    % Model the channel which models Doppler, delay, path loss, and noise on
    % the signal

    % Introduce Doppler to the signal as a frequency offset
    dopsig = pfo(iqsig,fShift(:,istep).');
    % si introduce un offset di frequenza doppler al segnale I/Q 

    % Introduce variable fractional delay
    staticDelayedSignal = staticDelayObj(dopsig);
    leftoutDelay = delays(:,istep)*samplingRate - staticdelay; % Value must always be positive
    delayedSig = vfd(staticDelayedSignal,leftoutDelay.');
    % si introduce un ritardo frazionario variabile al segnale dopsig

    % Scale the delayed signal as per received power calculated
    rmsPow = rms(delayedSig);
    rmsPow(rmsPow==0) = 1;                               % To avoid division by zero
    scaledsig = (SqrtPr(:,istep).').*delayedSig./rmsPow;
    % si scala il segnale ritardato in base alla potenza ricevuta

    % Add up all the signals at the receiver
    % è il segnale risultante dalla somma di tutti i segnali scalati per
    % ciascun satellite
    resultsig = sum(scaledsig,2);

    % Generate noise
    % segnale di rumore bianco gaussiano
    noisesig = (wgn(numSamples,1,10*log10(Nr)) + 1j*wgn(numSamples,1,10*log10(Nr)))./sqrt(2);

    % Add constant thermal noise to the composite signal
    rxwaveform = resultsig + noisesig;

    % Scale the received signal for having unit power
    waveform = rxwaveform/rms(rxwaveform);
    % segnale ricevuto scalato per avere potenza unitaria

    %% Receiver
    % Because there are large delays (70 to 80 milliseconds) modeled on the
    % signal, start the receiver after some time only so as to process
    % valid signal instead of pure noise.
    if istep > rxWaitTime
        if performInitSync == 1
            performInitSync = 0;
            [y,corrval] = initialsync(waveform,1:32); % Initial synchronization
            display(y);


            PRNIDsToSearch = y(y(:,4).IsDetected==1,1).PRNID.';
            PRNIDsToSearch = sort(PRNIDsToSearch);
            % Offset di frequenza doppler per i satelliti rilevati
            doppleroffsets = y(y(:,4).IsDetected==1,2).FrequencyOffset;
            % offset di fase del codice per i satelliti rilevati
            codephoffsets = y(y(:,4).IsDetected==1,3).CodePhaseOffset;

            numdetectsat = length(PRNIDsToSearch);
            if numdetectsat > maxNumTrackingChannels
                % Set this value to limit the number of tracking channels
                numdetectsat = maxNumTrackingChannels;
            end

            % Perform the required initialization of tracking modules for each
            % channel and buffers to store the data.
            disp("The detected satellite PRN IDs: " + num2str(PRNIDsToSearch))

            % Plot the correlation plot for the first satellite
            figure(4);
            mesh(-10e3:500:10e3, 0:size(corrval,1)-1, corrval(:,:,PRNIDsToSearch(1)));%index
            xlabel("Doppler Offset",'FontSize',16)
            ylabel("Code Phase Offset","FontSize",16)
            zlabel("Correlation",'FontSize',16)
            msg = ["Correlation Plot for PRN ID: " num2str(PRNIDsToSearch(1))];%index
            title(msg,'FontSize',20)

             % Initialize frame synchronization objects for each detected satellite
            framesyncbuffer = cell(1,numdetectsat);

            % Create a cell array, where each element corresponds to a carrier
            % tracking object.
            carrierCodeTrack = cell(numdetectsat,1);
            framesync = cell(numdetectsat,1);

            % Update properties for each tracking loop
            for isat = 1:numdetectsat
                carrierCodeTrack{isat} = HelperGPSCACodeCarrierTracker;
                carrierCodeTrack{isat}.SampleRate = samplingRate;
                carrierCodeTrack{isat}.CenterFrequency = 0;
                carrierCodeTrack{isat}.PLLNoiseBandwidth = PLLNoiseBW;
                carrierCodeTrack{isat}.FLLNoiseBandwidth = FLLNoiseBW;
                carrierCodeTrack{isat}.DLLNoiseBandwidth = DLLNoiseBW;
                carrierCodeTrack{isat}.PLLIntegrationTime = 1; % In milliseconds
                carrierCodeTrack{isat}.PRNID = PRNIDsToSearch(isat);
                carrierCodeTrack{isat}.InitialDopplerShift = doppleroffsets(isat);
                carrierCodeTrack{isat}.InitialCodePhaseOffset = codephoffsets(isat);

                % Initialize frame synchronization object
                framesync{isat} = HelperGPSLNAVFrameSynchronizer;
            end
        end
       
        % Because it would be sufficient to get receiver position after
        % running the simulation for 50 seconds of data, stop the loop
        % after executing for 50 seconds of data. In the default run of the
        % example, only 3 seconds of data is processed and this line is not
        % needed. This check is when running the example for large data
        % duration only (for at least 50 seconds of data).
        if istep > maxSimSteps
            break;
        end

        for isat = 1:numdetectsat % Perform tracking for each satellite

            [integwave,fqyerr,fqyoffset,pherr,phoffset,derr,dnco] = ...
                carrierCodeTrack{isat}(waveform);

            % Accumulate the values to see the results at the end
            accuFqyErr(rxcntr,isat) = fqyerr;
            accuFqy(rxcntr,isat) = fqyoffset;
            accuPhErr(rxcntr,isat) = pherr;
            accuPh(rxcntr,isat) = phoffset;
            accuIntegWave(rxcntr,isat) = sum(integwave);
            accuDelayErr(rxcntr,isat) = derr;
            accuDelay(rxcntr,isat) = dnco;
        end
        % Perform bit synchronization, frame synchronization, and data
        % decoding if numWaitingStepsForBitSync of receiver steps are
        % complete.
        if rxcntr > numWaitingStepsForBitSync
            % For each detected satellite, perform bit synchronization,
            % frame synchronization, and data decoding
            for isat = 1:numdetectsat
                if ~isBitSyncComplete(isat)
                    maxTransitionLocation(isat) = ...
                        gnssBitSynchronize( ...
                        imag(accuIntegWave(1:numWaitingStepsForBitSync,isat)), ...
                        numCACodeBlocksPerBit);
                    isBitSyncComplete(isat) = 1;
                    sampleCounter(isat) = rxcntr - maxTransitionLocation(isat) + 1;
                    framesyncbuffer{isat} = accuIntegWave( ...
                        maxTransitionLocation(isat):end,isat);
                else % Perform frame synchronization and data decoding
                    sampleCounter(isat) = sampleCounter(isat) + 1;
                    framesyncbuffer{isat}(sampleCounter(isat)) = accuIntegWave(rxcntr,isat);
                    if mod(sampleCounter(isat),numStepsPerSubframe) == 0
                        samples = framesyncbuffer{isat}(sampleCounter(isat) - ...
                            numStepsPerSubframe+1:sampleCounter(isat));
                        sym = mean(reshape(samples,numCACodeBlocksPerBit,[]));
                        bits = imag(sym)<0;
                        [syncidx(isat),rxsubframes,subframeIDs] = framesync{isat}(bits(:));
                        if ~isempty(rxsubframes) % Then perform data decoding
                            deccfg{isat}.PRNID = PRNIDsToSearch(isat);
                            deccfg{isat} = HelperGPSLNAVDataDecode(rxsubframes,deccfg{isat});
                        end
                    end
                end
            end
        end

        if mod(rxcntr,1000) == 0
            disp("Processed " + (rxcntr/1000) + " sec of data at the receiver.")
            rxconstellation(accuIntegWave(rxcntr-999:rxcntr,1)/ ...
                rms(accuIntegWave(rxcntr-999:rxcntr,1)))
        end

        % Update rxcntr
        rxcntr = rxcntr + 1;


    end
   
 
end

 ifscope = spectrumAnalyzer(SampleRate = samplingRate, ...
        PlotAsTwoSidedSpectrum = true, ...
        SpectrumType = "Power", ...
        SpectrumUnits = "dBW", ...
        Title = "IF Spectrum Comparison of GPS Signal with Thermal Noise", ...
        ShowLegend = true, ...
        ChannelNames = ["GPS IF waveform spectrum" "Noise spectrum"], ...
        YLimits = [-190 -155]);
    ifscope([resultsig, noisesig]);

rxscope = spectrumAnalyzer(SampleRate = samplingRate, ...
        PlotAsTwoSidedSpectrum = true, ...
        SpectrumType = "Power", ...
        SpectrumUnits = "dBW", ...
        Title = "Received signal IF spectrum after scaling");
    rxscope(waveform);


caChipRate = 1.023e6; % In Hz
codeOffsetTime = codephoffsets(1:numdetectsat)/caChipRate;
rxcntr(rxcntr<=1) = 2; % So that rxcntr-1 is at least 1
trackingOffsetTime = accuDelay(rxcntr-1,1:numdetectsat)/caChipRate;
bitsyncTime = (maxTransitionLocation(1:numdetectsat) - 1)*caCodeBlockDuration;
framesyncTime = (syncidx(1:numdetectsat)-1)*numCACodeBlocksPerBit*caCodeBlockDuration;

% Calculate transmission time from these parameters
tt1 = codeOffsetTime(:) - trackingOffsetTime(:) + bitsyncTime + framesyncTime;

%% Tracking loop results
figure(5)
plot(abs(accuFqy(:,1)))
