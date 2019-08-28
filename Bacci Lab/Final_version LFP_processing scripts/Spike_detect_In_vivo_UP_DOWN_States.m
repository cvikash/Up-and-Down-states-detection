%Spike_detect_In_vivo_UP_DOWN_states.m takes detetcted threshold for UP and
%DOWN states, LFP and juxtacellular recordings, minimum period for states
%and between states, bandwidths to use for detection and which methods to
%use
warning off all
clear all
close all
clc

% Load pClamp Files into MATLAB and plot and define important parameters
channel = [1 2];
threshold = 1.5; %Set the threshold to detect spikes
peakdistance = 50; 
[dataFilename,pathName] = uigetfile('*.abf','Select your data file'); % e.g. B100519_0005
[~,expName] = fileparts(dataFilename);
filename = dataFilename;
[data,si,h] = abf2load([expName '.abf']); %si is the sampling interval (in us); h is the structure containing all info about the pClamp file
SPKs = squeeze(data(:,channel(1),:)); %This is the channel containing the spikes (intracellular or juxta)
LFP = squeeze(data(:,channel(2),:)); %This is the LFP channel
time = (0:si:(length(data)-1)*si)*1e-6; % time in seconds
numchannels = size(channel, 2);
Fs=1/si*1e+6; %Original sampling frequncy in Hz
dt = si*1e-6;     % Original sampling interval [s]
Nsamples = size(LFP,1);    % Original lenght of the waveform, in points
T = Nsamples*dt;      % Original length of the waveform, in s
bandwidth_LFP = input('Please enter the LFP bandwidths from the avalaible bands (e.g. [0 1; 0 2]): '); %these bands can be extracted from sample intacellular recordings analyzed with vikash_intra_UP_DOWN.m
minInterv = input ('Please enter the minimum interval between states (in samples; e.g. 50): ');
minDuration = input ('Please enter the minimum duration of a states (in samples; e.g. 100): ');

%Highpass filter the Spike signal to remove drift
Fc  = 100;                  % Higher cut-off frequency [Hz]
[bspks,aspks] = butter(8, Fc/(Fs/2.), 'high');
fSPKs = filtfilt(bspks,aspks,double(SPKs)); %do this filtering by the filter that we just built

% Highpass filter the LFP signal from for MUAs 
MUA_H  = 2000;                 % Higher cut-off frequency [Hz]
[MUA_b,MUA_a] = butter(8, MUA_H/(Fs/2.), 'high');
MUA = filtfilt(MUA_b,MUA_a,double(LFP)); %do the highpass filtering
invMUA = -(MUA);


% Let's decimate the LFP signal from Fs to Fs_dec (Fs_dec < Fs) and filter
% the signal below 200 Hz
Fp = 200; % cut off frequency
[bLow,aLow] = ellip(2,0.1,40,Fp/(0.5*Fs));  % A lowpass filter at the cutoff frequncy 200 Hz
data = filtfilt(bLow,aLow,LFP);
data = data-mean(data);
Fs_dec  = 1000;         % [Hz] New sampling rate
R  = Fs / Fs_dec;  % This is needed for the function decimate()
dLFP = decimate(data, R); % From Fs Hz to Fs_dec Hz
si_dec = (1./Fs_dec);        % Sampling interval, after decimation
dtime = 0:si_dec:(T-si_dec); % new time axis [s]   ...

%This should fix the extra sample bug 
if size(dtime,2)>size(dLFP,1)
    dtime = 0:si_dec:(T-si_dec);
    elseif size(dtime,2)<size(dLFP,1)
    dtime = 0:si_dec:T;
end

%Let's band-pass filter the decimated LFP from LOW to HIGH 
LOW   = 20.;                  % Lower cut-off frequency [Hz]
HIGH  = 100.;                 % Higher cut-off frequency [Hz]
[b,a] = butter(8, [LOW/(Fs_dec/2.) HIGH/(Fs_dec/2.)]);
filtered_LFP = filtfilt(b,a,double(dLFP)); %do this filtering by the filter that we just built

% Detect spikes according to the threshold, and visualize detected spikes
figure(1)
A = subplot(411);
plot(time, fSPKs)
[pks, locs]=findpeaks(fSPKs, 'MinPeakHeight', threshold, 'MinPeakDistance', peakdistance); %Detect spikes in the intra or juxta channel using findpeaks
peaktimes = locs.*dt; %Times of spikes (in s)
hold on
scatter(peaktimes, pks,'r');%Place a red marlker where each spike is found
plot([0 T], [threshold threshold], 'r') %visualzie the threshold with a red line
xlabel('Time (s)');
ylabel('mV');

B = subplot(412);
plot(dtime, filtered_LFP, 'k') %Plot the decimated and filtered LFP (in the 20-100 Hz range)
xlabel('Time (s)');
ylabel('Filtered LFP (mV)');

C = subplot(413);
plot(time, LFP, 'r') %Plot the unfiltered LFP
xlabel('Time (s)');
ylabel('LFP (mV)');

D = subplot(414);
% build the spectrogram 
movingwin       = [0.5 0.05]; % Size [s] and overlap [s] of the moving window
params.Fs       = 1000;       % Sampling rate [Hz] of the (decimated) LFP
params.fpass    = [10 100];    % Range of frequency to plot [Hz]
params.tapers   = [5 9];      % Multi-taper mysterious parameters 
params.trialave = 0;          % Whether (1) or not (0) to perform an average across all trials
params.err      = 0;	      % Whether (1) or not (0) produce the statistical error in the output

odLFP = dLFP - mean(dLFP);    % Remove offset 
[S1,t,f]= mtspecgramc(odLFP, movingwin, params);
plot_matrix(S1,t,f); 
colormap jet;
colorbar off;
linkaxes([A, B, C, D], 'x');

figure(2)
plot_matrix(S1,t,f); 
colormap jet;

spikefreq = length(pks)/(time(end));
x = peaktimes(2:end);
ISI = x-peaktimes(1:end-1);
instfreq = 1./ISI;

%display LFP and MUAs
figure(3)
MA = subplot(211);
plot(time, LFP, 'r')
xlabel('time')
ylabel('LFP')
MB = subplot(212);
plot(time, MUA, 'k')
xlabel('time')
ylabel('MUA')
linkaxes([MA, MB], 'x');

%collect all data in a structure
Result_Spikes=struct;
Result_Spikes.Peak = pks;
Result_Spikes.PeakTime = peaktimes;
Result_Spikes.SpikeFreq = spikefreq;
Result_Spikes.InterSpikeInterval = ISI;
Result_Spikes.InstFrequency = instfreq;

save([filename '_preprocessed.mat'],'dt', 'LFP', 'peaktimes');

figure(4)
plot(time, SPKs)
hold on
plot(time, fSPKs)
xlabel('time')
ylabel('SPKs & fSPKS ')

%% Analyze UP and DOWN states (from Valentina Pasquale and Tommaso Fellin)
sf = Fs_dec;
selectedBands = bandwidth_LFP;
phaseMethod = 'hilbert';
bands = {'1-2','2-4'};
theta = [3.1048,3.0931];
nSTD_UP = 2;
nSTD_DOWN = 2;
fHigh = [20, 100];
combFlag = 0;
rmsWin = 5;
smoothWin = 51;
nfft = 1024;
welchWin = 2; %s
welchWinSamples = welchWin.*sf;
welchOverlap = round(0.5*welchWinSamples);
nsamples_dec = length(dLFP);

[S_comb,S_LFP,S_2,mu,sigma,thUP,thDOWN,GMmodel] = estimateUPDOWNdetTh(dLFP,sf,selectedBands,phaseMethod,bands,theta,nSTD_UP,nSTD_DOWN,fHigh,combFlag,rmsWin,smoothWin); %Set the threshold for UP and DOWN state detection

[UP_states_DET, DOWN_states_DET] = UP_DOWN_DET_detectStates(S_LFP, thUP, thDOWN, minInterv, minDuration); %Detect UP and DOWN states

figure(5)
plot(S_comb)
xlabel('time')
ylabel('S comb (evidence variable)')
legend('S_comb','thUP','thDOWN')
hold on
plot([0 size(S_comb,1)], [thDOWN thDOWN], 'b')
plot([0 size(S_comb,1)], [thUP thUP], 'r')
%% LFP with Up and Down states
%Plot the decimated LFP and mark the beginning and end of each UP and DOWN state. 
%Red and green horizontal bars are UP and DOWN states, respectively
nsamples = nsamples_dec*(Fs/Fs_dec);
% Plot LFP and UP & DOWN states
UP_states_DET_org = round(UP_states_DET.*(Fs/Fs_dec));
DOWN_states_DET_org = round(DOWN_states_DET.*(Fs/Fs_dec));
time_org = 1/Fs.*(1:1:nsamples);
figure(6)
plot(data, 'k')
xlabel('time')
ylabel('LFP (orignal)')
hold on

for iUP = 1: size(UP_states_DET_org, 1)
    x1 = UP_states_DET_org(iUP,1);
    x2 = UP_states_DET_org(iUP,2);
    x = [x1, x2];
    y1 = line([x1 x1],[min(data)-0.5 max(data)],'Color','red','LineWidth', 1);
    y2 = line([x2 x2],[min(data)-0.5 max(data)],'Color','r', 'LineWidth', 1);
    patch([x1 x2 x2 x1], [min(data)-0.5  min(data)-0.5 max(data) max(data)], [255/255,182/255,193/255],'FaceAlpha',0.3)
end


for iDOWN = 1: size(DOWN_states_DET_org, 1)
    x1 = DOWN_states_DET_org(iDOWN,1);
    x2 = DOWN_states_DET_org(iDOWN,2);
    x = [x1, x2];
    y1 = line([x1 x1],[min(data)-0.5 max(data)],'Color','b','LineWidth', 1);
    y2 = line([x2 x2],[min(data)-0.5 max(data)],'Color','b', 'LineWidth', 1);
    patch([x1 x2 x2 x1], [min(data)-0.5  min(data)-0.5 max(data) max(data)], [130/255 200/255 250/255],'FaceAlpha',0.3)
end


%% Here we can see the distribution of UP state durations
UP_state_duration = UP_states_DET(:,2)-UP_states_DET(:,1);
figure(7)
hist(UP_state_duration)
xlabel('UP state duration (ms)');

%% Compute and plot PSDs
[psdLFP,fLFP] = pwelch(dLFP,welchWinSamples,welchOverlap,nfft,sf);
psdLFP_UP = periodogram(dLFP(UP_states_DET(1,1):UP_states_DET(1,2)),[],nfft,sf);
for uu = 2:size(UP_states_DET,1)
    psdLFP_UP = psdLFP_UP + periodogram(dLFP(UP_states_DET(uu,1):UP_states_DET(uu,2)),[],nfft,sf);
end
psdLFP_UP = psdLFP_UP./size(UP_states_DET,1);
psdLFP_DOWN = periodogram(dLFP(DOWN_states_DET(1,1):DOWN_states_DET(1,2)),[],nfft,sf);
for uu = 2:size(DOWN_states_DET,1)
    psdLFP_DOWN = psdLFP_DOWN + periodogram(dLFP(DOWN_states_DET(uu,1):DOWN_states_DET(uu,2)),[],nfft,sf);
end
psdLFP_DOWN = psdLFP_DOWN./size(DOWN_states_DET,1);
%
figure(8)
loglog(fLFP,psdLFP,'k');
hold all
loglog(fLFP,psdLFP_UP,'r');
loglog(fLFP,psdLFP_DOWN,'b');
xlabel('Frequency')
ylabel('PSD')
legend('LFP','UP states','DOWN states')
figure(9)
semilogx(psdLFP_UP./psdLFP_DOWN)
xlabel('Frequency')
xlim([1 200])
ylabel('Ratio UP/DOWN')

%% Detect and quantify spikes in each UP states
UP_states_DET_sec = UP_states_DET_org.*(si*1e-6);
DOWN_states_DET_sec = DOWN_states_DET_org.*(si*1e-6);
Spikes_UP_states = struct;
figure(10)
UPa = subplot(211);
plot(time,fSPKs)
for iUP = 1: size(UP_states_DET_org, 1)
    x1 = UP_states_DET_sec(iUP,1);
    x2 = UP_states_DET_sec(iUP,2);
    x = [x1, x2];
    y1 = line([x1 x1],[min(fSPKs) max(fSPKs)-2],'Color','red','LineWidth', 1);
    y2 = line([x2 x2],[min(fSPKs) max(fSPKs)-2],'Color','r', 'LineWidth', 1);
    patch([x1 x2 x2 x1], [min(fSPKs)  min(fSPKs) max(fSPKs)-2 max(fSPKs)-2], [255/255,182/255,193/255],'FaceAlpha',0.3)
end

for iDOWN = 1: size(DOWN_states_DET_sec, 1)
    x1 = DOWN_states_DET_sec(iDOWN,1);
    x2 = DOWN_states_DET_sec(iDOWN,2);
    x = [x1, x2];
    y1 = line([x1 x1],[min(fSPKs) max(fSPKs)-2],'Color','b','LineWidth', 1);
    y2 = line([x2 x2],[min(fSPKs) max(fSPKs)-2],'Color','b', 'LineWidth', 1);
    patch([x1 x2 x2 x1], [min(fSPKs)  min(fSPKs) max(fSPKs)-2 max(fSPKs)-2], [130/255 200/255 250/255],'FaceAlpha',0.3)
end
for iUP=1: size(UP_states_DET, 1)
    peaktimes_UP = peaktimes(peaktimes>UP_states_DET_sec(iUP,1) & peaktimes<UP_states_DET_sec(iUP,2)); %get the spike times (detected above) within each UP state
    pks_UP = pks(find(peaktimes_UP)); %get the peak of the each spike within each UP state (used to place a marker in the figure)
    
    hold on
    scatter(peaktimes_UP, pks_UP,'r');%Place a red marlker where each spike is found
    plot([0 T], [threshold threshold], 'r') %visualzie the threshold with a red line
    xlabel('Time (s)');
    ylabel('fSPKS (mV)');
    Spikes_UP_states(iUP).PeakTimes_UP_states = peaktimes_UP;
end

UPb = subplot(212);
plot(time, LFP, 'k')
for iUP = 1: size(UP_states_DET_org, 1)
    x1 = UP_states_DET_sec(iUP,1);
    x2 = UP_states_DET_sec(iUP,2);
    x = [x1, x2];
    y1 = line([x1 x1],[min(LFP)-0.5 max(LFP)],'Color','red','LineWidth', 1);
    y2 = line([x2 x2],[min(LFP)-0.5 max(LFP)],'Color','r', 'LineWidth', 1);
    patch([x1 x2 x2 x1], [min(LFP)-0.5  min(LFP)-0.5 max(LFP) max(LFP)], [255/255,182/255,193/255],'FaceAlpha',0.3)
end


for iDOWN = 1: size(DOWN_states_DET_sec, 1)
    x1 = DOWN_states_DET_sec(iDOWN,1);
    x2 = DOWN_states_DET_sec(iDOWN,2);
    x = [x1, x2];
    y1 = line([x1 x1],[min(LFP)-0.5 max(LFP)],'Color','b','LineWidth', 1);
    y2 = line([x2 x2],[min(LFP)-0.5 max(LFP)],'Color','b', 'LineWidth', 1);
    patch([x1 x2 x2 x1], [min(LFP)-0.5  min(LFP)-0.5 max(LFP) max(LFP)], [130/255 200/255 250/255],'FaceAlpha',0.3)
end
xlabel('Time (s)');
ylabel('LFP (mV)');
linkaxes([UPa, UPb], 'x');
%% Back calculation of spike numbers from detected UP states
%Calculate the number of spikes in each UP state
for i = 1:size(Spikes_UP_states, 2)
Nspks_per_UPstate (i) = numel(Spikes_UP_states(i).PeakTimes_UP_states);
end

Nspks_per_UPstate = Nspks_per_UPstate'; %output of the number of spikes per UP State
edges = [0:1:max(Nspks_per_UPstate)];
figure(11)
subplot(121)
hist(Nspks_per_UPstate, edges)
xlabel('Number of spikes per UP state');
ylabel('count');
subplot(122)
avgNumSpikesUPstate = mean(Nspks_per_UPstate); %average number of spikes per UP state
stdevSpikesUPstate = std(Nspks_per_UPstate); %standard deviation of the mean number of spikes per UP state
bar(avgNumSpikesUPstate)
hold on
errorbar(avgNumSpikesUPstate, stdevSpikesUPstate)

%% phase analysis
arg = input("Do you want to continue the monte carlo analysis?\n yes\n no\n", 's');
switch arg
    case 'yes'
        disp("Continuing the execution of the script :-) ")
    case 'no'
        return 
end
postprocessing_ppc_and_phaseHistogram_2;
save([filename '_Spikes_Phase.mat'],'Result_Spikes', '-append');






