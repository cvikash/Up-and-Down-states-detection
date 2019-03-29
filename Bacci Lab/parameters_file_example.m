% Common parameters
dataFilename = 'D:\Valentina\current_activities\Fellin\full_dataset\UP&DOWNstateDet\Interneuron Analysis\PV-positive\090114a\090114a.mat';
sf = 10000;                     % sampling frequency
sweepDur = 5.001;               % s
% Select sweeps to analyze for this file
selectWinFlag = true;           % TRUE if one wants to analyze a subset of sweeps, FALSE if one wants to analyze all sweeps
sweepSaveFlag = true;           % TRUE if one wants to save all figures 
seriesNumbers = [6 7 8 9];      % array containing series indexes to analyze
seriesSweeps = {[1:6,8:12,14,15,17,20],...
    [8,9,12,13,14,16,17,18,19,20,24,25,27:34,37:39,41,56,58:63,65,66,68,76:81,85,87,88]...
        [5,9,11:13,15,22:25,32,34,37:41,44:44,59,61,66,71:84,86,88,89]...
        [27,32,33,35:44,46,47,48,50:56,58,59,60,62,63,64,66,67,68,70:73,77,78,87]}; % cell array containing sweep indexes to analyze 
% LFP traces pre-processing
inverted = true;                % set to TRUE if LFP signal is inverted wrt to normal (UP state --> negative; DOWN state --> positive)
sfLFP = 1000;                   % LFP decimated to 1KHz
lowPassCutOff = 500;            % LFP low pass filtering cut-off [Hz]
lowPassCutOffNorm = lowPassCutOff/sf;
[bLow,aLow] = ellip(2,0.1,40,lowPassCutOffNorm);
% UP & DOWN states detection parameters
stateDurTh_ms = 100;            % minimum state duration
stateIntervTh_ms = 50;          % minimum state interval
nSTD_UP = 2;                    % number of standard deviations for UP state threshold identification
nSTD_DOWN = 2;                  % number of standard deviations for DOWN state threshold identification
stateDurTh = stateDurTh_ms*1e-3*sfLFP;
stateIntervTh = stateIntervTh_ms*1e-3*sfLFP;
% Frequency bands for LFP analysis
fBands = [0 1; 1 3];            % LFP frequency bands for UP state identification
fHigh = [20 100];               % high-frequency bandwidth (gamma range)
% Likelihood computation
phaseMethod = 'hilbert';        % method for instantaneous phase computation
if strcmp(phaseMethod,'hilbert')
    theta = [193.33289; 189.12116]./360*2*pi; % come from analysis of intracellular traces --> 0-1, 1-3 Hz
else if strcmp(phaseMethod,'interpol')
        theta = [192.05628; 189.01212]./360*2*pi; % come from analysis of intracellular traces --> 0-1, 1-3 Hz
    end
end
% Combined
combFlag = 1;                   % this parameter must be set to 1 if you want to include info about high frequencies
fAlphaBeta = [10 40];           % high-frequency range
winRmsFilt = 5*1e-3;            % parameters for processing of alpha-beta frequency range
winSmoothFilt = 50*1e-3;
win_samples = winRmsFilt*sfLFP;
winSmooth_samples = winSmoothFilt*sfLFP;
if mod(winSmooth_samples,2)==0
    winSmooth_samples = winSmooth_samples+1;
end
% Juxtasomal spike detection parameters
HPcutOff = 300;     % high-pass filtering cut-off
multCoeff = 6;      % number of standard deviations for spike detection in juxtasomal
% Transition-triggered analysis of LFP and spike data
halfTransTriggWin_s = 1;         % half-window [s] for transition-triggered analysis
halfTransTriggWin = halfTransTriggWin_s*sfLFP;      % samples
binSize_s = 25e-3;                  % bin size for transition-triggered AFR analysis
winSize = binSize_s*sf; % samples