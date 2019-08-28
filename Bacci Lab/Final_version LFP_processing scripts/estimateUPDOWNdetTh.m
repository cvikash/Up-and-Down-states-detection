function [S_comb,S_LFP,S_2,mu,sigma,thUP,thDOWN,GMmodel] = estimateUPDOWNdetTh(LFPdata,sf,selectedBands,phaseMethod,bands,theta,nSTD_UP,nSTD_DOWN,fHigh,combFlag,rmsWin,smoothWin)
% estimateUPDOWNdetTh.m takes the LFP signal as input and returns decision
% variables and thresholds for up/down state detection.
% License: V. Pasquale, 2018, Optical Approaches to Brain Function, CC BY 4.0
% ----Input arguments----:
% *     LFPdata: column vector, time series of LFP, preferably sampled at 1 KHz and
%           low-pass filtered below 500 Hz, size [numSamples x 1]
% *     sf: sampling frequency in Hz (e.g. 1000 Hz)
% *     selectedBands: matrix of selected frequency bands (in Hz), size [numBands x 2], e.g. [0, 1; 1, 2];
% *     phaseMethod: 'hilbert' (preferred) or 'interpol'
% *     bands: available bands (in Hz), cell arrays of strings, e.g. 1×4 cell array {'0-1','1-2','0-2','2-4'}
% *     theta: vector of angles (in radians, not degrees) per each
%           frequency band, size 1 x numBands, e.g. [3.1148,2.9931,3.0941,2.9925]
% *     nSTD_UP: scalar, number of standard deviations for estimating up state
%           detection threshold, e.g. 2
% *     nSTD_DOWN: scalar, number of standard deviations for estimating down state
%           detection threshold, e.g. 2
% *     fHigh: high-frequency bandwidth (in Hz), size 1 x 2, e.g. [10, 50] or [20, 100] Hz
% *     combFlag: 0 or 1, 0 = only Saleem, 1 = Saleem + Mukovski
% *     rmsWin: root mean square filter (see Mukovski) window, in samples at the sampling frequency used, e.g. 5 at 1 Khz = 5 ms
% *     smoothWin: smooth filter (see Mukovski) window, in samples at the
%           sampling frequency used, e.g. 51 at 1 Khz = 51 ms (better if an odd
%           number)
% ----Output arguments----:
% *     S_comb: combined decision variable, Saleem + Mukovski, A.U.
% *     S_LFP: Saleem decision variable, A.U.
% *     S_2: Mukovski decision variable, A.U.
% *     mu: means of fitted Gaussians
% *     sigma: standard deviations of fitted Gaussians
% *     thUP: threshold for up state
% *     thDOWN: threshold for down state:
% *     GMmodel: Gaussian mixture fitted model
%% Compute k_high (for S_LFP computation)
[b_high,a_high] = ellip(2,0.1,40,fHigh./(0.5*sf));
LFPdataHigh = filtfilt(b_high,a_high,LFPdata);
[k_high,~] = UP_DOWN_DET_compInstPhaseAmpl(LFPdataHigh,'hilbert');
%% Compute MUA
if combFlag==1
    %% Mukovski method
    rmsLFPdata = rmsFilt(LFPdataHigh,rmsWin);
    smoothLFPdata = smooth(rmsLFPdata,smoothWin);
    S_2 = (smoothLFPdata-min(smoothLFPdata))./(max(smoothLFPdata)-min(smoothLFPdata));
    S_2_p95 = prctile(S_2,95);
    S_2(S_2 > S_2_p95) = S_2_p95;
else
    S_2 = [];
end
%% Filter in the different frequency bands
LFPdataDecFilt = zeros(length(LFPdata),size(selectedBands,1));
for bb = 1:size(selectedBands,1)
    if selectedBands(bb,1) == 0 % lowpass
        [b_LFPBandFilt,a_LFPBandFilt] = ellip(2,0.1,40,selectedBands(bb,2)./(0.5*sf),'low');
    else    % bandpass
        [b_LFPBandFilt,a_LFPBandFilt] = ellip(2,0.1,40,selectedBands(bb,:)./(0.5*sf));
    end
    LFPdataDecFilt(:,bb) = filtfilt(b_LFPBandFilt,a_LFPBandFilt,LFPdata);
end
%% For each frequency band, approximate L_t
k = zeros(size(LFPdataDecFilt,1),size(selectedBands,1));
phi = zeros(size(LFPdataDecFilt,1),size(selectedBands,1));
L_t = zeros(size(LFPdataDecFilt,1),size(selectedBands,1));
for bb = 1:size(selectedBands,1)
    [k(:,bb),phi(:,bb)] = UP_DOWN_DET_compInstPhaseAmpl(LFPdataDecFilt(:,bb),phaseMethod);
    phi(phi(:,bb)<=0,bb) = phi(phi(:,bb)<=0,bb)+2*pi; % rescale from 0 to 2*pi
    phi(:,bb) = round(phi(:,bb)./(2*pi/360)).*(2*pi/360); % round to 1 degree resolution
    bandStr = [num2str(selectedBands(bb,1)),'-',num2str(selectedBands(bb,2))];
    currTheta = theta(strcmp(bands,bandStr));
    L_t(:,bb) = cos(phi(:,bb)-currTheta);
end
%% Derive S_LFP
K = k./repmat(k_high+sum(k,2),1,size(selectedBands,1));
S_LFP = (sum(K.*L_t,2)+1)./2;
%% Compute the combined signal between S_LFP and S_MUA
if combFlag==1
    S_comb = (S_LFP+S_2)./2;
else
    S_comb = S_LFP;
end
S_comb = (S_comb-min(S_comb))./(max(S_comb)-min(S_comb));
%% Fit distribution of S_comb by GM
options = statset('Display','final','MaxIter',500);
try
    GMmodel = gmdistribution.fit(S_comb(S_comb<=prctile(S_comb,95)),3,'Start','plus','Options',options);
catch err
    keyboard;
    %     continue;
end
mu = GMmodel.mu;
sigma = squeeze(GMmodel.Sigma);
%% Compute thresholds for UP & DOWN state detection
[~,idxUP] = max(mu);
[~,idxDOWN] = min(mu);
thUP = mu(idxUP)-nSTD_UP*sqrt(sigma(idxUP));
thDOWN = mu(idxDOWN)+nSTD_DOWN*sqrt(sigma(idxDOWN));
end