[d,si,h] = abf2load('C:\Users\vikash.choudhary\Desktop\scripts\vikash_UP_DOWN_01.abf');
sf = 1000000*(1/si);
pinkColor = [244 194 194]./255;
Fp = 200; % cut off frequency
[bLow,aLow] = ellip(2,0.1,40,Fp/(0.5*sf));  % A lowpass filter at the cutoff frequncy 500 Hz
data = filtfilt(bLow,aLow,d);
Fs_new = 1000;
data = decimate(data,sf/Fs_new);  % Here we are downsampling the 20 KHz recordings to 1KHz
LFPdata = data(:); %38000:500000
LFPdata = -LFPdata;
selectedBands = [1 2;2 4];bands = {'1-2','2-4'}; % 1-2 and 2-4
phaseMethod = 'hilbert';
theta = [2.9931,2.9925];
%[4.11898,3.7699,3.5674,3.3432]
nSTD_UP = 4;
nSTD_DOWN = 4;
fHigh = [10 60];
combFlag = 0;
rmsWin = 5;
smoothWin = 51;
minDuration =200; 
minInterv = 50;
nfft = 1024;
welchWin = 2; %s
welchWinSamples = welchWin.*sf;
welchOverlap = round(0.5*welchWinSamples);

[S_comb,S_LFP,S_2,mu,sigma,thUP,thDOWN,GMmodel,L_t,Phi,k_high,k] = estimateUPDOWNdetTh(LFPdata,Fs_new,selectedBands,phaseMethod,bands,theta,nSTD_UP,nSTD_DOWN,fHigh,combFlag,rmsWin,smoothWin);
switch combFlag
    case 0
    variable = S_LFP;
    case 1
    variable = S_comb;
end
[UP_states_DET, DOWN_states_DET] = UP_DOWN_DET_detectStates(S_comb, thUP, thDOWN, minInterv, minDuration);

figure(1)
plot(S_comb)
hold on
plot([0 size(S_comb,1)], [thDOWN thDOWN], 'b')
plot([0 size(S_comb,1)], [thUP thUP], 'r')

figure(2);
hold all
plot(LFPdata, 'k')
for iUP=1: size(UP_states_DET, 1)
    line([UP_states_DET(iUP,1) UP_states_DET(iUP,2)], [max(LFPdata), max(LFPdata)], 'Color', 'red', 'LineWidth', 3)        
end

 for iDOWN = 1: size(DOWN_states_DET, 1)
     line([DOWN_states_DET(iDOWN,1) DOWN_states_DET(iDOWN,2)], [max(LFPdata), max(LFPdata)], 'Color', 'blue', 'LineWidth', 3)
 end

for iUP = 1: size(UP_states_DET, 1)
    x1 = UP_states_DET(iUP,1);
    x2 = UP_states_DET(iUP,2);
    x = [x1, x2];
    y1 = line([x1 x1],[min(LFPdata)-0.5 max(LFPdata)],'Color','red','LineWidth', 1);
    y2 = line([x2 x2],[min(LFPdata)-0.5 max(LFPdata)],'Color','r', 'LineWidth', 1);
    patch([x1 x2 x2 x1], [min(LFPdata)-0.5  min(LFPdata)-0.5 max(LFPdata) max(LFPdata)], [255/255,182/255,193/255],'FaceAlpha',0.3)
end


for iDOWN = 1: size(DOWN_states_DET, 1)
    x1 = DOWN_states_DET(iDOWN,1);
    x2 = DOWN_states_DET(iDOWN,2);
    x = [x1, x2];
    y1 = line([x1 x1],[min(LFPdata)-0.5 max(LFPdata)],'Color','b','LineWidth', 1);
    y2 = line([x2 x2],[min(LFPdata)-0.5 max(LFPdata)],'Color','b', 'LineWidth', 1);
    patch([x1 x2 x2 x1], [min(LFPdata)-0.5  min(LFPdata)-0.5 max(LFPdata) max(LFPdata)], [130/255 200/255 250/255],'FaceAlpha',0.3)
end
UP_state_duration = UP_states_DET(:,2)-UP_states_DET(:,1);
figure
hist(UP_state_duration)
xlabel('UP state duration (ms)');

[psLFPdata,fLFP] = pwelch(LFPdata,welchWinSamples,welchOverlap,nfft,Fs_new);
psLFPdata_UP = periodogram(LFPdata(UP_states_DET(1,1):UP_states_DET(1,2)),[],nfft,Fs_new);
for uu = 2:size(UP_states_DET,1)
    psLFPdata_UP = psLFPdata_UP + periodogram(LFPdata(UP_states_DET(uu,1):UP_states_DET(uu,2)),[],nfft,Fs_new);
end
psLFPdata_UP = psLFPdata_UP./size(UP_states_DET,1);
psLFPdata_DOWN = periodogram(LFPdata(DOWN_states_DET(1,1):DOWN_states_DET(1,2)),[],nfft,Fs_new);
for uu = 2:size(DOWN_states_DET,1)
    psLFPdata_DOWN = psLFPdata_DOWN + periodogram(LFPdata(DOWN_states_DET(uu,1):DOWN_states_DET(uu,2)),[],nfft,Fs_new);
end
psLFPdata_DOWN = psLFPdata_DOWN./size(DOWN_states_DET,1);
%
figure(4)
loglog(fLFP,psLFPdata,'k');
hold all
loglog(fLFP,psLFPdata_UP,'r');
loglog(fLFP,psLFPdata_DOWN,'b');
xlabel('Frequency')
ylabel('PSD')
legend('LFP','UP states','DOWN states')
figure(5)
semilogx(psLFPdata_UP./psLFPdata_DOWN)
xlabel('Frequency')
xlim([1 200])
ylabel('Ratio UP/DOWN')
