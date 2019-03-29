   clear all
   close all
   channel = [1 2];
   [dataFilename,pathName] = uigetfile('*.abf','Select your data file'); % e.g. B100519_0005
   [~,expName] = fileparts(dataFilename);
   [data,si,h] = abf2load([expName '.abf']); %si is the sampling interval (in us); h is the structure containing all info about the pClamp file
   curSeries_mV = squeeze(data(:,channel(1),:)); %This is the channel containing the spikes (intracellular or juxta)
   LFPseries_mV = squeeze(data(:,channel(2),:)); %This is the LFP channel
%     [curSeries_mV,si,h] = abf2load('C:\Users\vikash.choudhary\Desktop\scripts\vikash_intra_09.abf');
%     [LFPseries_mV,si1,h1] = abf2load('C:\Users\vikash.choudhary\Desktop\scripts\vikash_up_down_09.abf');
    sf = 1e+6*(1/si);
    pinkColor = [244 194 194]./255;
    medFilterWin = 10e-3;
    bandCutOff = [0.1 20];    % 0.1 Hz has been chosen to cancel drift effects in the baseline of the membrane potential [0.1 20]
    nBins = 100;
    prctileTh = 99;
    bandCutOffNorm = bandCutOff/(sf);
    %[b1,a1] = butter(2,0.1/sf,'high');
    [b,a] = butter(2,bandCutOffNorm,'bandpass');
    boundaryWin = 1/bandCutOff(1).*sf; % samples
    % UP & DOWN states detection from membrane potential
    stateDurTh_ms = 100;
    stateIntervTh_ms = 50;
    stateDurTh = stateDurTh_ms*1e-3*sf;
    stateIntervTh = stateIntervTh_ms*1e-3*sf;
    % LFP traces pre-processing
    sfLFP = 1000;   % LFP decimated to 1KHz
    lowPassCutOff = 500; % Hz
    lowPassCutOffNorm = lowPassCutOff/(sf*0.5);
    [bLow,aLow] = ellip(2,0.1,40,lowPassCutOffNorm);
    % PSD computation parameters
    nfft = 1024;
    welchWin = 2; %s
    welchWinSamples = welchWin.*sfLFP;
    welchOverlap = round(0.5*welchWinSamples);
    % Transition-triggered analysis of LFP
    halfTransTriggWin = 0.5*sfLFP; % samples
    % Frequency bands for LFP analysis
    fBands = [0 1;1 2;0 2;2 4;1 3];
    phaseMethod = 'hilbert';
    curSeriesMedFilt = medfilt1(curSeries_mV,medFilterWin*sf);  
    %% Low-pass 20 Hz
    %curSeriesLPFilt = filtfilt(b1,a1,curSeriesMedFilt);
    curSeriesLPFilt = filtfilt(b,a,curSeriesMedFilt);
    curSeriesLPFilt = curSeriesLPFilt(boundaryWin+1:end);
    nSamples = length(curSeriesLPFilt);
    %% Plot
    time = (1:length(curSeries_mV)).*1/sf;
    figure(1)
    hold all
    plot(time(boundaryWin+1:end), curSeriesLPFilt,'k');%  time(boundaryWin+1:end),
    xlabel('Time [s]')
    ylabel('Voltage [mV]')
    title('Intracellular membrane potential')
    %% Exclude upper 1 prctile 
    outlierTh = prctile(curSeriesLPFilt,prctileTh);
    data2Fit = curSeriesLPFilt(curSeriesLPFilt<=outlierTh);
    figure(2)
    [hAmpl,bins] = hist(data2Fit,nBins); 
    bar(bins,hAmpl./length(data2Fit));
    ylabel('# occurrences')
    xlabel('Voltage [mV]')
    hold all
    %% EM fitting by gaussians mixture
    options = statset('Display','final','MaxIter',500);
    try
        GMmodel = gmdistribution.fit(data2Fit,2,'Options',options);   % ill-posed??
    catch err
        %keyboard;
%         continue; %#ok<BRKCONT>
    end
    mu = GMmodel.mu;
    sigma = squeeze(GMmodel.Sigma);
    figure(2)
    yFit = pdf(GMmodel,bins');
    plot(bins,yFit,'r');
    %% Compute thresholds for UP & DOWN state detection
    [~,idxUP] = max(mu);
    [~,idxDOWN] = min(mu);
    thUP = mu(idxUP)-sqrt(sigma(idxUP));
    thDOWN = mu(idxDOWN)+sqrt(sigma(idxDOWN));
    figure(2)
    line([thUP thUP],[0 max(yFit)],'LineStyle','--','LineWidth',1,'Color','c')
    line([thDOWN thDOWN],[0 max(yFit)],'LineStyle','--','LineWidth',1,'Color','g')
    legend('Membrane potential hist','GM fit','UP state th','DOWN state th');
    %% Detect UP and DOWN states from intra-cellular trace and save results
    UP_states = curSeriesLPFilt>thUP;
    UP_states_start = find(diff(UP_states)==1);
    UP_states_end = find(diff(UP_states)==-1);
    % correct boundaries
    if UP_states_end(1)<UP_states_start(1)
        UP_states_end = UP_states_end(2:end);
    end
    if UP_states_start(end)>UP_states_end(end)
        UP_states_start = UP_states_start(1:end-1);
    end
    % join periods of UP-state (above threshold) that are separated by less
    % than a pre-defined time threshold
    UP_states_DET = [UP_states_start UP_states_end]; % samples
    UP_states_interval = UP_states_DET(2:end,1)-UP_states_DET(1:end-1,2);
    UP_states_2join = find(UP_states_interval<=stateIntervTh);
    UP_states_DET(UP_states_2join,2)=UP_states_DET(UP_states_2join+1,2);
    UP_states_DET(UP_states_2join+1,:)=[];
    % delete UP states shorter than a pre-defined time threshold
    UP_state_dur = UP_states_DET(:,2)-UP_states_DET(:,1);
    UP_states_OK = UP_state_dur>=stateDurTh;
    UP_states_DET = UP_states_DET(UP_states_OK,:);
    %
    DOWN_states = curSeriesLPFilt<thDOWN;
    DOWN_states_start = find(diff(DOWN_states)==1);
    DOWN_states_end = find(diff(DOWN_states)==-1);
    if DOWN_states_end(1)<DOWN_states_start(1)
        DOWN_states_end = DOWN_states_end(2:end);
    end
    if DOWN_states_start(end)>DOWN_states_end(end)
        DOWN_states_start = DOWN_states_start(1:end-1);
    end
    % join periods of DOWN-state (below threshold) that are separated by less
    % than a pre-defined time threshold
    DOWN_states_DET = [DOWN_states_start DOWN_states_end]; % samples
    DOWN_states_interval = DOWN_states_DET(2:end,1)-DOWN_states_DET(1:end-1,2);
    DOWN_states_2join = find(DOWN_states_interval<=stateIntervTh);
    DOWN_states_DET(DOWN_states_2join,2)=DOWN_states_DET(DOWN_states_2join+1,2);
    DOWN_states_DET(DOWN_states_2join+1,:)=[];
    % delete DOWN states shorter than a pre-defined time threshold
    DOWN_state_dur = DOWN_states_DET(:,2)-DOWN_states_DET(:,1);
    DOWN_states_OK = DOWN_state_dur>=stateDurTh;
    DOWN_states_DET = DOWN_states_DET(DOWN_states_OK,:);
    figure(1)
    for uu = 1:length(UP_states_DET)
        plot(time(boundaryWin+(UP_states_DET(uu,1):UP_states_DET(uu,2))),curSeriesLPFilt(UP_states_DET(uu,1):UP_states_DET(uu,2)),'Color',pinkColor,'LineWidth',2)
    end
    for dd = 1:length(DOWN_states_DET)
        plot(time(boundaryWin+(DOWN_states_DET(dd,1):DOWN_states_DET(dd,2))),curSeriesLPFilt(DOWN_states_DET(dd,1):DOWN_states_DET(dd,2)),'Color','b','LineWidth',2)
    end
    
    %% Pre-processing of LFP trace
    % Low-pass filtering of LFP
    LFPseries_mV = LFPseries_mV(boundaryWin+1:end).*1e+3;
    LFPdataFilt = filtfilt(bLow,aLow,LFPseries_mV);
    % Downsampling to 1 KHz
    nSamplesDec = nSamples/(sf/sfLFP);
    timeDec = (1/sfLFP).*(1:1:(nSamplesDec+1))+boundaryWin/sf;
    LFPdata_dec = decimate(LFPdataFilt,sf/sfLFP);
    % Plot LFP and UP & DOWN states
    UP_states_DET_dec = round(UP_states_DET./(sf/sfLFP));
    DOWN_states_DET_dec = round(DOWN_states_DET./(sf/sfLFP));
    figure(3)
    plot(timeDec', LFPdata_dec,'k');
    hold all
    for uu = 1:length(UP_states_DET)
        plot(timeDec(UP_states_DET_dec(uu,1):UP_states_DET_dec(uu,2)),LFPdata_dec(UP_states_DET_dec(uu,1):UP_states_DET_dec(uu,2)),'Color',pinkColor,'LineWidth',2)
    end
    for dd = 1:length(DOWN_states_DET)
        plot(timeDec(DOWN_states_DET_dec(dd,1):DOWN_states_DET_dec(dd,2)),LFPdata_dec(DOWN_states_DET_dec(dd,1):DOWN_states_DET_dec(dd,2)),'Color','b','LineWidth',2)
    end
    title('LFP')
    xlabel('Time')
    ylabel('Voltage [mV]')
    %% Compute and plot PSDs
    [psdLFP,fLFP] = pwelch(LFPdata_dec,welchWinSamples,welchOverlap,nfft,sfLFP);
    psdLFP_UP = periodogram(LFPdata_dec(UP_states_DET_dec(1,1):UP_states_DET_dec(1,2)),[],nfft,sfLFP);
    for uu = 2:size(UP_states_DET_dec,1)
        psdLFP_UP = psdLFP_UP + periodogram(LFPdata_dec(UP_states_DET_dec(uu,1):UP_states_DET_dec(uu,2)),[],nfft,sfLFP);
    end
    psdLFP_UP = psdLFP_UP./size(UP_states_DET_dec,1);
    psdLFP_DOWN = periodogram(LFPdata_dec(DOWN_states_DET_dec(1,1):DOWN_states_DET_dec(1,2)),[],nfft,sfLFP);
    for uu = 2:size(DOWN_states_DET_dec,1)
        psdLFP_DOWN = psdLFP_DOWN + periodogram(LFPdata_dec(DOWN_states_DET_dec(uu,1):DOWN_states_DET_dec(uu,2)),[],nfft,sfLFP);
    end
    psdLFP_DOWN = psdLFP_DOWN./size(DOWN_states_DET_dec,1);
    %
    figure(4)
    loglog(fLFP,psdLFP,'k');
    hold all
    loglog(fLFP,psdLFP_UP,'r');
    loglog(fLFP,psdLFP_DOWN,'b');
    xlabel('Frequency')
    ylabel('PSD')
    legend('LFP','UP states','DOWN states')
    figure(5)
    semilogx(psdLFP_UP./psdLFP_DOWN)
    xlabel('Frequency')
    xlim([1 200])
    ylabel('Ratio UP/DOWN')
    %% Compute transition-triggered LFP average
    LFPtriggAvg_UP = zeros(halfTransTriggWin*2,size(UP_states_DET_dec,1));
    UP_states_DET_dec(UP_states_DET_dec(:,1)<=halfTransTriggWin|UP_states_DET_dec(:,1)>(length(LFPdata_dec)-halfTransTriggWin),:)=[];
    for uu = 1:size(UP_states_DET_dec,1)
        LFPtriggAvg_UP(:,uu) = LFPdata_dec(UP_states_DET_dec(uu,1)-halfTransTriggWin:UP_states_DET_dec(uu,1)+halfTransTriggWin-1);
    end
    LFPtriggAvg_UP_mean = mean(LFPtriggAvg_UP,2);
    LFPtriggAvg_UP_std = std(LFPtriggAvg_UP,[],2);
    %
    LFPtriggAvg_DOWN = zeros(halfTransTriggWin*2,size(DOWN_states_DET_dec,1));
    DOWN_states_DET_dec(DOWN_states_DET_dec(:,1)<=halfTransTriggWin|DOWN_states_DET_dec(:,1)>(length(LFPdata_dec)-halfTransTriggWin),:)=[];
    for uu = 1:size(DOWN_states_DET_dec,1)
        LFPtriggAvg_DOWN(:,uu) = LFPdata_dec(DOWN_states_DET_dec(uu,1)-halfTransTriggWin:DOWN_states_DET_dec(uu,1)+halfTransTriggWin-1);
    end
    LFPtriggAvg_DOWN_mean = mean(LFPtriggAvg_DOWN,2);
    LFPtriggAvg_DOWN_std = std(LFPtriggAvg_DOWN,[],2);
    figure(6)
    timeWin = (-halfTransTriggWin:(halfTransTriggWin-1))./sfLFP;
    shadedErrorBar(timeWin,LFPtriggAvg_UP_mean,LFPtriggAvg_UP_std,{'LineWidth',2,'Color','r'},1);
    hold all
    shadedErrorBar(timeWin,LFPtriggAvg_DOWN_mean,LFPtriggAvg_DOWN_std,{'LineWidth',2,'Color','b'},1);
    line([0 0],ylim,'LineWidth',1,'Color','k')
    xlabel('Time [s]')
    ylabel('Voltage [mV]')

    %% Filter LFP in different frequency bands
    LFPdataDecFilt = zeros(length(LFPdata_dec),size(fBands,1));
    for bb = 1:size(fBands,1)
        if fBands(bb,1) == 0 % lowpass
            [b_LFPBandFilt,a_LFPBandFilt] = ellip(2,0.1,40,fBands(bb,2)./sfLFP,'low');
        else    % bandpass
            [b_LFPBandFilt,a_LFPBandFilt] = ellip(2,0.1,40,fBands(bb,:)./sfLFP);
        end        
        LFPdataDecFilt(:,bb) = filtfilt(b_LFPBandFilt,a_LFPBandFilt,LFPdata_dec);
    end
    %% For each frequency band, determine instantaneous amplitude and phase
    UP_states_signal_dec = false(size(LFPdata_dec,1),1);
    for uu = 1:length(UP_states_DET_dec)
        UP_states_signal_dec(UP_states_DET_dec(uu,1):UP_states_DET_dec(uu,2))=true;
    end
    DOWN_states_signal_dec = false(size(LFPdata_dec,1),1);
    for uu = 1:length(DOWN_states_DET_dec)
        DOWN_states_signal_dec(DOWN_states_DET_dec(uu,1):DOWN_states_DET_dec(uu,2))=true;
    end
    K = zeros(size(LFPdataDecFilt,1),size(fBands,1));
    phi = zeros(size(LFPdataDecFilt,1),size(fBands,1));
    for bb = 1:size(fBands,1)
        [K(:,bb),phi(:,bb)] = UP_DOWN_DET_compInstPhaseAmpl(LFPdataDecFilt(:,bb),phaseMethod);       
    end
    phi(phi<=0) = phi(phi<=0)+2*pi; % rescale from 0 to 2*pi
    phi = round(phi./(2*pi/360)).*(2*pi/360);
    %% Plot phase distribution functions during UP & DOWN states in the different frequency bands
    K_UP = K(UP_states_signal_dec,:);
    K_DOWN = K(DOWN_states_signal_dec,:);
    phi_UP = phi(UP_states_signal_dec,:);
    phi_DOWN = phi(DOWN_states_signal_dec,:);
    figure(7)
    for bb = 1:size(fBands,1)
        subplot(size(fBands,1),1,bb);
        [counts,centers] = hist(phi_UP(:,bb),360);
        bar(centers,counts,'EdgeColor',pinkColor,'FaceColor',pinkColor)
        hold all
        [counts,centers] = hist(phi_DOWN(:,bb),360);
        bar(centers,-counts,'EdgeColor','b','FaceColor','b')
        xlim([0 2*pi]);
    end
    title('Phase histograms')
    %% Differential likelihood of states
    phi_values = unique(phi);
    P_UP = zeros(length(phi_values),size(fBands,1));
    P_DOWN = zeros(length(phi_values),size(fBands,1));
    L_UP_DOWN = zeros(length(phi_values),size(fBands,1));
    L_fit = zeros(length(phi_values),size(fBands,1));
    theta = zeros(1,size(fBands,1));
    legendStrings = cell(1,size(fBands,1));
    for bb = 1:size(fBands,1)
        for ff = 1:length(phi_values)
            P_UP(ff,bb) = mean(UP_states_signal_dec(phi(:,bb)==phi_values(ff)));
            P_DOWN(ff,bb) = mean(DOWN_states_signal_dec(phi(:,bb)==phi_values(ff)));
        end
        L_UP_DOWN(:,bb) = smooth(P_UP(:,bb)-P_DOWN(:,bb));
        cosFit = fittype('cos(x-theta)');
        fitObj = fit(phi_values,L_UP_DOWN(:,bb),cosFit);
        L_fit(:,bb) = feval(fitObj,phi_values);
        theta(1,bb) = mod(fitObj.theta./(2*pi/360),360);
        legendStrings{bb} = [num2str(fBands(bb,1)),'-',num2str(fBands(bb,2)),' Hz'];
    end
    figure(8)
    plot(phi_values,L_UP_DOWN)
    xlim([0 2*pi])
    legend(legendStrings);
    xlabel('Phase [radians]')
    ylabel('L_X(t)')
   %% Save results
    saveas(1,fullfile(pathName,[expName,'_INTRA_DET_Vm_filtTrace.fig']),'fig');
    close(1)  
    saveas(2,fullfile(pathName,[expName,'_INTRA_DET_Vm_distrib.fig']),'fig');
    close(2)
    save(fullfile(pathName,[expName,'_INTRA_DET_results.mat']),'UP_states_DET','DOWN_states_DET','boundaryWin');
    saveas(3,fullfile(pathName,[expName,'_INTRA_DET_LFP_trace.fig']),'fig');
    close(3)
    saveas(4,fullfile(pathName,[expName,'_INTRA_DET_LFP_PSD.fig']),'fig');
    close(4)
    saveas(5,fullfile(pathName,[expName,'_INTRA_DET_LFP_ratioPSD.fig']),'fig');
    close(5)
    saveas(6,fullfile(pathName,[expName,'_INTRA_DET_LFP_transTrigg.fig']),'fig');
    close(6)
    saveas(7,fullfile(pathName,[expName,'_INTRA_DET_LFP_phaseHist.fig']),'fig');
    close(7)
    saveas(8,fullfile(pathName,[expName,'_INTRA_DET_likelihood.fig']),'fig');
    close(8)
    save(fullfile(pathName,[expName, '_theta_X.txt']),'theta','-ASCII');