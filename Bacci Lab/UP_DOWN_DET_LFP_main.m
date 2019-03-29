% Script containing code to test automatic UP/DOWN state detection algorithm's performance
% on simultaneous patch clamp on cortical pyramidal neurons and local field potential recording.
% It relies on the parameters estimated by running UP_DOWN_DET_intra_main.m.
% Written by Valentina Pasquale (2017)
% Contact: valentina.pasquale@iit.it
%% Clear workspace
clear all
close all
%% Parameters
sf = 10000;     % sampling frequency
sweepDur = 5.001; % s
selectWinFlag = {true,true,false,false,true,true};
seriesNumbers = {[3];[2 3];[];[];[5];[2 3]};
seriesSweeps = {{[1:20 40:80 100:140]};{[0:40 50:80],[0:54 70:100 110:140]};{[]};{[]};{[0:40 80:120]};{[20:140], [1:122]}};
joinSeries = 1; % if this flag == 1, all series are joined in a unique trace; 
pinkColor = [244 194 194]./255; % fixed
% LFP traces pre-processing
sfLFP = 1000;   % LFP decimated to 1KHz
lowPassCutOff = 500; % Hz
lowPassCutOffNorm = lowPassCutOff/sf;
[bLow,aLow] = ellip(2,0.1,40,lowPassCutOffNorm);
% UP & DOWN states detection parameters
stateDurTh_ms = 100;    % 100 ms
stateIntervTh_ms = 50;  % 50 ms
nSTD_UP = 2;
nSTD_DOWN = 2;
stateDurTh = stateDurTh_ms*1e-3*sfLFP;
stateIntervTh = stateIntervTh_ms*1e-3*sfLFP;
% Frequency bands for LFP analysis
fBands = [0 1;1 3];
fHigh = [20 100];
% Likelihood computation
theta = [193.33289; 189.12116]./360*2*pi; % for hilbert: come from analysis of intracellular traces --> 0-1, 1-3 Hz
nBins = 100;
phaseMethod = 'hilbert';
% phaseMethod = 'interpol';
% Combined
combFlag = 1;   % this parameter must be set to 1 if you want to include info about high frequencies
fAlphaBeta = [10 40];
winRmsFilt = 5*1e-3;
winSmoothFilt = 50*1e-3;
win_samples = winRmsFilt*sfLFP;
winSmooth_samples = winSmoothFilt*sfLFP;
if mod(winSmooth_samples,2)==0
    winSmooth_samples = winSmooth_samples+1;
end
% Transition-triggered analysis of LFP
halfTransTriggWin = 0.5*sfLFP; % samples
% String to identify results
str = 'combined_0-1_1-3Hz_nSTD2_hilbert_gauss_20160301';
% Folders
parameters_struct = struct('stateDurationThreshold',stateDurTh_ms,'stateMinIntervalThreshold',stateIntervTh_ms,...
    'multCoeff_UP',nSTD_UP,'multCoeff_DOWN',nSTD_DOWN,'frequencyBands',fBands,...
    'alphaBetaBand',fAlphaBeta,'winRmsFilt',winRmsFilt*1e+3,'winSmoothFilt',winSmoothFilt*1e+3,...
    'phaseMethod',phaseMethod,'combFlag',combFlag,'stringResults',str);
%% D:\Valentina\current_projects\Fellin\dati_stefano_2photons\Patch+LFP\ must be changed
ROC_analysis_path = ['D:\Valentina\current_activities\Fellin\dati_stefano_2photons\Patch+LFP\ROC_analysis_',str];
perfChosenTh_path = ['D:\Valentina\current_activities\Fellin\dati_stefano_2photons\Patch+LFP\perfChosenTh_',str];
CoIn_path = ['D:\Valentina\current_activities\Fellin\dati_stefano_2photons\Patch+LFP\CoIn_',str];
%% D:\Valentina\current_projects\Fellin\dati_stefano_2photons\Patch+LFP\ must be changed with your path
sourceFolders = {...
    'D:\Valentina\current_activities\Fellin\dati_stefano_2photons\Patch+LFP\130214a\130214a.mat',...
    'D:\Valentina\current_activities\Fellin\dati_stefano_2photons\Patch+LFP\130214c\130214c.mat',...
    'D:\Valentina\current_activities\Fellin\dati_stefano_2photons\Patch+LFP\250314a\250314a.mat',...
    'D:\Valentina\current_activities\Fellin\dati_stefano_2photons\Patch+LFP\250314b\250314b.mat',...
    'D:\Valentina\current_activities\Fellin\dati_stefano_2photons\Patch+LFP\280314b\280314b.mat',...
    'D:\Valentina\current_activities\Fellin\dati_stefano_2photons\Patch+LFP\280314c\280314c.mat'};
for ss = 1:length(sourceFolders)
    %% Load data
    [pathName,dataFilename] = fileparts(sourceFolders{ss});
    [~,expName] = fileparts(dataFilename);
    saveFolderName = [pathName,filesep,expName,'_',datestr(now,'yyyymmdd_HHMM')];
    try
        if ~isdir(saveFolderName)
            mkdir(saveFolderName);
        else
            rmdir(saveFolderName,'s');
            mkdir(saveFolderName);
        end
        if ~isdir(ROC_analysis_path)
            mkdir(ROC_analysis_path);
        end
        if ~isdir(perfChosenTh_path)
            mkdir(perfChosenTh_path);
        end
        if ~isdir(CoIn_path)
            mkdir(CoIn_path);
        end
    catch ME
        errordlg(ME.message,ME.identifier)
        return
    end
    %% Join 1st and 2nd trace of the same sweep in a unique variable, join sweeps of the same series in a unique variable
    if selectWinFlag{ss} == true
        fullTrace = UP_DOWN_DET_organizeData(fullfile(pathName,[dataFilename,'.mat']),seriesNumbers{ss},seriesSweeps{ss});
    else
        [fullTrace,seriesNumbers{ss},~] = UP_DOWN_DET_organizeData(fullfile(pathName,dataFilename));
    end
    if isempty(fullTrace{1})
        return
    end
    %% Analyze all sweeps
    if joinSeries == 1
        fullTrace = cell2mat(fullTrace);
%         seriesNumbers{ss} = ['_',strrep(num2str(seriesNumbers{ss}),' ','_'),'_'];
        nSeries = 1;
    else
        nSeries = length(fullTrace);
    end
    for ii = 1:nSeries
        if joinSeries == 1
            curSeries = fullTrace;
            seriesText = ['_',strrep(num2str(seriesNumbers{ss}),' ','_'),'_'];
        else
            curSeries = fullTrace{ii};
            seriesText = num2str(seriesNumbers{ss}(ii));
        end          
        % 1st col: time
        % 2nd col: Trace 1 (intra- o juxta-cellular)
        % 3rd col: Trace 2 (LFP)
        %% Pre-processing of LFP trace
        % Low-pass filtering of LFP
        LFPseries_mV = curSeries(:,3).*1e+3;    % 3rd column
        LFPdataFilt = filtfilt(bLow,aLow,LFPseries_mV); % LP filter, cut-off 500 Hz
        % Downsampling to 1 KHz
        nSamples = length(LFPdataFilt);
        nSamplesDec = nSamples/(sf/sfLFP);
        timeDec = (1/sfLFP).*(1:1:nSamplesDec);
        LFPdata_dec = decimate(LFPdataFilt,sf/sfLFP);
        % Filter LFP in the delta frequency range
        [b_deltaFilt,a_deltaFilt] = ellip(2,0.1,40,[0.1 4]./sfLFP);
        deltaFiltLFP = filtfilt(b_deltaFilt,a_deltaFilt,LFPdata_dec);
        % Plot raw signals and synchronization index
        figure(1)
        subplot(3,1,1)
        intraSeries_mV = curSeries(:,2).*1e+3;    % 2nd column
        time = (1:length(intraSeries_mV)).*1/sf;
        plot(time,intraSeries_mV,'k');
        hold all
        xlim([0 max(time)])
        xlabel('Time [s]')
        ylabel('Voltage [mV]')
        title('Intracellular membrane potential')
        subplot(3,1,2);
        plot(timeDec,LFPdata_dec,'k');
        xlim([0 max(timeDec)])
        title('LFP')
        xlabel('Time [s]')
        ylabel('Voltage [mV]')
        hold all
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
        %% For each frequency band, approximate L_t
        k = zeros(size(LFPdataDecFilt,1),size(fBands,1));
        phi = zeros(size(LFPdataDecFilt,1),size(fBands,1));
        L_t = zeros(size(LFPdataDecFilt,1),size(fBands,1));
        for bb = 1:size(fBands,1)
            [k(:,bb),phi(:,bb)] = UP_DOWN_DET_compInstPhaseAmpl(LFPdataDecFilt(:,bb),phaseMethod);
            phi(phi(:,bb)<=0,bb) = phi(phi(:,bb)<=0,bb)+2*pi; % rescale from 0 to 2*pi
            phi(:,bb) = round(phi(:,bb)./(2*pi/360)).*(2*pi/360); % round to 1 degree resolution
            L_t(:,bb) = cos(phi(:,bb)-theta(bb));
        end
        %% Compute k_high
        [b_high,a_high] = ellip(2,0.1,40,fHigh./sfLFP);
        LFPdataHigh = filtfilt(b_high,a_high,LFPdata_dec);
        [k_high,phi_high] = UP_DOWN_DET_compInstPhaseAmpl(LFPdataHigh,'hilbert');
        %% Derive S_LFP
        K = k./repmat(k_high+sum(k,2),1,size(fBands,1));
        S_LFP = (sum(K.*L_t,2)+1)./2;
        %% Compute k_alphaBeta
        if combFlag==1
            [b_alphaBeta,a_alphaBeta] = ellip(2,0.1,40,fAlphaBeta./sfLFP);
            LFPdataAlphaBeta = filtfilt(b_alphaBeta,a_alphaBeta,LFPdata_dec);
            rmsLFPdata = rmsFilt(LFPdataAlphaBeta,win_samples);
            smoothLFPdata = smooth(rmsLFPdata,winSmooth_samples);
            S_alphaBeta = smoothLFPdata./max(smoothLFPdata);
            % Compute S_comb
            S_comb = (S_LFP+S_alphaBeta)./2;
        else
            S_comb = S_LFP;
        end
        %% Fit distribution of S_comb by GM
        options = statset('Display','final','MaxIter',500);
        try
            GMmodel = gmdistribution.fit(S_comb,3,'Options',options);
        catch err
            keyboard;
            continue;
        end
        mu = GMmodel.mu;
        sigma = squeeze(GMmodel.Sigma);
        figure(2)
        [hAmpl,bins] = hist(S_comb,nBins);
        bar(bins,hAmpl./length(S_comb).*100);
        ylabel('# occurrences')
        xlabel('S_{LFP}')
        hold all
        yFit = pdf(GMmodel,bins');
        plot(bins,yFit,'r');
        %% Compute thresholds for UP & DOWN state detection
        [~,idxUP] = max(mu);
        [~,idxDOWN] = min(mu);
        thUP = mu(idxUP)-nSTD_UP*sqrt(sigma(idxUP));
        thDOWN = mu(idxDOWN)+nSTD_DOWN*sqrt(sigma(idxDOWN));
        figure(2)
        line([thUP thUP],[0 max(yFit)],'LineStyle','--','LineWidth',1,'Color','c')
        line([thDOWN thDOWN],[0 max(yFit)],'LineStyle','--','LineWidth',1,'Color','g')
        %% Detect UP and DOWN states from LFP trace and save results
        [UP_states_DET, DOWN_states_DET] = UP_DOWN_DET_detectStates(S_comb, thUP, thDOWN, stateIntervTh, stateDurTh);
        UP_state_freq = length(UP_states_DET)./(length(LFPdata_dec)/sfLFP);
        DOWN_state_freq = length(DOWN_states_DET)./(length(LFPdata_dec)/sfLFP);
        mUPstateDur = mean((UP_states_DET(:,2)-UP_states_DET(:,1))./sfLFP);
        stdUPstateDur = std((UP_states_DET(:,2)-UP_states_DET(:,1))./sfLFP);
        mDOWNstateDur = mean((DOWN_states_DET(:,2)-DOWN_states_DET(:,1))./sfLFP);
        stdDOWNstateDur = std((DOWN_states_DET(:,2)-DOWN_states_DET(:,1))./sfLFP);
        UPDOWNstates_stat = [UP_state_freq DOWN_state_freq mUPstateDur stdUPstateDur mDOWNstateDur stdDOWNstateDur];
        figure(1)
        subplot(3,1,2)
        for uu = 1:length(UP_states_DET)
            plot(timeDec(UP_states_DET(uu,1):UP_states_DET(uu,2)),LFPdata_dec(UP_states_DET(uu,1):UP_states_DET(uu,2)),'Color',pinkColor,'LineWidth',2)
        end
        for dd = 1:length(DOWN_states_DET)
            plot(timeDec(DOWN_states_DET(dd,1):DOWN_states_DET(dd,2)),LFPdata_dec(DOWN_states_DET(dd,1):DOWN_states_DET(dd,2)),'Color','b','LineWidth',2)
        end
        % Compare with intracellularly detected states
        intraResults = load(fullfile(pathName,[expName, '_Trace_1_',seriesText,'_INTRA_DET_results.mat']));
        figure(1)
        subplot(3,1,1)
        for uu = 1:length(intraResults.UP_states_DET)
            plot(time(intraResults.boundaryWin+(intraResults.UP_states_DET(uu,1):intraResults.UP_states_DET(uu,2))),intraSeries_mV(intraResults.boundaryWin+(intraResults.UP_states_DET(uu,1):intraResults.UP_states_DET(uu,2))),'Color',pinkColor,'LineWidth',2)
        end
        for dd = 1:length(intraResults.DOWN_states_DET)
            plot(time(intraResults.boundaryWin+(intraResults.DOWN_states_DET(dd,1):intraResults.DOWN_states_DET(dd,2))),intraSeries_mV(intraResults.boundaryWin+(intraResults.DOWN_states_DET(dd,1):intraResults.DOWN_states_DET(dd,2))),'Color','b','LineWidth',2)
        end
        %% Plot S_comb & thresholds
        figure(1)
        subplot(3,1,3)
        plot(timeDec,S_comb)
        hold all
        line(xlim,[thUP thUP],'Color',pinkColor,'LineWidth',2,'LineStyle','--');
        line(xlim,[thDOWN thDOWN],'Color','b','LineWidth',2,'LineStyle','-.');
        xlim([0 max(timeDec)])
        %% LFP voltage relative to state start/end
        UP_states_INTRA_DET_dec = round((intraResults.UP_states_DET+intraResults.boundaryWin)./(sf/sfLFP));
        DOWN_states_INTRA_DET_dec = round((intraResults.DOWN_states_DET+intraResults.boundaryWin)./(sf/sfLFP));
        UP_states_signal_INTRA_DET = convert2stateSignal(UP_states_INTRA_DET_dec,length(LFPdata_dec));
        DOWN_states_signal_INTRA_DET = convert2stateSignal(DOWN_states_INTRA_DET_dec,length(LFPdata_dec));
        UP_states_DET(UP_states_DET(:,1)<=intraResults.boundaryWin/(sf/sfLFP) | UP_states_DET(:,2)<=intraResults.boundaryWin/(sf/sfLFP),:)=[];
        DOWN_states_DET(DOWN_states_DET(:,1)<=intraResults.boundaryWin/(sf/sfLFP) | DOWN_states_DET(:,2)<=intraResults.boundaryWin/(sf/sfLFP),:)=[];
        UP_states_signal_LFP_DET = convert2stateSignal(UP_states_DET,length(LFPdata_dec));
        DOWN_states_signal_LFP_DET = convert2stateSignal(DOWN_states_DET,length(LFPdata_dec));
        % UP_states_signal_INTRA_DET, DOWN_states_signal_INTRA_DET --> MP
        % detection
        % UP_states_signal_LFP_DET, DOWN_states_signal_LFP_DET --> LFP
        % detection
        TP_UP_signal = UP_states_signal_LFP_DET&UP_states_signal_INTRA_DET;
        FP_UP_negDOWN_signal = UP_states_signal_LFP_DET&DOWN_states_signal_INTRA_DET;   %LFP - UP
        FN_UP_negDOWN_signal = UP_states_signal_INTRA_DET&DOWN_states_signal_LFP_DET;   %INTRA - UP
        TP_DOWN_signal = DOWN_states_signal_LFP_DET&DOWN_states_signal_INTRA_DET;
        FP_DOWN_negUP_signal = DOWN_states_signal_LFP_DET&UP_states_signal_INTRA_DET;   %LFP - DOWN
        FN_DOWN_negUP_signal = DOWN_states_signal_INTRA_DET&UP_states_signal_LFP_DET;   %INTRA - DOWN
        states = {'UP','DOWN'};
        detection = {'LFP','INTRA'};
        bin_relative_hist = 0.1;
        hist_bins = -0.05:bin_relative_hist:1.05;
        color = {'r','k'};
        h_voltage_hist = figure('position',[50 50 1900 950]);
        h_FP_hist = figure('position',[50 50 1900 950]);
        for oo = 1:length(states)
            figure(h_voltage_hist)
            subplot(1,2,oo)
            hold all; linesHandle = zeros(2,1);
            figure(h_FP_hist)
            subplot(2,2,oo)
            hold all;  
            subplot(2,2,oo+2)
            hold all; barHandle = zeros(4,1);   
            for kk = 1:length(detection)
                signal = [states{oo},'_states_signal_',detection{kk},'_DET'];
                selected_bins = find(eval(signal)~=0);             
                state_voltage = LFPdata_dec(selected_bins);
                selected_bins_state_limits = cell(numel(selected_bins),1);
                selected_bins_relative_timing = zeros(numel(selected_bins),1);
                switch detection{kk}
                    case 'LFP'
                        up_extremes = UP_states_DET;
                        down_extremes = DOWN_states_DET;
                        % false positives
                        switch states{oo}
                            case 'UP'
                                state_false = FP_UP_negDOWN_signal(selected_bins);
                            case 'DOWN'
                                state_false = FP_DOWN_negUP_signal(selected_bins);
                        end                      
                    case 'INTRA'
                        up_extremes = UP_states_INTRA_DET_dec;
                        down_extremes = DOWN_states_INTRA_DET_dec;
                        % false negatives
                        switch states{oo}
                            case 'UP'
                                state_false = FN_UP_negDOWN_signal(selected_bins);
                            case 'DOWN'
                                state_false = FN_DOWN_negUP_signal(selected_bins);
                        end
                end
                points_counter = length(LFPdata_dec);
                for ind_bin=1:numel(selected_bins)
                    previous_events = [up_extremes(:);down_extremes(:)];
                    previous_events((previous_events-selected_bins(ind_bin))>=0) = [];
                    post_events = [up_extremes(:);down_extremes(:)];
                    post_events((post_events-selected_bins(ind_bin))<=0) = [];
                    selected_bins_state_limits{ind_bin} = [max([previous_events;1]) min([post_events;points_counter])]-selected_bins(ind_bin);
                    if any(selected_bins_state_limits{ind_bin}==0)
                        keyboard
                    end
                    selected_bins_relative_timing(ind_bin) = (selected_bins(ind_bin)-max([previous_events;1]))/(min([post_events;points_counter])-max([previous_events;1]));
                end
                relative_voltage = zeros(2,numel(hist_bins));
                relative_false = zeros(1,numel(hist_bins));
                for ind_hist=1:numel(hist_bins)
                    relative_voltage(1,ind_hist) = mean(state_voltage(hist_bins(ind_hist)-bin_relative_hist/2<selected_bins_relative_timing & hist_bins(ind_hist)+bin_relative_hist/2>selected_bins_relative_timing));
                    relative_voltage(2,ind_hist) = std(state_voltage(hist_bins(ind_hist)-bin_relative_hist/2<selected_bins_relative_timing & hist_bins(ind_hist)+bin_relative_hist/2>selected_bins_relative_timing));
                    relative_false(1,ind_hist) = sum(state_false(hist_bins(ind_hist)-bin_relative_hist/2<selected_bins_relative_timing & hist_bins(ind_hist)+bin_relative_hist/2>selected_bins_relative_timing));
                end 
                varName = [states{oo},'_states_',detection{kk},'_DET_relativeVoltage'];
                assignin('caller',varName,relative_voltage);
                varName2 = [states{oo},'_states_',detection{kk},'_DET_relativeFalse'];
                assignin('caller',varName2,relative_false);
                figure(h_voltage_hist)
                subplot(1,2,oo)
                hObject = shadedErrorBar(hist_bins,relative_voltage(1,:),relative_voltage(2,:),['-',color{kk}],1);
                linesHandle(kk) = hObject.mainLine;
                figure(h_FP_hist)
                subplot(2,2,sub2ind([2,2],oo,kk))
                barHandle(sub2ind([2,2],oo,kk)) = bar(hist_bins,relative_false,color{kk});
            end
            figure(h_voltage_hist)
            subplot(1,2,oo)
            ylabel(['Average LFP voltage during ',states{oo}])
            xlabel('Time relative to state extremes')
            legend(linesHandle,'LFP','MP')
            figure(h_FP_hist)
            subplot(2,2,oo)
            ylabel(['False positive ',states{oo}])
            xlabel('Time relative to LFP-detected state extremes')
            subplot(2,2,oo+2)
            ylabel(['False negative ',states{oo}])
            xlabel('Time relative to MP-detected state extremes')
        end
        saveas(h_voltage_hist,fullfile(saveFolderName,[expName, '_Trace_1_',seriesText,'_LFP_relativeVoltage.fig']),'fig');
        close(h_voltage_hist)
        saveas(h_FP_hist,fullfile(saveFolderName,[expName, '_Trace_1_',seriesText,'_LFP_relativeFP.fig']),'fig');
        close(h_FP_hist)
        save(fullfile(saveFolderName,[expName, '_Trace_1_',seriesText,'_LFP_relativeVoltage.mat']),'UP_states_LFP_DET_relativeVoltage',...
            'UP_states_INTRA_DET_relativeVoltage','DOWN_states_LFP_DET_relativeVoltage','DOWN_states_INTRA_DET_relativeVoltage')
        save(fullfile(saveFolderName,[expName, '_Trace_1_',seriesText,'_LFP_relativeFalse.mat']),'UP_states_LFP_DET_relativeFalse',...
            'UP_states_INTRA_DET_relativeFalse','DOWN_states_LFP_DET_relativeFalse','DOWN_states_INTRA_DET_relativeFalse')
        %% Compute performance: CoIn
        UP_state_freq_INTRA_DET = length(UP_states_INTRA_DET_dec)./(length(LFPdata_dec)/sfLFP);
        DOWN_state_freq_INTRA_DET = length(DOWN_states_INTRA_DET_dec)./(length(LFPdata_dec)/sfLFP);
        mUPstateDur_INTRA_DET = mean((UP_states_INTRA_DET_dec(:,2)-UP_states_INTRA_DET_dec(:,1))./sfLFP);
        stdUPstateDur_INTRA_DET = std((UP_states_INTRA_DET_dec(:,2)-UP_states_INTRA_DET_dec(:,1))./sfLFP);
        mDOWNstateDur_INTRA_DET = mean((DOWN_states_INTRA_DET_dec(:,2)-DOWN_states_INTRA_DET_dec(:,1))./sfLFP);
        stdDOWNstateDur_INTRA_DET = std((DOWN_states_INTRA_DET_dec(:,2)-DOWN_states_INTRA_DET_dec(:,1))./sfLFP);
        UPDOWNstates_stat_INTRA_DET = [UP_state_freq_INTRA_DET DOWN_state_freq_INTRA_DET...
            mUPstateDur_INTRA_DET stdUPstateDur_INTRA_DET...
            mDOWNstateDur_INTRA_DET stdDOWNstateDur_INTRA_DET];
        % CoIn
        CoIn_UP = UP_DOWN_DET_computeCoIn(UP_states_signal_INTRA_DET,UP_states_signal_LFP_DET);
        CoIn_DOWN = UP_DOWN_DET_computeCoIn(DOWN_states_signal_INTRA_DET,DOWN_states_signal_LFP_DET);
        CoIn = (CoIn_UP+CoIn_DOWN)/2;
        % Compute TPR and FPR
        TPR_UP_chosenTh = sum(UP_states_signal_LFP_DET&UP_states_signal_INTRA_DET)./sum(UP_states_signal_INTRA_DET);
        FPR_UP_chosenTh_negNotUP = sum(UP_states_signal_LFP_DET&(~UP_states_signal_INTRA_DET))./sum(~UP_states_signal_INTRA_DET);
        TPR_DOWN_chosenTh = sum(DOWN_states_signal_LFP_DET&DOWN_states_signal_INTRA_DET)./sum(DOWN_states_signal_INTRA_DET);
        FPR_DOWN_chosenTh_negNotDOWN = sum(DOWN_states_signal_LFP_DET&(~DOWN_states_signal_INTRA_DET))./sum(~DOWN_states_signal_INTRA_DET);
        %
        percFP_DOWN_chosenTh_negNotDOWN = sum(DOWN_states_signal_LFP_DET&(~DOWN_states_signal_INTRA_DET))./sum(DOWN_states_signal_LFP_DET);
        percFP_DOWN_chosenTh_negUP = sum(DOWN_states_signal_LFP_DET&(UP_states_signal_INTRA_DET))./sum(DOWN_states_signal_LFP_DET);
        %
        FPR_UP_chosenTh_negDOWN = sum(UP_states_signal_LFP_DET&DOWN_states_signal_INTRA_DET)./sum(DOWN_states_signal_INTRA_DET);
        FPR_DOWN_chosenTh_negUP = sum(DOWN_states_signal_LFP_DET&UP_states_signal_INTRA_DET)./sum(UP_states_signal_INTRA_DET);
        %
        figure(7)
        subplot(2,1,1)
        plot(timeDec,UP_states_signal_LFP_DET,'b',timeDec,FP_UP_negDOWN_signal,'r');
        %
        subplot(2,1,2)
        plot(timeDec,DOWN_states_signal_LFP_DET,'b',timeDec,FP_DOWN_negUP_signal,'r');
        %
        figure(8)
        subplot(2,1,1)
        [nelements,centers] = hist(LFPdata_dec(TP_UP_signal),50);
        plot(centers,nelements,'b.-')
        hold all
        [nelements,centers] = hist(LFPdata_dec(FP_UP_negDOWN_signal),50);
        plot(centers,nelements,'r.-')
        xlim([-1 1])
        legend('TP_{UP}','FP_{UP} - negDOWN')
        xlabel('Voltage amplitude [mV]')
        ylabel('Frequency')
        subplot(2,1,2)
        [nelements,centers] = hist(LFPdata_dec(TP_DOWN_signal),50);
        plot(centers,nelements,'b.-')
        hold all
        [nelements,centers] = hist(LFPdata_dec(FP_DOWN_negUP_signal),50);
        plot(centers,nelements,'r.-')
        xlim([-1 1])
        legend('TP_{DOWN}','FP_{DOWN} - negUP')
        xlabel('Voltage amplitude [mV]')
        ylabel('Frequency')
        %
        figure(9)
        fLow = [0.5 4];
        [b_lowFreqFilt,a_lowFreqFilt] = ellip(2,0.1,40,fLow./sfLFP);
        LFPdataLow = filtfilt(b_lowFreqFilt,a_lowFreqFilt,LFPdata_dec);
        [k_low,phi_low] = UP_DOWN_DET_compInstPhaseAmpl(LFPdataLow,'hilbert');
        subplot(2,1,1)
        [nelements,centers] = hist(k_low(TP_UP_signal),50);
        plot(centers,nelements,'b.-')
        hold all
        [nelements,centers] = hist(k_low(FP_UP_negDOWN_signal),50);
        plot(centers,nelements,'r.-')
        xlim([0 1])
        legend('TP_{UP}','FP_{UP} - negDOWN')
        xlabel('Instantaneous low power content')
        ylabel('Frequency')
        subplot(2,1,2)
        [nelements,centers] = hist(k_low(TP_DOWN_signal),50);
        plot(centers,nelements,'b.-')
        hold all
        [nelements,centers] = hist(k_low(FP_DOWN_negUP_signal),50);
        plot(centers,nelements,'r.-')
        xlim([0 1])
        legend('TP_{DOWN}','FP_{DOWN} - negUP')
        xlabel('Instantaneous low power content')
        ylabel('Frequency')
        %
        figure(10)
        [b_alphaBetaFreqFilt,a_alphaBetaFreqFilt] = ellip(2,0.1,40,fAlphaBeta./sfLFP);
        LFPdataAlphaBeta = filtfilt(b_alphaBetaFreqFilt,a_alphaBetaFreqFilt,LFPdata_dec);
        [k_alphaBeta,phi_alphaBeta] = UP_DOWN_DET_compInstPhaseAmpl(LFPdataAlphaBeta,'hilbert');
        subplot(2,1,1)
        [nelements,centers] = hist(k_alphaBeta(TP_UP_signal),50);
        plot(centers,nelements,'b.-')
        hold all
        [nelements,centers] = hist(k_alphaBeta(FP_UP_negDOWN_signal),50);
        plot(centers,nelements,'r.-')
        xlim([0 1])
        legend('TP_{UP}','FP_{UP} - negDOWN')
        xlabel('Instantaneous alpha/beta power content')
        ylabel('Frequency')
        subplot(2,1,2)
        [nelements,centers] = hist(k_alphaBeta(TP_DOWN_signal),50);
        plot(centers,nelements,'b.-')
        hold all
        [nelements,centers] = hist(k_alphaBeta(FP_DOWN_negUP_signal),50);
        plot(centers,nelements,'r.-')
        xlim([0 1])
        legend('TP_{DOWN}','FP_{DOWN} - negUP')
        xlabel('Instantaneous alpha/beta power content')
        ylabel('Frequency')
        %% Compute performance: ROC curve 1
        thArray = 0:0.01:1;
        TPR_UP = zeros(length(thArray),1);
        FPR_UP_negNotUP = zeros(length(thArray),1);
        FPR_UP_negDOWN = zeros(length(thArray),1);
        TPR_DOWN = zeros(length(thArray),1);
        FPR_DOWN_negNotDOWN = zeros(length(thArray),1);
        FPR_DOWN_negUP = zeros(length(thArray),1);
        for tt = 1:length(thArray)
            [cur_UP_states_DET, cur_DOWN_states_DET] = UP_DOWN_DET_detectStates_wBound(S_comb, thArray(tt), thArray(tt), 0, 0);
            %         [cur_UP_states_DET, cur_DOWN_states_DET] = UP_DOWN_DET_detectStates(S_comb, thArray(tt), thArray(tt), stateIntervTh, 0);
            cur_UP_states_signal = convert2stateSignal(cur_UP_states_DET,length(LFPdata_dec));
            cur_DOWN_states_signal = convert2stateSignal(cur_DOWN_states_DET,length(LFPdata_dec));
            TPR_UP(tt) = sum(cur_UP_states_signal & UP_states_signal_INTRA_DET)./sum(UP_states_signal_INTRA_DET);
            FPR_UP_negNotUP(tt) = sum(cur_UP_states_signal & (~UP_states_signal_INTRA_DET))./sum(~UP_states_signal_INTRA_DET);
            FPR_UP_negDOWN(tt) = sum(cur_UP_states_signal & DOWN_states_signal_INTRA_DET)./sum(DOWN_states_signal_INTRA_DET);
            TPR_DOWN(tt) = sum(cur_DOWN_states_signal & DOWN_states_signal_INTRA_DET)./sum(DOWN_states_signal_INTRA_DET);
            FPR_DOWN_negNotDOWN(tt) = sum(cur_DOWN_states_signal & (~DOWN_states_signal_INTRA_DET))./sum(~DOWN_states_signal_INTRA_DET);
            FPR_DOWN_negUP(tt) = sum(cur_DOWN_states_signal & UP_states_signal_INTRA_DET)./sum(UP_states_signal_INTRA_DET);
        end
        figure(5)
        subplot(1,2,1)
        [FPR_UP_sort,FPR_UP_sort_idx]=sort(FPR_UP_negNotUP,'ascend');
        plot(FPR_UP_sort,TPR_UP(FPR_UP_sort_idx),'r*-');
        AUC_UP_NegNotUP = sum([FPR_UP_sort(1); diff(FPR_UP_sort)].*TPR_UP(FPR_UP_sort_idx));
        hold all
        [FPR_DOWN_sort,FPR_DOWN_sort_idx]=sort(FPR_DOWN_negNotDOWN,'ascend');
        plot(FPR_DOWN_sort,TPR_DOWN(FPR_DOWN_sort_idx),'b*-');
        AUC_DOWN_NegNotDOWN = sum([FPR_DOWN_sort(1); diff(FPR_DOWN_sort)].*TPR_DOWN(FPR_DOWN_sort_idx));
        legend('UP state ROC','DOWN state ROC')
        title('Negative state: NotUP/NotDOWN')
        subplot(1,2,2)
        [FPR_UP_sort,FPR_UP_sort_idx]=sort(FPR_UP_negDOWN,'ascend');
        plot(FPR_UP_sort,TPR_UP(FPR_UP_sort_idx),'r*-');
        AUC_UP_NegDOWN = sum([FPR_UP_sort(1); diff(FPR_UP_sort)].*TPR_UP(FPR_UP_sort_idx));
        hold all
        [FPR_DOWN_sort,FPR_DOWN_sort_idx]=sort(FPR_DOWN_negUP,'ascend');
        plot(FPR_DOWN_sort,TPR_DOWN(FPR_DOWN_sort_idx),'b*-');
        AUC_DOWN_NegUP = sum([FPR_DOWN_sort(1); diff(FPR_DOWN_sort)].*TPR_DOWN(FPR_DOWN_sort_idx));
        legend('UP state ROC','DOWN state ROC')
        title('Negative state: DOWN/UP')
        %% Compute performance: ROC curve 2
        [X_ROC_UP,Y_ROC_UP,T_ROC_UP,AUC_UP,optOprPt_UP] = perfcurve(UP_states_signal_INTRA_DET,S_comb,1);
        [X_ROC_DOWN,Y_ROC_DOWN,T_ROC_DOWN,AUC_DOWN,optOprPt_DOWN] = perfcurve(DOWN_states_signal_INTRA_DET,-S_comb,1);
        optThUP = T_ROC_UP(X_ROC_UP==optOprPt_UP(1)&Y_ROC_UP==optOprPt_UP(2));
        optThDOWN = abs(T_ROC_DOWN(X_ROC_DOWN==optOprPt_DOWN(1)&Y_ROC_DOWN==optOprPt_DOWN(2)));
        figure(6)
        plot(X_ROC_UP,Y_ROC_UP,'r.-');
        hold all
        plot(X_ROC_DOWN,Y_ROC_DOWN,'b.-');
        legend('UP state ROC','DOWN state ROC')
        %%
        [~,phi_deltaFiltLFP] = UP_DOWN_DET_compInstPhaseAmpl(deltaFiltLFP,phaseMethod);
        phi_deltaFiltLFP(phi_deltaFiltLFP<=0) = phi_deltaFiltLFP(phi_deltaFiltLFP<=0)+2*pi; % rescale from 0 to 2*pi
        UPstart = UP_states_DET(:,1);
        UPend = UP_states_DET(:,2);
        DOWNstart = DOWN_states_DET(:,1);
        DOWNend = DOWN_states_DET(:,2);
        phaseValues_UP = phi_deltaFiltLFP(UP_states_signal_INTRA_DET~=0);
        phaseValues_DOWN = phi_deltaFiltLFP(DOWN_states_signal_INTRA_DET~=0);
        phaseValues_UPstart = phi_deltaFiltLFP(UPstart);
        phaseValues_DOWNstart = phi_deltaFiltLFP(DOWNstart);
        phaseValues_UPend = phi_deltaFiltLFP(UPend);
        phaseValues_DOWNend = phi_deltaFiltLFP(DOWNend);
        %% Plot phase values' distribution during UP & DOWN states
        phaseValues_UP = phaseValues_UP./(2*pi/360);
        phaseValues_DOWN = phaseValues_DOWN./(2*pi/360);
        phaseValues_UPstart = phaseValues_UPstart./(2*pi/360);
        phaseValues_DOWNstart = phaseValues_DOWNstart./(2*pi/360);
        phaseValues_UPend = phaseValues_UPend./(2*pi/360);
        phaseValues_DOWNend = phaseValues_DOWNend./(2*pi/360);
        hHist = figure;
        subplot(2,2,1)
        [yHist,xHist] = hist(phaseValues_UP,50);
        plot(xHist, yHist, 'Color', pinkColor, 'LineWidth', 2)
        xlim([0 360])
        xlabel('Phase [degrees]')
        ylabel('Number of elements')
        hold all
        [yHist,xHist] = hist(phaseValues_DOWN,50);
        plot(xHist, yHist, 'Color', 'b', 'LineWidth', 2)
        xlim([0 360])
        xlabel('Phase [degrees]')
        ylabel('Number of elements')
        subplot(2,2,2)
        [yHist,xHist] = hist(phaseValues_UPstart,50);
        plot(xHist, yHist, 'Color', pinkColor, 'LineWidth', 2)
        hold all
        [yHist,xHist] = hist(phaseValues_DOWNstart,50);
        plot(xHist, yHist, 'Color', 'b', 'LineWidth', 2)
        subplot(2,2,3)
        [yHist,xHist] = hist(phaseValues_UPend,50);
        plot(xHist, yHist, 'Color', pinkColor, 'LineWidth', 2)
        hold all   
        [yHist,xHist] = hist(phaseValues_DOWNend,50);
        plot(xHist, yHist, 'Color', 'b', 'LineWidth', 2)
        saveas(hHist,fullfile(saveFolderName,[expName, '_Trace_1_',seriesText,'_phaseValueDistrib_LFP.fig']),'fig');
        close(hHist)
        %% Save results
        save(fullfile(saveFolderName,[expName, '_Trace_1_',seriesText,'_parameters.mat']),'parameters_struct');
        saveas(1,fullfile(saveFolderName,[expName, '_Trace_1_',seriesText,'_LFP_DET_traceComparison.fig']),'fig');
        close(1)
        saveas(2,fullfile(saveFolderName,[expName, '_Trace_1_',seriesText,'_LFP_DET_S_comb_distrib.fig']),'fig');
        close(2)
        saveas(5,fullfile(saveFolderName,[expName, '_Trace_1_',seriesText,'_LFP_DET_ROC1.fig']),'fig');
        close(5)
        saveas(6,fullfile(saveFolderName,[expName, '_Trace_1_',seriesText,'_LFP_DET_ROC2.fig']),'fig');
        close(6)
        saveas(7,fullfile(saveFolderName,[expName, '_Trace_1_',seriesText,'_LFP_DET_TP_FP_signals.fig']),'fig');
        close(7)
        saveas(8,fullfile(saveFolderName,[expName, '_Trace_1_',seriesText,'_LFP_DET_TP_FP_amplDistrib.fig']),'fig');
        close(8)
        saveas(9,fullfile(saveFolderName,[expName, '_Trace_1_',seriesText,'_LFP_DET_TP_FP_instLowFreqPowerDistrib.fig']),'fig');
        close(9)
        saveas(10,fullfile(saveFolderName,[expName, '_Trace_1_',seriesText,'_LFP_DET_TP_FP_instAlphaBetaFreqPowerDistrib.fig']),'fig');
        close(10)
        save(fullfile(saveFolderName,[expName,'_Trace_1_',seriesText,'_UPDOWNstates_stat_LFP_DET.txt']),'UPDOWNstates_stat','-ASCII');
        save(fullfile(saveFolderName,[expName,'_Trace_1_',seriesText,'_UPDOWNstates_stat_INTRA_DET.txt']),'UPDOWNstates_stat_INTRA_DET','-ASCII');
        save(fullfile(saveFolderName,[expName,'_Trace_1_',seriesText,'_LFPdet_results.mat']),'UP_states_DET','DOWN_states_DET');
        save(fullfile(saveFolderName,[expName,'_Trace_1_',seriesText,'_phaseValueDistribData_LFP.mat']),...
            'phaseValues_UP','phaseValues_DOWN','phaseValues_UPstart','phaseValues_DOWNstart','phaseValues_UPend','phaseValues_DOWNend');
        CoIn_2save = [CoIn CoIn_UP CoIn_DOWN];
        save(fullfile(CoIn_path,[expName, '_Trace_1_',seriesText,'_CoIn.txt']),'CoIn_2save','-ASCII');
        perf_chosenTh = [TPR_UP_chosenTh FPR_UP_chosenTh_negNotUP FPR_UP_chosenTh_negDOWN ...
            TPR_DOWN_chosenTh FPR_DOWN_chosenTh_negNotDOWN FPR_DOWN_chosenTh_negUP percFP_DOWN_chosenTh_negNotDOWN percFP_DOWN_chosenTh_negUP];
        save(fullfile(perfChosenTh_path,[expName, '_Trace_1_',seriesText,'_performanceChosenTh.txt']),'perf_chosenTh','-ASCII');
        ROC_values_1 = [AUC_UP_NegNotUP AUC_UP_NegDOWN AUC_DOWN_NegNotDOWN AUC_DOWN_NegUP];
        save(fullfile(ROC_analysis_path,[expName, '_Trace_1_',seriesText,'_ROC1_analysis.txt']),'ROC_values_1','-ASCII');
        ROC_values_2 = [AUC_UP AUC_DOWN optThUP optThDOWN thUP thDOWN];
        save(fullfile(ROC_analysis_path,[expName, '_Trace_1_',seriesText,'_ROC2_analysis.txt']),'ROC_values_2','-ASCII');
    end
end