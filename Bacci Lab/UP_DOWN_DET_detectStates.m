function [UP_states_DET, DOWN_states_DET] = UP_DOWN_DET_detectStates(variable, thUP, thDOWN, minInterv, minDuration)
% UP_DOWN_DET_detectStates.m takes the decision variable, the thresholds, and parameters of min. state interval and min. state duration as input, and 
% and returns the detected UP and DOWN states as output.
% License: V. Pasquale, 2018, Optical Approaches to Brain Function, CC BY 4.0
% ----Input arguments----:
% *     variable: decision variable (output of estimateUPDOWNdetTh.m)
% *     thUP: up state threshold (output of estimateUPDOWNdetTh.m)
% *     thDOWN: down state threshold (output of estimateUPDOWNdetTh.m)
% *     minInterv: minimum interval between states, in samples
% *     minDuration: minimum duration of a state, in samples
% ----Output arguments----:
% *     UP_states_DET: matrix numStates x 2, containing in the first column
%                   the start sample and in the second column the end sample of each state
% *     DOWN_states_DET: matrix numStates x 2, containing in the first column
%                   the start sample and in the second column the end sample of each state
if thDOWN >= thUP
    errordlg('DOWN state detection threshold should be lower than UP state detection threshold.');
    UP_states_DET = [];
    DOWN_states_DET = [];
    return
end
if minInterv >= minDuration
    errordlg('The minimum interval allowed between states of the same type should be lower than the minimum required duration of a state.');
    UP_states_DET = [];
    DOWN_states_DET = [];
    return
end    
%% UP state detection
UP_states = variable>thUP;
if all(UP_states)
    UP_states_DET = [1 length(variable)];
elseif any(UP_states)
    UP_states_start = find(diff([0; UP_states])==1);
    UP_states_end = find(diff([UP_states; 0])==-1);
    % exclude border states
    if UP_states_start(1)==1 % then discard the first one
        UP_states_start = UP_states_start(2:end);
        UP_states_end = UP_states_end(2:end);
    end
    if UP_states_end(end)==length(variable)
        UP_states_start = UP_states_start(1:end-1);
        UP_states_end = UP_states_end(1:end-1);
    end
    if ~isempty(UP_states_start) && ~isempty(UP_states_end)
        %% join periods of UP-state (above threshold) that are separated by less than a pre-defined time threshold
        UP_states_DET = [UP_states_start UP_states_end]; % samples
        UP_states_interval = UP_states_DET(2:end,1)-UP_states_DET(1:end-1,2);
        UP_states_2join = find(UP_states_interval<=minInterv); % Min interval minInterv
        joinPeriods = [UP_states_DET(UP_states_2join,2) UP_states_DET(UP_states_2join+1,1)];
        UP_states_signal = convert2stateSignal(UP_states_DET,length(variable));
        joinPeriodsSignal = convert2stateSignal(joinPeriods,length(variable));
        UP_states_signal_join = UP_states_signal | joinPeriodsSignal;
        UP_states_start = find(diff([0; UP_states_signal_join])==1);
        UP_states_end = find(diff([UP_states_signal_join; 0])==-1);
        UP_states_DET = [UP_states_start UP_states_end]; % samples
        %% delete UP states shorter than a pre-defined time threshold
        UP_state_dur = UP_states_DET(:,2)-UP_states_DET(:,1);
        UP_states_OK = UP_state_dur>= minDuration; %minimum duration minDuration
        UP_states_DET = UP_states_DET(UP_states_OK,:);
    else
        UP_states_DET = [];
    end
else
    UP_states_DET = [];
end
%% DOWN state detection
DOWN_states = variable<thDOWN;
if all(DOWN_states)
    DOWN_states_DET = [1 length(variable)];
elseif any(DOWN_states)
    DOWN_states_start = find(diff([0; DOWN_states])==1);
    DOWN_states_end = find(diff([DOWN_states; 0])==-1);
    % exclude borders
    if DOWN_states_start(1)==1 % then discard the first one
        DOWN_states_start = DOWN_states_start(2:end);
        DOWN_states_end = DOWN_states_end(2:end);
    end
    if DOWN_states_end(end)==length(variable)
        DOWN_states_start = DOWN_states_start(1:end-1);
        DOWN_states_end = DOWN_states_end(1:end-1);
    end
    if ~isempty(DOWN_states_start) && ~isempty(DOWN_states_end)
        %% join periods of DOWN-state (below threshold) that are separated by less than a pre-defined time threshold
        DOWN_states_DET = [DOWN_states_start DOWN_states_end]; % samples
        DOWN_states_interval = DOWN_states_DET(2:end,1)-DOWN_states_DET(1:end-1,2);
        DOWN_states_2join = find(DOWN_states_interval<=minInterv); %minimum interval minInterv
        joinPeriods = [DOWN_states_DET(DOWN_states_2join,2) DOWN_states_DET(DOWN_states_2join+1,1)];
        DOWN_states_signal = convert2stateSignal(DOWN_states_DET,length(variable));
        joinPeriodsSignal = convert2stateSignal(joinPeriods,length(variable));
        DOWN_states_signal_join = DOWN_states_signal | joinPeriodsSignal;
        DOWN_states_start = find(diff([0; DOWN_states_signal_join])==1);
        DOWN_states_end = find(diff([DOWN_states_signal_join; 0])==-1);
        DOWN_states_DET = [DOWN_states_start DOWN_states_end]; % samples
        %% delete DOWN states shorter than a pre-defined time threshold
        DOWN_state_dur = DOWN_states_DET(:,2)-DOWN_states_DET(:,1);
        DOWN_states_OK = DOWN_state_dur>=minDuration; %minDuration minDuration
        DOWN_states_DET = DOWN_states_DET(DOWN_states_OK,:);
    else
        DOWN_states_DET = [];
    end
else
    DOWN_states_DET = [];
end