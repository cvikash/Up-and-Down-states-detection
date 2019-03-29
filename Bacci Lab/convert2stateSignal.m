function states_signal = convert2stateSignal(states_DET,len)
states_signal = zeros(len,1);
if ~isempty(states_DET)
    idx = cell2mat((cellfun(@(x,y) [str2double(x):str2double(y)], cellstr(num2str(states_DET(:,1))),cellstr(num2str(states_DET(:,2))),'UniformOutput',false))')';
    states_signal(idx) = 1;
end