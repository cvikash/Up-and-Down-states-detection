function rmsFiltData = rmsFilt(data,win)
nSamples = length(data);
rmsFiltData = zeros(nSamples,1);
for ii = 1:((floor(nSamples/win)-1)*win)
    rmsFiltData(ii) = rms(data((1:win)+(ii-1)));    
end