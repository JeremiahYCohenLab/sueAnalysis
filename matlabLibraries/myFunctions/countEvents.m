function respLats = countEvents(spiketimes, alignTime, window)
    respLats = cell(length(alignTime), 1);
    for i = 1:length(alignTime)
        spikeInds = (spiketimes > alignTime(i) + window(1)) & (spiketimes <alignTime(i) + window(2));
        respLats{i} = spiketimes(spikeInds) - alignTime(i); 
    end
end