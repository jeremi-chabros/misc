function erank = effectiveRank(spikeTimes, spikeDetectionResult, binMs, method, repNum)

fs = spikeDetectionResult.params.fs;
duration_s = spikeDetectionResult.params.duration;
duration = duration_s*fs;
unit = spikeDetectionResult.params.unit;
spikeMatrix = times2matrix(spikeTimes, duration_s, fs, method, unit);
dsMatrix = downSample(spikeMatrix, '', binMs*fs/1000, 'sum', 1);

erankReal = efrank(dsMatrix);

for rep = 1:repNum
    synthMatrix = zeros(size(dsMatrix));
    
    % Random shuffling
    for i = 1:60
        synthMatrix(i,:) = circshift(dsMatrix(i,:), randi(length(dsMatrix),1));
    end

    erankSynth(rep) = efrank(synthMatrix);
end

erank = erankReal/mean(erankSynth);
end