function downsampledMatrix = event_based_downsampling(spikeMatrix, dsFactor)

% Determine the number of channels (n) and the number of datapoints (t) in the input matrix
[n, t] = size(spikeMatrix);

% Initialize the output matrix
downsampledMatrix = zeros(n, ceil(t/dsFactor));

% Iterate over each channel
for i = 1:n
    % Find the indices of the non-zero elements (spike times) for this channel
    spikeIndices = find(spikeMatrix(i,:));
    
    % Determine the number of spikes in this channel
    numSpikes = length(spikeIndices);
    
    % If there are no spikes in this channel, move on to the next channel
    if numSpikes == 0
        continue;
    end
    
    % Initialize the downsampled spike train for this channel
    downsampledSpikeTrain = zeros(1, ceil(t/dsFactor));
    
    % Initialize the variable that tracks the current downsampled index
    downsampledIndex = 1;
    
    % Iterate over the indices of the spikes in this channel
    for j = 1:numSpikes
        % Determine the downsampled index for this spike time
        downsampledSpikeIndex = ceil(spikeIndices(j)/dsFactor);
        
        % If this downsampled index is the same as the previous one, skip this spike
        if downsampledSpikeIndex == downsampledIndex
            continue;
        end
        
        % Otherwise, mark the downsampled index as containing a spike
        downsampledSpikeTrain(downsampledSpikeIndex) = 1;
        
        % Update the current downsampled index
        downsampledIndex = downsampledSpikeIndex;
    end
    
    % Add the downsampled spike train for this channel to the output matrix
    downsampledMatrix(i,:) = downsampledSpikeTrain;
end
