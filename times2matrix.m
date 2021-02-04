function spikeMatrix = times2matrix(spikeTimes, duration, method)
% Description: Helper function converting spike times cell to (binary) spike matrix
%----------
% INPUT
% spikeTimes  - [n x 1]  cell with spike time structures; spikeTimes{}.(method)
% duration    - [scalar] duration of the recordings in seconds
% method      - [string] spike detection method
%----------
% OUTPUT
% spikeMatrix - [n x t] binary matrix; 1 - spike, 0 - no spike
%----------
% @author JJChabros (jjc80@cam.ac.uk), January 2021

if ~exist('method', 'var')
    method = 'mea';
end

spikeMatrix = zeros(length(spikeTimes), duration);
for i = 1:length(spikeTimes)
    spikeMatrix(i, spikeTimes{i}.(method)) = 1;
end

end