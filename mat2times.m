function spike_times = mat2times(mat, fs, unit, method)

% Description: Helper function converting binary spike matrix to spike times
%----------
% INPUT
% mat         - [num_nodes x num_samples] binary matrix (1:= event)
% fs          - [scalar] sampling frequency in Hz
% unit        - [string] unit of the spike times ('s', 'ms', 'frames')
% method      - [string] spike detection method
%----------
% OUTPUT
% spike_times - [cell] cell containing spike (event) times
%----------
% @author JJChabros (jjc80@cam.ac.uk), Feb 2021

num_nodes = min(size(mat));

if ~exist('unit', 'var')
    unit = 'frames';
end
if ~exist('fs', 'var')
    fs = 25000;
end
spike_times = cell([1,60]);
for node = 1:num_nodes
    switch unit
        case 's'
            spike_times{node}.(method) = find(mat(node,:)==1)/fs;
        case 'ms'
            spike_times{node}.(method) = find(mat(node,:)==1)/(fs/1000);
        case 'frames'
            spike_times{node}.(method) = find(mat(node,:)==1);
    end
end

end

