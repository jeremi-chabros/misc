function [all_spikes, intersect_matrix, unique_idx] = mergeSpikes(spike_times, option)
% Description: Merges spike times across detection methods.
%----------
% INPUT
% spike_times - [n x 1 cell] containing for each electrode a structure
% option      - [string] specifies merging option: 'all' or 'wavelets'; if
%                        not specified, defaults to 'all'
%----------
% OUTPUT
% all_spikes  - [struct] with merged spike times appended to spike_times{channel}.('all')
% intersect_matrix - [table] containing for each method a binary vector specifying which
%                            spikes were detected by a given method
% unique_idx - [vector] containing the indexes of spikes detected uniquely by a given method

% TIP: to get spike waveforms corresponding to the merged spike times call e.g.:
%      [~, spikes] = alignPeaks(spike_times.all, filtered_trace, 25, 0, 1, 10, 10);
%----------
% @author JJChabros (jjc80@cam.ac.uk), January 2021

% Define a list of admissible wavelets
wavelets = {'mea','bior1.5','bior1.3','db2'};
if ~exist('option', 'var')
    option = 'all';
end
methods = fieldnames(spike_times);
all_spikes = [];

switch option
    case 'all'
        for method = 1:numel(methods)
            all_spikes = union(all_spikes, spike_times.(methods{method}));
        end
    case 'wavelets'
        for m = 1:length(wavelets)
            method = strrep(wavelets{m},'.','p');
            if isfield(spike_times, method)
                all_spikes = union(all_spikes, spike_times.(method));
            end
        end
end
all_spikes = all_spikes';

intersect_matrix = table;
for method = 1:length(methods)
    spk_vec = spike_times.(methods{method});
    intersect_vec = zeros(length(all_spikes),1);
    for spikeIndex = 1:length(all_spikes)
        if ismember(all_spikes(spikeIndex), spk_vec)
            intersect_vec(spikeIndex) = 1;
        end
    end
    intersect_matrix.(methods{method}) = intersect_vec;
end

unique_idx = zeros(1, height(intersect_matrix));
for spike = 1:height(intersect_matrix)
    ff = find(intersect_matrix{spike, :} == 1);
    if length(ff) == 1 && ff ~=0
        unique_idx(spike) = ff;
    end
end

end
