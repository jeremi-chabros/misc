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
%%
% Vectorized operations + Logical indexing strike again...
% Fast but hard to follow
all_spikes_bin = zeros(1, max(all_spikes));
all_spikes_bin(all_spikes)=1;

intersect_matrix = table;
for method = 1:length(methods)
    spk_vec = zeros(size(all_spikes_bin));
    intersect_vec = zeros(size(all_spikes_bin));
    spk_vec(spike_times.(methods{method}))=1;
    intersect_vec = spk_vec.*all_spikes_bin;
    intersect_matrix.(methods{method}) = intersect_vec';
end

unique_idx = zeros(1, height(intersect_matrix));
rows = sum(intersect_matrix{:,:},2);
[~,unique_spks] = find(intersect_matrix{rows==1,:} == 1);
unique_idx(rows==1) = unique_spks;
end
