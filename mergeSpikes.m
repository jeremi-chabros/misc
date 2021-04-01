function [all_spikes, unique_idx, intersect_matrix, methods] = mergeSpikes(spike_times, option)
% Description: Merges spike times across detection methods.
%----------
% INPUT
% spike_times - [n x 1 cell] containing for each electrode a structure
% option      - [string] specifies merging option: 'all' or 'wavelets'; if
%                        not specified, defaults to 'all'
%----------
% OUTPUT
% all_spikes  - [struct] with merged spike times appended to spike_times{channel}.('all')
% unique_idx - [vector] containing the indexes of spikes detected uniquely by a given method

% TIP: to get spike waveforms corresponding to the merged spike times call e.g.:
%      [~, spikes] = alignPeaks(spike_times.all, filtered_trace, 25, 0, 1, 10, 10);
%----------
% @author JJChabros (jjc80@cam.ac.uk), January 2021

% NOTE: the number in unique_idx corresponds to the position in fieldnames
%       of spike_times, e.g. if field no. 1 is 'mea', then 1 in unique_idx
%       will represent a spike detected uniquely by 'mea' method

switch option
    case 'all'
        methods = fieldnames(spike_times);
    case 'wavelets'
        methods = {'mea','bior1p5','bior1p3','db2'};
end

if ~exist('option', 'var')
    option = 'all';
end

all_spikes = [];
for method = 1:numel(methods)
    all_spikes = union(all_spikes, spike_times.(methods{method}));
end
all_spikes = all_spikes';
% Vectorized operations + Logical indexing = Fast but hard to follow

% Pre-allocate
all_spikes_bin = zeros(max(all_spikes),1);
all_spikes_bin(all_spikes) = 1;
intersect_matrix = zeros(length(all_spikes_bin),length(methods));

for method = 1:length(methods)
    spk_vec = zeros(size(all_spikes_bin));
    intersect_vec = zeros(size(all_spikes_bin));
    spk_vec(spike_times.(methods{method})) = 1;
    intersect_vec = spk_vec .* all_spikes_bin;
    intersect_matrix(:,method) = intersect_vec';
end

idx = sum(intersect_matrix,2)>0;
J = intersect_matrix(idx,:);
unique_idx = sum(J,2);

intersect_matrix = intersect_matrix(any(intersect_matrix,2),:);
end

% TODO: make it work
% all_unique_idx = intersect_matrix(sum(intersect_matrix,2)==1,:);
% unique_counts = sum(all_unique_idx);
% f = zeros(length(all_unique_idx),1);
% [r,c] = find(all_unique_idx);
% f(r)=c;