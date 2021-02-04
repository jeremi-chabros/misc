function [all_spikes, intersect_matrix, unique_idx] = mergeSpikes(spike_times)

methods = fieldnames(spike_times);
all_spikes = [];
for method = 1:numel(methods)
    all_spikes = union(all_spikes, spike_times.(methods{method}));
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