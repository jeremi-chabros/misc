% Requires customHeatmap.m
% Download from: https://github.com/jeremi-chabros/misc/blob/main/customHeatmap.m

%% Load data

clearvars; clc;

spike_file_name = '200708_slice1_1_L_-0.25071_spikes (1).mat';
raw_file_name = spike_file_name(1:strfind(spike_file_name, '_L_')-1);

load(spike_file_name);
load(raw_file_name);

duration = spikeDetectionResult.params.duration;
fs = spikeDetectionResult.params.fs;

%% Filtering
lowpass = 600;
highpass = 8000;
wn = [lowpass highpass] / (fs / 2);
filterOrder = 3;
[b, a] = butter(filterOrder, wn);
filteredData = filtfilt(b, a, double(dat));

methods = sort(fieldnames(spikeTimes{1}));

bin_s = 10;

%% Plot moving average summed across all electrodes

for i = 1:length(methods)
    spk_vec_all = zeros(1, duration*fs);
    method = methods{i};
    for j = 1:60
        spk_times = round(spikeTimes{j}.(method)*fs);
        spike_count(j) = length(spk_times);
        spk_vec = zeros(1, duration*fs);
        spk_vec(spk_times) = 1;
        spk_vec_all = spk_vec_all+spk_vec;
    end
    spike_counts.(method) = spike_count;
    spike_freq.(method) = spike_count/duration;
    plot(movmean(spk_vec_all, bin_s*fs))
    hold on
end
xticks(linspace(1, duration*fs,10))
xticklabels(round(linspace(1, duration, 10)));
pbaspect([2, 1, 1])
box off
legend(strrep(methods,'p','.'))
xlabel('Time (s)')

%% Plot spike count heatmaps
% TODO: use just the merged method
for i = 1:length(methods)
    close all;
    method = methods{i};
    [f, c] = customHeatmap(spike_freq.(method), 'c_map', 'thermal',...
        'cbarTitle', 'Spiking frequency (Hz)');
end