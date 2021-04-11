clearvars; clc;

%% GET FILES
dataPath = '/Users/jeremi/mea/spikes';
addpath(dataPath);
files = dir([dataPath filesep '*.mat']);

%% PARAMS
params = struct;

params.binMs = 25;
params.lag_ms = 25;
params.method = 'merged';
params.repNumErank = 1;
params.minSpikingFrequency = 1; % in Hz
params.tail = 0.05;
params.repNumAdjm = 500;

%%

for file = 1:1
    load(files(file).name, 'spikeTimes', 'spikeDetectionResult');
    duration_s = spikeDetectionResult.params.duration;
    fs = spikeDetectionResult.params.fs;
    merged_spikes = [];
    adjM = [];
    adjMci = [];
    erank = [];
    spiking_frequency = zeros(60,1);
    
    for i = 1:60
        [merged_spikes, ~, ~] = mergeSpikes(spikeTimes{i});
        if length(merged_spikes)/duration_s > params.minSpikingFrequency
            spikeTimes{i}.merged = merged_spikes;
            spiking_frequency(i) = length(merged_spikes)/duration_s;
        else
            spikeTimes{i}.merged = [];
            spiking_frequency(i) = 0;
        end
    end
    
    erank(file) = effectiveRank(spikeTimes, spikeDetectionResult, params.binMs, params.method, params.repNumErank);
    
    [adjM, adjMci] = adjM_thr_parallel(spikeTimes, params.method, params.lag_ms, params.tail, fs,...
    duration_s, params.repNumAdjm);
    
    fname = files(file).name;
    saveName = fname(1:strfind(fname, '_L_')-1);
    vars2save = {'erank', 'adjM', 'adjMci', 'spiking_frequency',...
        'spikeDetectionResult', 'params'};
    save([saveName '_results.mat'], vars2save);
end