clearvars; clc;
dataPath = '/Users/jeremi/mea/data/PV-ArchT/';
addpath(dataPath);
addpath(genpath('/Users/jeremi/SpikeDetection-Toolbox'));
file = load('PAT200219_2C_DIV170002.mat');
%%

global fs duration TP FP FN nSamples
fs = 25000;
TP = 0;
FP = 0;
FN = 0;
nSamples = 0;
duration = length(file.dat)/fs;

% channel = randi(60,1);
channel = 12;
trace_raw = file.dat(1:fs*60, channel);
spikeTimes = struct;

Ns = 10;
multiplier = 2.5;
wnames = {'mea','bior1.5','bior1.3', 'db2'};
% wnames = {'mea'};
Wid = [0.5 0.8];
L = -0.25;
nSpikes = 200;
ttx = 0;
minPeakThrMultiplier = 2;
maxPeakThrMultiplier = 10;
posPeakThrMultiplier = 10;

% Run CWT spike detection
% tic
for wname = wnames
    good_wname = strrep(wname,'.','p');
    
    [spikeTimes.(good_wname{1}), ~, trace] = detectSpikesCWT(...
    trace_raw, fs, Wid, wname{1}, L, Ns, multiplier, nSpikes, ttx, ...
    minPeakThrMultiplier, maxPeakThrMultiplier, posPeakThrMultiplier);
%     length(spikeTimes.(good_wname{1}))
end
% toc
% %%
% threshold = mean(trace) - median(abs(trace-mean(trace)));
% trace_cropped = trace;
% plot3(zeros(1,length(trace_cropped)),1:length(trace_cropped),trace_cropped,...
%     'color', 'k',...
%     'linewidth', 1);
% 
% hold on
% 
% lineStyles=linspecer(4,'colorblind');
% spikepos = spikeTimes.mea;
% % spikepos = spikepos(intersect(spikepos>10000,spikepos<3000));
% scatter3( zeros(1,length(spikepos)),spikepos, repmat(4*std(trace)-(1*threshold/-7), length(spikepos), 1), 'v', 'filled',...
%     'markerfacecolor', lineStyles(1,:));
% 
% zlim([-20 20]);
% ylim([1000, 2000]);
% ylabel('time');
% zlabel('voltage');
% set(gca, 'ydir','reverse',...
%     'xcolor','none');
% grid on;

%%
% Run SWTTEO spike detection
params.filter = 1;
in.M = trace_raw;
in.SaRa = fs;
params.method = 'auto';
[spikepos1, ~] = SWTTEO(in,params);
[spikeTimes.('swtteo'), ~] = alignPeaks(spikepos1, trace, 10,...
    1, minPeakThrMultiplier, maxPeakThrMultiplier, posPeakThrMultiplier);
%%
% Threshold spike detection
[frames, ~, threshold] = detectSpikesThreshold(trace, 2.5, 0.1, fs, 0);
[spikeTimes.('threshold'), spike_waves_thr] = alignPeaks(find(frames==1), trace, 10,...
    1, minPeakThrMultiplier, maxPeakThrMultiplier, posPeakThrMultiplier);

% Merge spikes from all methods
unique_idx = [];
spikeTimes.all = [];
spikeTimes.wavelets = [];
[spikeTimes.('all'), unique_idx] = mergeSpikes(spikeTimes, 'all');

[~, spikeWaveforms] = alignPeaks(spikeTimes.('all'), trace, 10,...
    0, minPeakThrMultiplier, maxPeakThrMultiplier, posPeakThrMultiplier);

% Merge spikes from wavelets
if numel(wnames) > 1
    [spikeTimes.('wavelets'), ~] = mergeSpikes(spikeTimes, 'wavelets');
end

%% Plot all spike markers

global bin_ms
bin_ms = 100;

figure
set(gcf,'color','w');
h = plot(trace-mean(trace), 'k', 'linewidth', 1);
hold on

if numel(wnames) > 1
admissible = {'threshold', 'wavelets', 'swtteo', 'all'};
else
    admissible = {'threshold','swtteo','mea', 'all'};
end
lineStyles=linspecer(4,'colorblind');

for m = 1:length(admissible)
    spikepos = spikeTimes.(admissible{m});
    scatter(spikepos, repmat(4*std(trace)-(m*threshold/-7), length(spikepos), 1), 'v', 'filled',...
        'markerfacecolor', lineStyles(m,:));
end

hold on
yl = yline(threshold, 'magenta--', "Threshold = " + round(threshold,2) + "\muV");
yl.LabelVerticalAlignment = 'bottom';
yl.LineWidth = 1.5;
yl.FontSize = 10;

% yl = yline(minPeakThrMultiplier*(median(trace)-median(abs(trace-mean(trace))))/0.6745, 'cyan--', "Threshold = " + round(3/3.5*threshold) + "\muV");
% yl.LabelVerticalAlignment = 'top';
% yl.LineWidth = 1.5;
% yl.FontSize = 10;

% Aesthetics
st = 1;
en = st+(bin_ms*25);
xlim([st en]);
ylim([-6*std(trace) 5*std(trace)])
pbaspect([5 1 1])
box off
title("Electrode " + channel)
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend(admissible, 'location', 'northeastoutside')
set(gcf, 'position', [300 300 1500 350]);

%----------------------------------------
% Create buttons and callbacks
set(gcf,'KeyPressFcn',@keys);

nextButton = uicontrol('Position',[700 1 150 30],'String','Next',...
    'Callback', @Next);

prevButton = uicontrol('Position',[500 1 150 30],'String','Previous',...
    'Callback', @Prev);

saveButton = uicontrol('Position',[1325, 120, 100, 30],'String','Save',...
    'Callback', @Save);
%----------------------------------------
%% Plot unique spike markers

bin_ms = 100;
figure
set(gcf,'color','w');
h = plot(trace-mean(trace), 'k', 'linewidth', 1);
hold on
methods = fieldnames(spikeTimes);
lineStyles=linspecer(length(methods));
clear unique_counts
% unique_counts = table('rownames', methods);
unique_counts = table;

for m = 1:length(methods)
    if ~(strcmp(methods{m}, 'all') || strcmp(methods{m}, 'wavelets'))
        spks = spikeTimes.all;
        spikepos = find(unique_idx == m);
        
        scatter(spks(spikepos), repmat(5*std(trace)-(m*threshold/-8), length(spikepos), 1), 'v', 'filled',...
            'markerfacecolor', lineStyles(m,:));
        
        unique_counts.(methods{m}) = length(spikepos);
    elseif strcmp(methods{m}, 'all')
        spikepos = spikeTimes.(methods{m});
        scatter(spikepos, repmat(5*std(trace)-(m*threshold/-7), length(spikepos), 1), 'v', 'filled',...
            'markerfacecolor', lineStyles(m,:));
    end
end
% unique_counts.Properties.VariableNames = {'No. unique spikes'};
unique_counts(end-1:end-2,:) = [];
hold on
yl = yline(threshold, 'magenta--', "Threshold = " + round(threshold,2) + "\muV");
yl.LabelVerticalAlignment = 'bottom';
yl.LineWidth = 1.5;
yl.FontSize = 10;
clc
unique_counts

% yl = yline(3/3.5*threshold, 'cyan--', "Threshold = " + round(3/3.5*threshold) + "\muV");
% yl.LabelVerticalAlignment = 'top';
% yl.LineWidth = 1.5;
% yl.FontSize = 10;

% Aesthetics
st = 1;
en = st+(bin_ms*25);
xlim([st en]);
ylim([-6*std(trace) 6*std(trace)])
pbaspect([5 1 1])
box off
title("Electrode " + channel)
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend(methods, 'location', 'northeastoutside')
set(gcf, 'position', [300 300 1500 350]);


%----------------------------------------
% Create buttons and callbacks
set(gcf,'KeyPressFcn',@keys);

nextButton = uicontrol('Position',[700 1 150 30],'String','Next',...
    'Callback', @Next);

prevButton = uicontrol('Position',[500 1 150 30],'String','Previous',...
    'Callback', @Prev);

saveButton = uicontrol('Position',[1325, 120, 100, 30],'String','Save',...
    'Callback', @Save);

%% Plot unique spike waveforms
methods = fieldnames(spikeTimes);
t = tiledlayout(1, length(methods)-2);
t.Title.String = 'Spikes detected uniquely by each method';
for i = 1:length(methods)-2
    method = methods{i};
    if ~strcmp(method, 'all')
        
        spk_method = find(unique_idx == i);
        spk_waves_method = spikeWaveforms(spk_method, :);
        nexttile
        plot(spk_waves_method', 'linewidth', 0.1, 'color', [0.7 0.7 0.7])
        hold on
        plot(mean(spk_waves_method), 'linewidth', 1.5, 'color', [0 0 0])
        title({[method],["No. unique spikes: " + length(spk_method)]})
        box off;
        axis tight
        ylim([-6*std(trace) 5*std(trace)]);
        ylabel('Voltage [\muV]')
        set(gca, 'xcolor', 'none');
        yl = yline(threshold, 'r--');
        yl.LineWidth = 1.5;
    end
end
set(gcf, 'color', 'w');
%% UI functions
function keys(src, event)
global TP FP FN nSamples bin_ms
switch event.Key
    case 'rightarrow'
        Next;
    case 'leftarrow'
        Prev;
    case 'z'
        TP = TP + 1;
        nSamples = nSamples + 1;
    case 'x'
        FP = FP + 1;
        nSamples = nSamples + 1;
    case 'c'
        FN = FN + 1;
        nSamples = nSamples + 1;
    case 'r'
        FN = 0;
        TP = 0;
        FP = 0;
        nSamples = 0;
end
delete(findall(gcf,'type','annotation'));

annotation('textbox', [0.91, 0.5, 0.1, 0.1],...
    'String', "Sensitivity: " + round(TP/(TP + FN), 2)*100 + "%",...
    'FontSize', 14, 'edgecolor', 'none');

annotation('textbox', [0.91, 0.4, 0.1, 0.1],...
    'String', "Precision: " + round(TP/(TP + FP), 2)*100 + "%",...
    'FontSize', 14, 'edgecolor', 'none')
end

function Next(nextButton, EventData)
global duration fs bin_ms
lims = get(gca, 'xlim');
st = lims(2);
en = st+(bin_ms*25);
if en <= duration*fs
    set(gca,'xlim', [st en]);
end
end

function Prev(prevButton, EventData)
global bin_ms
lims = get(gca, 'xlim');
en = lims(1);
st = en-(bin_ms*25);
if st >= 1
    set(gca,'xlim', [st en]);
end
end

function Save(saveButton, EventData)
% save("MPT200220_3A_DIV21" + "stats.mat", 'TP', 'TN', 'FP', 'nSamples');
end
