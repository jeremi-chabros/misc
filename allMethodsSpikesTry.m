clearvars; clc;
dataPath = '/Users/jeremi/mea/data/PV-Ai32/';
addpath(dataPath);
addpath(genpath('/Users/jeremi/SpikeDetection-Toolbox'));
file = load('PVAi32190904_2E_DIV26.mat');
%%

global fs duration TP FP FN nSamples
fs = 25000;
TP = 0;
FP = 0;
FN = 0;
nSamples = 0;
duration = length(file.dat)/fs;

% channel = randi(60,1);
channel=15;
trace_raw = file.dat(1:fs*60, channel);
spikeTimes = struct;

Ns = 2;
multiplier = 3.5;
wnames = {'mea','bior1.5','bior1.3', 'db2'};
Wid = [0.5 1];
L = -0.43;
nSpikes = 200;
ttx = 0;
minPeakThrMultiplier = 3;
maxPeakThrMultiplier = 10;
posPeakThrMultiplier = 10;

% Run CWT spike detection
for wname = wnames
    good_wname = strrep(wname,'.','p');
    
    [spikeTimes.(good_wname{1}), ~, trace] = detectSpikesCWT(...
    trace_raw, fs, Wid, wname{1}, L, Ns, multiplier, nSpikes, ttx, ...
    minPeakThrMultiplier, maxPeakThrMultiplier, posPeakThrMultiplier);

end
%%
% Run SWTTEO spike detection
params.filter = 0;
in.M = trace;
in.SaRa = fs;
params.method = 'lambda';
params.lambda = 430;
tic
[spikepos1, ~] = SWTTEO(in,params);
toc
[spikeTimes.('swtteo'), ~] = alignPeaks(spikepos1, trace, 10,...
    1, minPeakThrMultiplier, maxPeakThrMultiplier, posPeakThrMultiplier);
%%
% Threshold spike detection
[frames, ~, threshold] = detectSpikesThreshold(trace, 3.5, 0.1, fs, 0);
[spikeTimes.('threshold'), ~] = alignPeaks(find(frames==1), trace, 10,...
    1, minPeakThrMultiplier, maxPeakThrMultiplier, posPeakThrMultiplier);

% Merge spikes from all methods
[spikeTimes.('all'), intersect_matrix, unique_idx] = mergeSpikes(spikeTimes, 'all');

[spikeTimes.('all'), spikeWaveforms] = alignPeaks(spikeTimes.('all'), trace, 10,...
    1, minPeakThrMultiplier, maxPeakThrMultiplier, posPeakThrMultiplier);

% Merge spikes from wavelets
[spikeTimes.('wavelets'), ~, ~] = mergeSpikes(spikeTimes, 'wavelets');


%% Plot all spike markers
figure
set(gcf,'color','w');
h = plot(trace-mean(trace), 'k', 'linewidth', 1);
hold on

admissible = {'threshold', 'wavelets', 'swtteo', 'all'};
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

% Aesthetics
st = 1;
en = st+(60*25);
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
figure
set(gcf,'color','w');
h = plot(trace-mean(trace), 'k', 'linewidth', 1);
hold on



methods = fieldnames(spikeTimes);
% unique_counts = table('rownames',methods);
lineStyles=linspecer(length(methods));
clear unique_counts
unique_counts = table('rownames', methods);


for m = 1:length(methods)
    if ~strcmp(methods{m}, 'all')
    spks = spikeTimes.(methods{m});
    spikepos = find((unique_idx == m)==1);
    scatter(spikepos, repmat(4*std(trace)-(m*threshold/-7), length(spikepos), 1), 'v', 'filled',...
        'markerfacecolor', lineStyles(m,:));
    unique_counts{m,1} = length(spikepos);
    else
    spikepos = spikeTimes.(methods{m});
    scatter(spikepos, repmat(4*std(trace)-(m*threshold/-7), length(spikepos), 1), 'v', 'filled',...
        'markerfacecolor', lineStyles(m,:));
    end
end
unique_counts.Properties.VariableNames = {'No. unique spikes'};
unique_counts(end-1:end-2,:) = [];
hold on
yl = yline(threshold, 'magenta--', "Threshold = " + round(threshold,2) + "\muV");
yl.LabelVerticalAlignment = 'bottom';
yl.LineWidth = 1.5;
yl.FontSize = 10;

% yl = yline(3/3.5*threshold, 'cyan--', "Threshold = " + round(3/3.5*threshold) + "\muV");
% yl.LabelVerticalAlignment = 'top';
% yl.LineWidth = 1.5;
% yl.FontSize = 10;

% Aesthetics
st = 1;
en = st+(60*25);
xlim([st en]);
ylim([-6*std(trace) 5*std(trace)])
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
t = tiledlayout(1, length(methods));
t.Title.String = 'Spikes detected uniquely by each method';
for i = 1:length(methods)-1
    method = methods{i};
    if ~strcmp(method, 'all')
        spk_method = find((unique_idx == m)==1);
        spk_waves_method = spikeWaveforms(:, spk_method);
        nexttile
        plot(spk_waves_method, 'linewidth', 0.1, 'color', [0.7 0.7 0.7])
        hold on
        plot(mean(spk_waves_method,2), 'linewidth', 1.5, 'color', [0 0 0])
        title({[method],["No. unique spikes: " + length(spk_method)]})
        box off;
        axis tight
        ylim([-6*std(trace) 5*std(trace)]);
        ylabel('Voltage [\muV]')
        set(gca, 'xcolor', 'none');
        yline(threshold, 'r--');
    end
end

%% UI functions
function keys(src, event)
global TP FP FN nSamples
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
global duration fs
lims = get(gca, 'xlim');
st = lims(2);
en = st+(60*25);
if en <= duration*fs
    set(gca,'xlim', [st en]);
end
end

function Prev(prevButton, EventData)
lims = get(gca, 'xlim');
en = lims(1);
st = en-(60*25);
if st >= 1
    set(gca,'xlim', [st en]);
end
end

function Save(saveButton, EventData)
% save("MPT200220_3A_DIV21" + "stats.mat", 'TP', 'TN', 'FP', 'nSamples');
end
