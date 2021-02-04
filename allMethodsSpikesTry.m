clearvars; clc;
file = load('200708_slice1_1.mat');
ttx_file = load('200708_slice1_3_TTX.mat');

%%

global fs duration TP FP FN nSamples
fs = 25000;
TP = 0;
FP = 0;
FN = 0;
nSamples = 0;
duration = length(ttx_file.dat)/fs;
channel = 14;
trace_raw = file.dat(1:fs*60, channel);
% trace = ttx_file.dat(1:fs*60, 12);

Ns = 10;
multiplier = 3.0;

 [spikeTimes.('mea'), ~, trace] = detectSpikesCWT(...
    trace_raw, fs, [0.5 1], 'mea', -0.4, Ns, multiplier, 200, 0, ...
    1, 10, 10);

params.filter = 0;
in.M = trace;
in.SaRa = fs;
%%
params.filter = 0;
in.M = trace;
in.SaRa = fs;
params.method = 'auto';
[spikepos1, ~] = SWTTEO(in,params);
spike_times = [];
spikeWaveforms = [];
for i = 1:length(spikepos1)
    if spikepos1(i)>25
        bin = trace(spikepos1(i)-25:spikepos1(i)+25);
        if trace(spikepos1(i))<0
            spikeWaveforms(:,end+1) = bin;
            spike_times(end+1) = spikepos1(i);
        end
    end
end
 [spikeTimes.('swtteo'), ~] = alignPeaks(spike_times, trace, 25,...
                                                    1, 1, 5, 10);

%%

[frames, ~, threshold] = ...
    detectSpikesThreshold(trace, 3.5, 0.1, fs, 0);
[spikeTimes.('threshold'), ~] = alignPeaks(find(frames==1), trace, 10,...
                                                    0, 0, 5, 10);
%%
spikeTimes.all = [];
[spikeTimes.('all'), intersect_matrix, unique_idx] = mergeSpikes(spikeTimes);
[~, spikeWaveforms] = alignPeaks(spikeTimes.('all'), trace, 25,...
                                                    0, 1, 10, 10);
%%
h = plot(trace-mean(trace), 'k');
hold on
methods = fieldnames(spikeTimes);
col = parula(length(methods)+1);
for m = 1:length(methods)
    spikepos = spikeTimes.(methods{m});
    scatter(spikepos, repmat(4*std(trace)-(m+1)*3, length(spikepos), 1), 'v', 'filled',...
        'markerfacecolor', col(m,:))
end
hold on
yline(threshold, 'r--', "Threshold = " + round(threshold,2) + "\muV");
% Aesthetics
st = 1;
en = st+(50*25);
xlim([st en]);
ylim([-6*std(trace) 5*std(trace)])
pbaspect([5 1 1])
box off
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend(methods, 'location', 'northeastoutside')
set(gcf, 'position', [300 300 1500 350]);

% Create buttons and callbacks
set(gcf,'KeyPressFcn',@keys);

nextButton = uicontrol('Position',[700 1 150 30],'String','Next',...
    'Callback', @Next);

prevButton = uicontrol('Position',[500 1 150 30],'String','Previous',...
    'Callback', @Prev);

saveButton = uicontrol('Position',[1325, 120, 100, 30],'String','Save',...
    'Callback', @Save);
%%

t = tiledlayout(1, length(methods)-1);
t.Title.String = 'Spikes detected uniquely by each method';
for i = 1:length(methods)
    method = methods{i};
    if ~strcmp(method, 'all')
        spk_method = find(unique_idx == i);
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
%%
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
en = st+(50*25);
if en <= duration*fs
    set(gca,'xlim', [st en]);
end
end

function Prev(prevButton, EventData)
lims = get(gca, 'xlim');
en = lims(1);
st = en-(50*25);
if st >= 1
    set(gca,'xlim', [st en]);
end
end

function Save(saveButton, EventData)
% save("MPT200220_3A_DIV21" + "stats.mat", 'TP', 'TN', 'FP', 'nSamples');
end
