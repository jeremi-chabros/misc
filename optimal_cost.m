
Ns = 10;
multiplier = 3.5;
wnames = {'mea','bior1.5','bior1.3', 'db2'};
wnames = {'mea'};
Wid = [0.5 1];
L = -0.25;
nSpikes = 200;
ttx = 0;
minPeakThrMultiplier = 3.5;
maxPeakThrMultiplier = 10;
posPeakThrMultiplier = 10;

channel=12;
trace_raw = file.dat(1:fs*60, channel);
trace_ref = file.dat(1:fs*60, 15);
spikeTimes = struct;

cost = [-0.35:0.025:0];
for c = 1:length(cost)
    L = cost(c);
    [spikeTimes(c).mea, ~, ~] = detectSpikesCWT(...
        trace_raw, fs, Wid, wname{1}, L, Ns, multiplier, nSpikes, ttx, ...
        minPeakThrMultiplier, maxPeakThrMultiplier, posPeakThrMultiplier);
   
    [spikeTimes(c).ref, ~, ~] = detectSpikesCWT(...
        trace_ref, fs, Wid, wname{1}, L, Ns, multiplier, nSpikes, ttx, ...
        minPeakThrMultiplier, maxPeakThrMultiplier, posPeakThrMultiplier);
    
    counts_ref(c) = length(spikeTimes(c).ref);
    counts(c) = length(spikeTimes(c).mea);
end
%%
rs = [1, 10, 20, 50, 100, 500, 1000];
for r = rs
    ratio = (counts+r)./(counts_ref+r);
    plot(ratio*r, 'linewidth',1);
    hold on;
    [x,y] = max(ratio);
    s = scatter(y,x*r,'or', 'linewidth',2);
    set(get(get(s,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on;
end
pbaspect([2,1,1]);
% title("Best cost parameter: " + cost(8-y+1));
xlabel('Cost parameter');
xlim([1 inf])
ylabel('$\frac{\textrm{rec}+r}{\textrm{ref}+r}$','interpreter','latex', 'fontsize', 16);
l = legend(num2str(rs(:)));
l.Location = 'bestoutside';
l.Title.String = 'Coefficient';
set(gca, 'xtick', 1:length(cost),...
    'xticklabels', num2str(sort(cost(:),'descend')));
box off;
set(gcf,'color','w');
