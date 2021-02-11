clearvars; clc; close all;
files = dir('*-0.3*.mat');

for file = 1:length(files)
    
    spike_file_name = files(file).name;
    raw_file_name = [spike_file_name(1:strfind(spike_file_name, '_L_')-1) '.mat'];
    culture_name = strrep(raw_file_name(1:end-4),'_',' ');
    
    load(spike_file_name);
    load(raw_file_name);
    
    lowpass = 600;
    highpass = 8000;
    wn = [lowpass highpass] / (fs / 2);
    filterOrder = 3;
    [b, a] = butter(filterOrder, wn);
    filtered_data = filtfilt(b, a, double(dat));
    
    methods = fieldnames(spikeTimes{1});
    methods = sort(methods);
    
    bin_s = 10;
    fs = spikeDetectionResult.params.fs;
    duration = spikeDetectionResult.params.duration;
    channel = 15;
    while channel == 15
        channel = randi([1,60],1);
    end
    trace = filtered_data(:, channel);
    %%
    close all
    dSampF = 250;
    for i = 1:length(methods)
        spk_vec_all = zeros(1, duration*fs);
        spk_matrix = zeros(60, duration*fs/dSampF);
        method = methods{i};
        for j = 1:60
            spk_times = round(spikeTimes{j}.(method));
            spike_count(j) = length(spk_times);
            spk_vec = zeros(1, duration*fs);
            spk_vec(spk_times) = 1;
            spk_vec_all = spk_vec_all+spk_vec;
            spk_matrix(j,:) = nansum(reshape([spk_vec(:); nan(mod(-numel(spk_vec),dSampF),1)],dSampF,[]));
        end
        spike_counts.(method) = spike_count;
        spike_freq.(method) = spike_count/duration;
        spk_matrix_all.(method) = spk_matrix;
        plot(movmean(spk_vec_all, bin_s*fs))
        hold on
    end
    xticks(linspace(1, duration*fs,10));
    xticklabels(round(linspace(1, duration, 10)));
    pbaspect([2,1,1])
    box off
    legend(strrep(methods, 'p','.'));
    xlabel('Time (s)');
    yticklabels(get(gca, 'ytick')*bin_s)
    ylabel('Spiking frequency (Hz)');
    title(culture_name)
    print(culture_name,'-dpdf','-fillpage');
    %%
    for l = 1:15 %times 4
        close all;
        tiledlayout(4,1)
        for p = 1:4
            nexttile
            cmap = jet;
            colors = round(linspace(1,256,length(methods)));
            plot(trace, 'k-')
            hold on
            
            for m = 1:length(methods)
                method = methods{m};
                spike_train = spikeTimes{channel}.(method);
                
                switch spikeDetectionResult.params.unit
                    case 's'
                        spike_train = spike_train * fs;
                    case 'ms'
                        spike_train = spike_train * fs/1000;
                    case 'frames'
                end
                color = parula(colors(m));
                scatter(spike_train, repmat(max(trace)+length(methods)-m-2, length(spike_train), 1), 15, 'v', 'filled','markerfacecolor',cmap(colors(m),:), 'markeredgecolor', 'k', 'linewidth',0.1);
                
                
            end
            methodsl = strrep(methods, 'p','.');
            legend('Filtered voltage trace', methodsl{:}, 'location','bestoutside');
%             st = randi([1 length(trace)-(25*30)]);
            st = randi([1 length(spike_train)]);
            st = spike_train(st);
            if st+15*25 < length(trace)
            xlim([st-15*25 st+15*25]);
            else
                xlim([st-15*25 inf]);
            end
            ylim([min(trace) max(trace)+6])
            box off
            set(gca,'xcolor','none');
            set(gcf,'position', [0 0 1200 1200]);
            ylabel('Amplitude (\muV)');
            axis fill
            
        end
        print([culture_name ' ' num2str(l)],'-dpdf','-fillpage');
        append_pdfs([culture_name '.pdf'], [culture_name ' ' num2str(l) '.pdf']);
        delete([culture_name ' ' num2str(l) '.pdf']);
    end
    
    
    %%
    for m = 1:length(methods)
        method = methods{m};
        maxF(m) = max(spike_freq.(method));
    end
    maxF = max(maxF);
    %%
    close all;
    tiledlayout(3,3)
    for m = 1:length(methods)
        method = methods{m};
        nexttile
        [f,cbar] = customHeatmap(spike_freq.(method), 'markersize', 100,...
            'cbarLimits', [1, maxF]);
        %     axis square
        title({[''],[strrep(method,'p','.')]});
    end
    set(gcf,'position',[0 0 1000 1000])
    print([culture_name num2str(l+1)],'-dpdf','-fillpage');
    append_pdfs([culture_name '.pdf'], [culture_name num2str(l+1) '.pdf']);
    delete([culture_name ' ' num2str(l+1) '.pdf']);
    %%
    close all;
    t = tiledlayout(length(methods),1);
    
    t.Title.String = culture_name;
    
    t.Title.Interpreter = 'none';
    % cmap = bone;
    % cmap = sort(cmap, 'descend');
    % colormap(cmap)
    for m = 1:length(methods)
        method = methods{m};
        nexttile
        imagesc(spk_matrix_all.(method));
        xlim([1 length(spk_matrix_all.(method))/2])
        axis fill
        box off
        ylabel(strrep(method,'p','.'));
    end
    set(gcf,'position', [0 0 1000 1410]);
    
    
    print([culture_name num2str(l+1)],'-dpdf','-fillpage');
    append_pdfs([culture_name '.pdf'], [culture_name num2str(l+2) '.pdf']);
    delete([culture_name ' ' num2str(l+2) '.pdf']);
    
end