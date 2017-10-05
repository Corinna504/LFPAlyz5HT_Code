function [fig_h, fig_spect] = LFPGui( ex1, ex2, fname)
%[fig_h, fig_spect] = LFPGui( ex1, ex2, fname)
% gui showing different LFP related plots:
% raw lfp vs. filtered time series
% raw lfp vs. filtered power
% spike triggered average
% stimulus triggered average
% spike field coherence
%
%
%
% ex1 containes the baseline experiment, 
% ex2 containes the drug experiment
% fname is the figure name

%%

tag_number = randi(1000,1);

fig_h = figure('Position', [57 188 1100 442], 'Name', fname);

fig_spect = findobj('tag', ['spectogram' num2str(tag_number)], 'Name', fname);
if isempty(fig_spect)
    fig_spect = figure('tag', ['spectogram' num2str(tag_number)], 'Position', [101   279   989   701], 'Name', fname);
else
    figure(fig_spect, 'Position', [101   279   989   701], 'Name', fname);
end
            
params.Fs = 1000; params.trialave = 1; 

%% specifications
bg = uibuttongroup(fig_h, 'Title','Filter and Multitaper', ...
    'Position', [0.03 0.5 0.18 0.45]);

% notch filter
uicontrol(bg, 'Style','text', 'String','Notch', 'Position', [5 160 50 15], ...
    'HorizontalAlignment', 'left');
notchfreqlb_h = uicontrol(bg, 'Style','edit', 'String', '49', 'Position', [45 160 25 15]);
notchfrequb_h = uicontrol(bg, 'Style','edit', 'String', '51', 'Position', [80 160 25 15]);


uicontrol(bg, 'Style','text', 'String', 'order', 'Position', [120 160 50 15], ...
    'HorizontalAlignment', 'left');
notchorder_h = uicontrol(bg, 'Style','edit', 'String', '2', 'Position', [150 160 30 15]);


% lowpass filter
uicontrol(bg, 'Style','text', 'String','Lowpass', 'Position', [5 130 50 15], ...
    'HorizontalAlignment', 'left');
lowpassfreq_h = uicontrol(bg, 'Style','edit', 'String', '100', 'Position', [55 130 30 15]);

uicontrol(bg, 'Style','text', 'String','order', 'Position', [120 130 50 15], ...
    'HorizontalAlignment', 'left');
lowpassorder_h = uicontrol(bg, 'Style','edit', 'String', '2', 'Position', [150 130 30 15]);


% multitaper
uicontrol(bg, 'Style', 'text', 'String','Multitaper Method', 'Position', [5 80 150 15], ...
    'HorizontalAlignment', 'left');
uicontrol(bg, 'Style', 'text', 'String','nw', 'Position', [5 60 50 15], ...
    'HorizontalAlignment', 'left');
nw_h = uicontrol(bg, 'Style','edit', 'String', '3', 'Position', [55 60 30 15]);

%%frequency range to show
uicontrol(bg, 'Style', 'text', 'String','Frequency Range', ...
    'Position', [5 25 100 15], 'HorizontalAlignment', 'left');
freqrange_lb_h = uicontrol(bg, 'Style','edit', 'String', '1', ...
    'Position', [105 25 30 15]);
freqrange_ub_h = uicontrol(bg, 'Style','edit', 'String', '100', ...
    'Position', [140 25 30 15]);


%%window length and step size for spectogram plots
uicontrol(bg, 'Style', 'text', 'String','Specgram Window', ...
    'Position', [5 5 100 15], 'HorizontalAlignment', 'left');
winlength_h = uicontrol(bg, 'Style','edit', 'String', '0.1', ...
    'Position', [105 5 30 15]);
winstep_h = uicontrol(bg, 'Style','edit', 'String', '0.005', ...
    'Position', [140 5 30 15]);


%% update button
[stimparam, valsB] = getStimParam(ex1);
[~, valsD] = getStimParam(ex2);
vals = intersect(valsD, valsB); % find the common stimuli 

stimcond_h = uicontrol(fig_h, 'Style', 'popupmenu', 'String', {'all'; num2str(vals')},...
    'Position', [120 120 100 30]);


%% popupmenu for comparison panels
uicontrol(fig_h, 'Style', 'Text', 'String','only comparison (lower plots)',...
    'Position', [5 120 100 30]);

uicontrol(fig_h, 'Style', 'pushbutton', 'String','Update', ...
    'Position', [120 80 100 30], 'Callback', @UpdateAxes)



%% axes
ax_lfp_t = axes(fig_h,'Position', [0.3 0.8 0.1 0.15]); 
ax_lfp_f = axes(fig_h,'Position', [0.45 0.8 0.1 0.15]); 
ax_spktrigavg_t = axes(fig_h,'Position', [0.6 0.8 0.1 0.15]);
ax_spkfieldcoh = axes(fig_h,'Position', [0.75 0.8 0.1 0.15]);
% ax_spectogram = axes(fig_h,'Position', [0.9 0.8 0.1 0.15]);


ax_lfp_t2 = axes(fig_h,'Position', [0.3 0.55 0.1 0.15]); 
ax_lfp_f2 = axes(fig_h,'Position', [0.45 0.55 0.1 0.15]);
ax_spktrigavg_t2 = axes(fig_h,'Position', [0.6 0.55 0.1 0.15]);
ax_spkfieldcoh2 = axes(fig_h,'Position', [0.75 0.55 0.1 0.15]);
% ax_spectogram2 = axes(fig_h,'Position', [0.9 0.55 0.1 0.15]);


ax_lfp_t_cmp = axes(fig_h,'Position', [0.3 0.35 0.1 0.1]); 
ax_lfp_f_cmp = axes(fig_h,'Position', [0.45 0.35 0.1 0.1]);
ax_spktrigavg_t_cmp = axes(fig_h,'Position', [0.6 0.35 0.1 0.1]);
ax_spkfieldcoh_cmp = axes(fig_h,'Position', [0.75 0.35 0.2 0.1]);
% ax_spectogram_cmp = axes(fig_h,'Position', [0.9 0.3 0.1 0.15]);



ax_stimulus_lfp_t  = axes(fig_h,'Position', [0.3 0.15 0.1 0.1]);
ax_stimulus_lfp_f  = axes(fig_h,'Position', [0.45 0.15 0.1 0.1]);
ax_stimulus_stafield_f = axes(fig_h,'Position', [0.6 0.15 0.1 0.1]);
ax_stimulus_coh = axes(fig_h,'Position', [0.75 0.15 0.2 0.1]);



UpdateAxes([], [])
%%
    function UpdateAxes(src, evt)
        
        
        k = 2; % #tapers, as in Martin et al. (2016)
        params.tapers = [str2double(nw_h.String), k];
        params.fpass =  [str2double(freqrange_lb_h.String), str2double(freqrange_ub_h.String)];
        win = [str2double(winlength_h.String), str2double(winstep_h.String)];
        CallSpectogramPlot;
        
%         delete(fig_h);
        if ~isempty(ax_lfp_f.Children)
            fprintf('\n\n Updated %s \n', datetime('now'));
            return
        end
        %------------------------------- baseline condition
        ex_base = ex1;
        
        ex_base = frequAnalysis(ex_base,  ...
            'notchf', [str2double(notchfreqlb_h.String)  str2double(notchfrequb_h.String)],...
            'notchord', str2double(notchorder_h.String), ...
            'lowpf', str2double(lowpassfreq_h.String), ...
            'lowpord', str2double(lowpassorder_h.String), ...
            'nw', str2double(nw_h.String));
        
        % LFP time domain
        axes(ax_lfp_t); hold off;   
        lfpave_base = lfpTimeDomain(ex_base);
        
        % LFP frequency domain
        axes(ax_lfp_f); hold off;   
        powave_base = lfpFreqDomain(ex_base, params.fpass);
        
        % Spike Triggered Average
        axes(ax_spktrigavg_t); hold off; 
        sta_base = spktriglfp(ex_base, 'plot', true);
        ylim auto

        % Spike Field Coherence
        axes(ax_spkfieldcoh); hold off; 
        [coh_base, fcoh_base] = spkfieldcoh(ex_base, params); xlim(params.fpass);
        
       
        
        %------------------------------- drug condition
        
        ex_drug = ex2;
        
        ex_drug = frequAnalysis(ex_drug, ...
            'notchf', [str2double(notchfreqlb_h.String)  str2double(notchfrequb_h.String)],...
            'notchord', str2double(notchorder_h.String), ...
            'lowpf', str2double(lowpassfreq_h.String), ...
            'lowpord', str2double(lowpassorder_h.String), ...
            'nw', str2double(nw_h.String));
        
        % LFP time domain
        axes(ax_lfp_t2); hold off;   
        lfpave_drug = lfpTimeDomain(ex_drug);
        
        % LFP frequency domain
        axes(ax_lfp_f2); hold off;   
        powave_drug = lfpFreqDomain(ex_drug, params.fpass);
        
        % Spike Triggered Average
        axes(ax_spktrigavg_t2); hold off; 
        sta_drug = spktriglfp( ex_drug, 'plot', true); ylim auto
        
        
        % Spike Field Coherence
        axes(ax_spkfieldcoh2); hold off; 
        [coh_drug, fcoh_drug] = spkfieldcoh(ex_drug, params); xlim(params.fpass);
        
        
        
        
        %---------------------------------- comparison
        
        if strcmp(stimcond_h.String(stimcond_h.Value), 'all')
            par = ones(length(ex_base.Trials),1); par = num2cell(par);
            [ex_base.Trials.(stimparam)] = deal( par{:});
            
            par = ones(length(ex_drug.Trials),1); par = num2cell(par);
            [ex_drug.Trials.(stimparam)] = deal( par{:});
            
        else
            stimval = str2double(stimcond_h.String{stimcond_h.Value});
            ex_base.Trials = ex_base.Trials([ex_base.Trials.(stimparam)]==stimval);
            ex_drug.Trials = ex_drug.Trials([ex_drug.Trials.(stimparam)]==stimval);
        end
        
        
        % LFP time domain
        axes(ax_lfp_t_cmp); hold off;   
        lfpTimeDomain(ex_base); 
        delete(findobj(ax_lfp_t_cmp, 'LineWidth', 2));
        set(findobj(ax_lfp_t_cmp, 'type', 'line'), 'Color', 'k');
        lfpTimeDomain(ex_drug);
        
        % LFP frequency domain
        axes(ax_lfp_f_cmp); hold off;   
        lfpFreqDomain(ex_base, params.fpass); 
        set(findobj(ax_lfp_f_cmp,'type', 'line'), 'Color', 'k');
        lfpFreqDomain(ex_drug, params.fpass);
        
        
        % Spike Field Coherence
        axes(ax_spkfieldcoh_cmp); hold off; 
        spkfieldcoh(ex_base, params); hold on;
        set(findobj(ax_spkfieldcoh_cmp,'type', 'line'), 'Color', 'k');
        spkfieldcoh(ex_drug, params); 
        xlim(params.fpass);
        legend('baseline', 'drug', 'Location', 'EastOutside');
        
        
        % Spike triggered average
        axes(ax_spktrigavg_t_cmp); hold off; 
        spktriglfp( ex_base, 'plot', true);hold on
        set(findobj(ax_spktrigavg_t_cmp,'type', 'line'), 'Color', 'k');
        hold on
        spktriglfp( ex_drug, 'plot', true);
        delete(findobj(ax_spktrigavg_t_cmp, 'Type', 'Patch'));
        crossl;
        title('blue: raw     filtered:red'); ylim auto
        
        %----------------------------- stimulus averaged values
        
        % stimulus vs lfp signal
        axes(ax_stimulus_lfp_t);
        x = 1:length(valsB);
        xticklabel = num2str(valsB');
        plot(x, nanmean(lfpave_base, 2), '.-k', x, nanmean(lfpave_drug, 2), '.-b');
        set(gca, 'XTick', x, 'XTickLabel', xticklabel, 'XTickLabelRotation', 45 );
        xlabel('stimulus'); 
        ylabel('averaged signal');
        
        % stimulus vs lfp pow
        axes(ax_stimulus_lfp_f);
        plot(x, nanmean(powave_base, 2), '.-k', ...
            x, nanmean(powave_drug, 2), '.-b');
        set(gca, 'XTick', x, 'XTickLabel', xticklabel, 'XTickLabelRotation', 45 );
        xlabel('stimulus'); 
        ylabel('averaged signal');

        % coherence
        axes(ax_stimulus_coh);
        plot(x, nanmean(coh_base, 2), '.-k', ...
            x, nanmean(coh_drug, 2), '.-b');
        set(gca, 'XTick', x, 'XTickLabel', xticklabel, 'XTickLabelRotation', 45 );

        legend('baseline', 'drug', 'Location', 'EastOutside');

        %----------------------------- 
        set(findobj('Type', 'Axes'), 'FontSize', 8);

        fprintf('\n\n Updated %s \n', datetime('now'));
        
    end


    function CallSpectogramPlot
        
        %----------------------------- LFP spectogram
        win = [str2double(winlength_h.String), str2double(winstep_h.String)];

        figure(fig_spect);
        
        exlfp_base = ex1;
        exlfp_base = frequAnalysis(exlfp_base, ...
            'notchf', [str2double(notchfreqlb_h.String)  str2double(notchfrequb_h.String)],...
            'notchord', str2double(notchorder_h.String), ...
            'lowpf', str2double(lowpassfreq_h.String), ...
            'lowpord', str2double(lowpassorder_h.String), ...
            'nw',str2double(nw_h.String));
        ex_drug = ex2;
%         ex_drug.Trials = getPartialTrials(ex_drug.Trials);

        ex_drug = frequAnalysis(ex_drug, ...
            'notchf', [str2double(notchfreqlb_h.String)  str2double(notchfrequb_h.String)],...
            'notchord', str2double(notchorder_h.String), ...
            'lowpf', str2double(lowpassfreq_h.String), ...
            'lowpord', str2double(lowpassorder_h.String), ...
            'nw', str2double(nw_h.String));
        
                   
       % other stimuli
        for i = 1:length(vals)
          
            % baseline
            s(1,i) = subplot(3, length(vals), i);
            S = lfpspecgram(exlfp_base, win, params, i);

               
            % drug
            s(2,i) = subplot(3, length(vals), i+length(vals));
            [S_drug,t,f] = lfpspecgram(ex_drug, win, params, i);
        
            % time spectrum
            s(3,i) = subplot(3, length(vals), i+length(vals)*2);
            
            % I am not so sure whether I should use the difference of
            % log-spaced power or the ratio of linearly spaced power here.
            % Since the power is shown in dB, I think, the
            % difference is more appropriate.
            delta_S(:,:,i) = mag2db(S_drug) - mag2db(S); 
            plot_matrix(delta_S(:,:,i), t, f, 'n');   
            ylabel('Frequency'); xlabel('time [s]');
            colormap('jet');
            
            smax = max(max(abs(delta_S(:,:,i))));
            caxis([-smax smax]); % center the axis
            xlim([t(1) t(end)]);

        end
        
        s(1,1).YLabel.String = char('Baseline stimulus triggered', 'Spectrogram (dB)');
        s(2,1).YLabel.String = char('Drug stimulus triggered', 'Spectrogram (dB)');
        s(3,1).Title.String = 'pow drug (dB) - pow base (dB)';
        
        set(findobj('Type', 'Axes'), 'FontSize', 8);
        
    end


end

