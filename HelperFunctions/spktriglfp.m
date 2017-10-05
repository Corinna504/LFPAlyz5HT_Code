function [avg_stlfp, sem_stlfp, accspk] = spktriglfp(ex, varargin)
% spike triggered average lfp signal 
% [avg_stlfp, sem_stlfp, nspk] = spktriglfp( exSpk, exLFP)
%
% optional arguments are:
% 'time' - with following argument. restricting the lfp to the time before
%           and after a spike event. This also affects the number of spikes
%           because the surrounding LFP must be within the stimulus
%           presentation time. The default setting is +/-300ms.
% 'plot' - without following argument. Plots the spike triggered average
%           LFP
%
% @CL 09.04.2017


% warning(['I narrowed the analysis window to the window beginning 600ms after stimulus onset'...
%     'to 500ms before stimulus end in spktriglfp\getSpks. \n'...
%     'I do this to avoid dominant slow fluctuations in the beginning and the end for RC data. \n'])


p_flag = false;
wnd = 0.30; % window before and after spike event to consider

%%% parse input
k = 1; 
while k<=length(varargin)
    switch varargin{k}
        case 'time'
            wnd = varargin{k+1};
        case 'plot'
            p_flag = true;
    end
    k=k+1;
end

%%% find spikes and estimate the lfp +/- time t around it
accstlfp = []; accspk = 0;

for t = 1:length(ex.Trials)
    [stlfp, nspk] = getstlfp4trial(ex.Trials(t), wnd);
    accstlfp = [accstlfp; stlfp];
    accspk = accspk+length(nspk);
end

%%% compute statistics
avg_stlfp = nanmean(accstlfp, 1); % average of the spike triggered lfp (stLFP)
sem_stlfp = nanstd(accstlfp, 0, 1)./sqrt(accspk); % SEM of the stLFP


%%% plot results
if p_flag
    
    [stimparam, vals] = getStimParam(ex);
    
    a1 = fill([-wnd:0.001:wnd, fliplr(-wnd:0.001:wnd)], ...
        [avg_stlfp - sem_stlfp , fliplr(avg_stlfp + sem_stlfp)], 'b');
    a1.FaceColor = [0.5 0.5 0.5]; a1.FaceAlpha = 0.4;
    a1.EdgeColor = 'w'; a1.EdgeAlpha = 0; hold on;
    
    plot(-wnd:0.001:wnd, avg_stlfp, ...
     'ButtonDownFcn', {@PlotAllSTA, -wnd:0.001:wnd, accstlfp}, ...
     'LineWidth', 2, ...
     'DisplayName', sprintf([stimparam '= %1.3f \n #spk: %1.0f'], vals, accspk),...
     'UserData', getFname(ex));
 
    xlabel('time rel:spike [s]');
    ylabel('avg LFP +/- SEM (\muV)');
    xlim([-wnd wnd]);
    legend('show');
    
end

end


%% Helper
function [stlfp, spikes] = getstlfp4trial(trial, wnd)
% lfp at spk +/- window wnd 

spikes = getSpks(trial, wnd); % spikes within the stimulus presentation time
stlfp = nan(length(spikes), length(-wnd:1/1000:wnd));

for i = 1:length(spikes)
    
    tspk = find(trial.LFP_prepro_time <= spikes(i), 1, 'last');
    tstrt = tspk-(wnd*1000);
    tend = tstrt+size(stlfp,2)-1;
    
    try
        % assign the spike centered LFP to the final matrix
        stlfp(i,:) = trial.LFP_prepro(tstrt:tend);
    catch
        disp('');
        
    end
end

end


function spk = getSpks(trials, wnd)
% spikes within the stimulus presentation time


% time offset for spikes 
awnd_strt = wnd+0.3; %<- 600ms after stimulus onset
awnd_end = -wnd-0.3; %<- 500ms before stimulus end


t_strt = trials.Start - trials.TrialStart;
t_end = t_strt(end)+mean(diff(t_strt));

spk = trials.Spikes( trials.Spikes >= t_strt(1)+awnd_strt & ...
    trials.Spikes <= t_end-awnd_end ) - t_strt(1);
spk = round(spk*1000)/1000;

end


%% Callback Function
function PlotAllSTA(source, ~,time, sta)

fprintf('plotting all spike triggered lfp \n')
nspikes = size(sta,1);

figure('Name', source.UserData);

for i = 1:nspikes
    plot(time, sta(i,:)); hold on;
end
set(findobj(gca, 'Type', 'Line'), 'Color', [0 0.2 0.2 0.2], 'LineWidth', 0.5);
plot(time, nanmean(sta), 'k'); hold on;
xlabel('time rel:spike [s]');
ylabel('\muV');

title(sprintf('single spike triggered LFP signals, #spikes %1d', nspikes));
crossl
end
