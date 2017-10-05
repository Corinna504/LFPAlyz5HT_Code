function ex = frequAnalysis( ex, varargin )
% filters and interpolates lfp traces and returns their the ex file also
% containing the filtered, interpolated and power spectra
%
% optional input variables are:
% fs : sampling frequency
% notchf : notch filtered frequency
% notchord : notch filter order
% bpf:   band pass filtered frequency
% bpord: band pass filter order
% nw    : time bandwidth product
%
% @CL 
%


%% define variables
stimdur = getStimDur(ex); % stimulus presentation duration

% notch filter variables
Fs = 1000;              % sampling frequency
notchf = [49 51];     % notch filter frequency1
notchf2 = [99 101];
notchf3 = [149 151];
notchord = 2;           % filter order

bpf = [1 100];     % bandpass filter cutoff frequency
bpord = 5;            % filter order

% power spectrum in Chalk et al. time halfbandwidth was 3 and k was 5
nw = 2;         % time half bandwidth product

nfft = 1024; % number of frequency samples - I think

% time relative to stimulus onset to filter and interpolate LFP
t_off = -0.2;



% TODO: how to handle files with missing LFP entries (such as ma0014,
% 5.44pm)
if any(cellfun(@isempty, {ex.Trials.LFP}))
    warning('TODO: in frequAnalysis line 44. How to proceed with empty LFP entries.');
    ex.Trials = ex.Trials(~cellfun(@isempty, {ex.Trials.LFP}));
end

%% parse input
k = 1;
while k<=length(varargin)
    switch varargin{k}
        case 'fs'
            Fs = varargin{k+1};
        case 'notchf'
            notchf = varargin{k+1};
        case 'notchord'
            notchord = varargin{k+1};
        case 'bpf'
            bpf = varargin{k+1};
        case 'bpord'
            bpord = varargin{k+1};
        case 'nw'
            nw = varargin{k+1};
    end
    k=k+2;
end
clearvars k;

%% generate filters
% define notch filter
[b_notch,a_notch] = butter(notchord, notchf/(Fs/2), 'stop' );

[b_notch2,a_notch2] = butter(bpord, notchf2/(Fs/2), 'stop');
[b_notch3,a_notch3] = butter(bpord, notchf3/(Fs/2), 'stop');

% define bandpass filter
[b_bpass,a_bpass] = butter(bpord, bpf/(Fs/2), 'bandpass');


%% perform functions on each trial lfp
for ind = 1:length(ex.Trials)
   
    t_frame = ex.Trials(ind).Start - ex.Trials(ind).TrialStart; % time of frame onsets
    t_lfp = ex.Trials(ind).LFP_ts - t_frame(1) ; % time rel:stimulus onset

    %%% bandpass filter
    filt1 = filtfilt(b_bpass, a_bpass, ex.Trials(ind).LFP);        

    %%% notch filter - filters line noise
    filt2 = filtfilt(b_notch, a_notch, filt1);
    filt2 = filtfilt(b_notch2, a_notch2, filt2);
    filt2 = filtfilt(b_notch3, a_notch3, filt2);

    
    %     filt2 = ex.Trials(ind).LFP; % for debugging only
    

    %%% reduce the lfp signal to the period of stimulus presentation
    time = t_frame(1)+t_off : 1/Fs : t_frame(1)+stimdur;
    time = time - t_frame(1);
    LFP_proc = interp1( t_lfp, filt2, time );
    
    %%% we are only interested in the stimulu induced fluctuations
    LFP_proc = LFP_proc-mean(LFP_proc);
    
    %%% mutli taper function
    % the PSD returned by pmtm is normalized per frequency unit
    [ex.Trials(ind).POW, ex.Trials(ind).FREQ] = ...
        pmtm( LFP_proc(time>=0), nw, nfft, Fs);
%     [ex.Trials(ind).POW, ex.Trials(ind).FREQ] = ...
%         pwelch( LFP_proc(time>=0), [], [], [], Fs);
     
    % the preprocessed lfp and the corresponding time vector
    ex.Trials(ind).LFP_prepro = LFP_proc; 
    ex.Trials(ind).LFP_prepro_time = time;

end

ex.time = time; % time corresponding to saved lfp signal

end




function avgstimdur = getStimDur(ex)
% returns the averaged and rounded stimulus presentation duration across
% trials


t_frame = cellfun(@(x, y) x-y, {ex.Trials.Start}, {ex.Trials.TrialStart}, ...
    'UniformOutput', 0); % time of frame onsets
stimdur = cellfun(@(x) x(end)+mean(diff(x)) - x(1), t_frame);


% also round the average stimulus duration to 2 digits precision
avgstimdur = round(mean(stimdur), 2);

end

