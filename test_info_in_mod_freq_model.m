clear;
clc;


allSPLs= 55:15:130;
ax=nan(length(allSPLs),1);

% stim params
fs=100e3;
[stim.sig,stim.fs]= audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_08_15-Q361_AN_PTS/Signals/MH/SNRenv/SNR_0/FLN_Stim_S_P.wav');
if size(stim.sig,1)~=1
    stim.sig=stim.sig';
end

if stim.fs~=fs
    stim= gen_resample(stim, fs_stim, fs);
    stim.fs= fs;
end
stim.time= (1:length(stim.sig))/stim.fs;
stim.dur= stim.time(end);

% model parameters
CF    = 2e3;   % CF in Hz;
AN.spont = 40;   % spontaneous firing rate
AN.tabs   = 0.6e-3; % Absolute refractory period
AN.trel   = 0.6e-3; % Baseline mean relative refractory period
AN.cohc  = 1.0;    % normal ohc function
AN.cihc  = 1.0;    % normal ihc function
AN.species = 1;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
AN.noiseType = 1;  % 1 for variable fGn; 0 for fixed (frozen) fGn
AN.implnt = 0;     % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
AN.respDur= 1.1*stim.dur;
ax= nan(length(allSPLs), 1);

% PSTH parameters
nrep = 75;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = .1e-3; % binwidth in seconds;
psthbins = round(psthbinwidth*fs);  % number of psth bins per psth bin
dt=1/fs; %  time step

allPsthPos= cell(length(allSPLs), 1);
allPsthNeg= cell(length(allSPLs), 1);

parfor splVar=1:length(allSPLs)
    SPL= allSPLs(splVar);
    sig= gen_rescale(stim.sig, SPL);
    
    vihc = ANModelBEZ2018.model_IHC_BEZ2018(sig,CF,nrep,dt,AN.respDur,AN.cohc,AN.cihc,AN.species);
    [psth,~, ~, ~, ~,~] = ANModelBEZ2018.model_Synapse_BEZ2018(vihc,CF,nrep,dt,AN.noiseType,AN.implnt,AN.spont,AN.tabs,AN.trel);
    Psth_pos = sum(reshape(psth,psthbins,length(psth)/psthbins));
    
    vihc = ANModelBEZ2018.model_IHC_BEZ2018(-sig,CF,nrep,dt,AN.respDur,AN.cohc,AN.cihc,AN.species);
    [psth,meanrate, ~, ~, ~,~] = ANModelBEZ2018.model_Synapse_BEZ2018(vihc,CF,nrep,dt,AN.noiseType,AN.implnt,AN.spont,AN.tabs,AN.trel);
    Psth_neg= sum(reshape(psth,psthbins,length(psth)/psthbins));
    
    %     figure(1); clf;
    %     ax(1)=subplot(3,1,1);
    %     plot(stim.time, stim.sig);
    %
    %     ax(2)=subplot(3,1,2);
    %     plot((1:length(meanrate))/fs, meanrate);
    
    allPsthPos{splVar}= Psth_pos;
    allPsthNeg{splVar}= Psth_neg;
    
end
psthTime= (1:length(allPsthPos{1}))*psthbinwidth;

%% plot uRs
lw= 2;
for splVar=1:length(allSPLs)
    SPL= allSPLs(splVar);
    ax(splVar)=subplot(3,2,splVar);
    plot(psthTime, allPsthPos{splVar}, psthTime, -allPsthNeg{splVar}, 'linew', lw);
    title(num2str(SPL));
end


subplot(3,2,5);
xlabel('time');
linkaxes(ax, 'x');
xlim([.38 .48]);