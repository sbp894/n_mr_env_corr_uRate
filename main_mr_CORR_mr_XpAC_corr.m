% Look at modulation filtered meanrates.

chinIDs=[321 322 325 338 341 343];
% function plot_mean_rates_per_chin(chinIDs)
hanModFilter=174;
hanAllModFreqCorr=19;
hanOneModFreqCorr=31;
AllmodFreqs=2.^(0:8);
fSize=14;
mrkSize=12;

saveAllFigs=1;
loading_Dir='/media/parida/DATAPART1/Matlab/SNRenv/n_sEPSM/OUTPUT/DataAnal/';
saving_Dir='/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/mr_corr_OUTPUT/';
Latex_Dir='/home/parida/Dropbox/Study Stuff/Presentations/TorstenPurdueVisit/Figures/';

plotWindowWeights=0;
plotAllSNRrate=0;
plot_modFiltCompare=0;
fontSize=16;

if ~isdir(saving_Dir)
    mkdir(saving_Dir);
end

meanrate_binwidth=.1e-3; % good for mod freq upto 500 Hz

modFreqWeights=[1  1  1  1  1   1   1    1    1    1];
% ------------- 1  2  4  8 16   32  64  128  256  512

if ~isdir(saving_Dir)
    mkdir(saving_Dir);
end

all_ChinSpikeData=[];
chin_snr_track_unit_mat=[];

for chinVar=1:length(chinIDs)
    curChinID=chinIDs(chinVar);
    curDirMeta=dir(sprintf('%s*%d*',loading_Dir,curChinID));
    if isempty(curDirMeta)
        warning('Skipping chin %d. No directory found. Run data_SNRenv_analysis on this chin.', curChinID);
    else
        if length(curDirMeta)>1
            [~, dir_ind2use]=sort([curDirMeta.datenum], 'descend');
            dir_ind2use=dir_ind2use(1);
            curDataDir=[loading_Dir curDirMeta(dir_ind2use).name filesep];
            warning('Multiple directories found. Using the latest dir (%s) for chin %d. ', curDataDir, curChinID);
        else
            curDataDir=[loading_Dir curDirMeta.name filesep];
        end
        
        load([curDataDir 'SpikeStimulusData.mat']);
        load([curDataDir 'ExpControlParams.mat']);
        all_ChinSpikeData=[all_ChinSpikeData, spike_data];  %#ok<*AGROW>
        chin_snr_track_unit_mat=[chin_snr_track_unit_mat; [repmat(curChinID, length(spike_data),1), [spike_data.SNR]', [spike_data.track]', [spike_data.unit]']];
    end
end

unique_chin_snr_track_unit_mat=unique(chin_snr_track_unit_mat, 'rows');
mr_corr_Data=repmat(struct('SSNcorr', nan, 'FLNcorr', nan, 'SSNcorrNF', nan, 'FLNcorrNF', nan, 'CF_Hz', nan, 'SR', nan, 'SNR', nan, 'FLNcorrFinal', nan, 'FLNcorrNFFinal', nan, ...
    'SSNcorrFinal', nan, 'SSNcorrNFFinal', nan), size(unique_chin_snr_track_unit_mat, 1), 1);

parfor plotVar=1:size(unique_chin_snr_track_unit_mat,1)
    cur_inds=find(sum(repmat(unique_chin_snr_track_unit_mat(plotVar,:), size(chin_snr_track_unit_mat,1), 1)==chin_snr_track_unit_mat,2)==size(chin_snr_track_unit_mat,2));
    
    for indVar=1:length(cur_inds)
        
        curMetaData=all_ChinSpikeData(cur_inds(indVar));
        curSpikeData=curMetaData.SpikeTrains;
        fName=sprintf('data_t%d_u%02d_CF_%1.1fk_SNR%i', curMetaData.track, curMetaData.unit, curMetaData.CF_Hz/1e3, curMetaData.SNR);
        minBinRate=max(curMetaData.SR*min(meanrate_binwidth), .1);
        
        
        [S_stim_plus, fs_stim]=audioread(curMetaData.StimsFNames{1,1}{1});
        t_stim=(1:length(S_stim_plus))/fs_stim;
        dur_stim=length(S_stim_plus)/fs_stim;
        
        %% Load stimuli
        [S_stim_plus, ~]=audioread(curMetaData.StimsFNames{1,1}{1});
        [N_stim_plus, ~]=audioread(curMetaData.StimsFNames{2,1}{1});
        [SN_stim_plus, fs_stim]=audioread(curMetaData.StimsFNames{3,1}{1});
        
        %% Load spikes
        S_rate_plus=sort(cell2mat(curSpikeData{1,1}));
        S_spike_minus=sort(cell2mat(curSpikeData{1,2}));
        
        N_spike_plus=sort(cell2mat(curSpikeData{2,1}));
        N_spike_minus=sort(cell2mat(curSpikeData{2,2}));
        
        SN_spike_plus=sort(cell2mat(curSpikeData{3,1}));
        SN_spike_minus=sort(cell2mat(curSpikeData{3,2}));
        
        %% Remove spikes outside of stimulus
        S_rate_plus(S_rate_plus>dur_stim)=[];
        S_spike_minus(S_spike_minus>dur_stim)=[];
        
        N_spike_plus(N_spike_plus>dur_stim)=[];
        N_spike_minus(N_spike_minus>dur_stim)=[];
        
        SN_spike_plus(SN_spike_plus>dur_stim)=[];
        SN_spike_minus(SN_spike_minus>dur_stim)=[];
        
        %% Create rate vector using histogram
        rate_binWidth=min(meanrate_binwidth);
        hist_edges=0:rate_binWidth:dur_stim;
        t_hist=mean([hist_edges(1:end-1); hist_edges(2:end)],1);
        
        S_rate_plus=histcounts(S_rate_plus,hist_edges);
        S_rate_minus=histcounts(S_spike_minus,hist_edges);
        S_rate_env=(S_rate_plus+S_rate_minus)/2;
        
        N_rate_plus=histcounts(N_spike_plus,hist_edges);
        N_rate_minus=histcounts(N_spike_minus,hist_edges);
        N_rate_env=(N_rate_plus+N_rate_minus)/2;
        
        SN_rate_plus=histcounts(SN_spike_plus,hist_edges);
        SN_rate_minus=histcounts(SN_spike_minus,hist_edges);
        SN_rate_env=(SN_rate_plus+SN_rate_minus)/2;
        
        if strcmp(curMetaData.noise, 'FLN')
            [mr_corr_Data(plotVar).FLNcorr, mr_corr_Data(plotVar).FLNcorrNF]=multires_modulation_XpAC_anal(S_rate_plus, S_rate_minus, N_rate_plus, N_rate_minus, ...
                SN_rate_plus, SN_rate_minus, 1/meanrate_binwidth, AllmodFreqs);
        elseif strcmp(curMetaData.noise, 'SSN')
            [mr_corr_Data(plotVar).SSNcorr, mr_corr_Data(plotVar).SSNcorrNF]=multires_modulation_XpAC_anal(S_rate_plus, S_rate_minus, N_rate_plus, N_rate_minus, ...
                SN_rate_plus, SN_rate_minus, 1/meanrate_binwidth, AllmodFreqs);
        end
        mr_corr_Data(plotVar).CF_Hz=curMetaData.CF_Hz;
        mr_corr_Data(plotVar).SR=curMetaData.SR;
        mr_corr_Data(plotVar).SNR=curMetaData.SNR;
        mr_corr_Data(plotVar).FLNcorrFinal=nanmean([mr_corr_Data(plotVar).FLNcorr]);
        mr_corr_Data(plotVar).FLNcorrNFFinal=nanmean([mr_corr_Data(plotVar).FLNcorrNF]);
        mr_corr_Data(plotVar).SSNcorrFinal=nanmean([mr_corr_Data(plotVar).SSNcorr]);
        mr_corr_Data(plotVar).SSNcorrNFFinal=nanmean([mr_corr_Data(plotVar).SSNcorrNF]);
    end
end

%%
lw2=2;
lw3=3;
xShift=.1;

figure(hanAllModFreqCorr);
clf;
hold on;
plot(nan, nan, '-dr', 'linew',2)
plot(nan, nan, 'm<--', 'linew',2)
plot(nan, nan, '-', 'color', .5*ones(1,3), 'linew',2)
plot(nan, nan, '-bd', 'linew',2)
plot(nan, nan, 'c>--', 'linew',2)
plot(nan, nan, 'k-', 'linew',2)

snrs=unique([mr_corr_Data.SNR]);
FLN_corr_NF_dPrime=zeros(length(snrs),1);
SSN_corr_NF_dPrime=zeros(length(snrs),1);
for snrVar=1:length(snrs)
    curSNR=snrs(snrVar);
    curSNRinds=find([mr_corr_Data.SNR]==curSNR);
    plot(curSNR*ones(1, length(curSNRinds))-2*xShift, [mr_corr_Data(curSNRinds).FLNcorrFinal], 'bd');
    plot(curSNR*ones(1, length(curSNRinds))+2*xShift, [mr_corr_Data(curSNRinds).SSNcorrFinal], 'rd');
    plot(curSNR*ones(1, length(curSNRinds))-xShift, [mr_corr_Data(curSNRinds).FLNcorrNFFinal], 'c>');
    plot(curSNR*ones(1, length(curSNRinds))+xShift, [mr_corr_Data(curSNRinds).SSNcorrNFFinal], 'm<');
    
    FLN_corr_NF_dPrime(snrVar,1)=mean([mr_corr_Data(curSNRinds).FLNcorrFinal]);
    FLN_corr_NF_dPrime(snrVar,2)=mean([mr_corr_Data(curSNRinds).FLNcorrNFFinal]);
    FLN_corr_NF_dPrime(snrVar,3)=calc_dPrime([mr_corr_Data(curSNRinds).FLNcorrFinal], [mr_corr_Data(curSNRinds).FLNcorrNFFinal]);
    
    SSN_corr_NF_dPrime(snrVar,1)=mean([mr_corr_Data(curSNRinds).SSNcorrFinal]);
    SSN_corr_NF_dPrime(snrVar,2)=mean([mr_corr_Data(curSNRinds).SSNcorrNFFinal]);
    SSN_corr_NF_dPrime(snrVar,3)=calc_dPrime([mr_corr_Data(curSNRinds).SSNcorrFinal], [mr_corr_Data(curSNRinds).SSNcorrNFFinal]);
end
set(gca, 'xtick', snrs);
plot(-2*xShift+snrs, FLN_corr_NF_dPrime(:, 1), 'b-o', 'linew', lw3, 'markersize', mrkSize);
plot(2*xShift+snrs, SSN_corr_NF_dPrime(:, 1), 'r-o', 'linew', lw3, 'markersize', mrkSize);
plot(-xShift+snrs, FLN_corr_NF_dPrime(:, 2), 'oc--', 'linew', lw2, 'markersize', mrkSize);
plot(xShift+snrs, SSN_corr_NF_dPrime(:, 2), 'om--', 'linew', lw2, 'markersize', mrkSize);
grid on;

yyaxis right;
set(gca, 'ycolor', 'k');
plot(snrs, FLN_corr_NF_dPrime(:, 3), 'd-k', 'linew', 2, 'markersize', mrkSize);
plot(snrs, SSN_corr_NF_dPrime(:, 3), 'd-', 'color', .5*ones(1,3), 'linew', 2, 'markersize', mrkSize);

xlabel('SNR (dB)');
yyaxis left, ylabel('Chi (avged across time)');
yyaxis right; ylabel('d-prime');
title('Using 1-256 Hz mod Freq/ mr-sEPSMcorr');
legend('SSN', 'SSN-NF', 'SSN-dPrime', 'FLN', 'FLN-NF', 'FLN-dPrime', 'location', 'northwest');

set(hanAllModFreqCorr, 'units', 'inches', 'position', [1 1 8 5]);
set(gca, 'fontsize', fSize);
if saveAllFigs
    saveas(hanAllModFreqCorr, [saving_Dir 'allModFreq_xPACMcorr'], 'tiff');
    saveas(hanAllModFreqCorr, [Latex_Dir 'xPACcorr_allModFreq'], 'epsc');
end
%% check single modualtion frequency
singleModFreq=128;
figure(hanOneModFreqCorr);
clf;
hold on;
plot(nan, nan, '-dr', 'linew',2)
plot(nan, nan, 'm<--', 'linew',2)
plot(nan, nan, '-', 'color', .5*ones(1,3), 'linew',2)
plot(nan, nan, '-bd', 'linew',2)
plot(nan, nan, 'c>--', 'linew',2)
plot(nan, nan, 'k-', 'linew',2)

snrs=unique([mr_corr_Data.SNR]);
FLN_corr_NF_dPrime=zeros(length(snrs),1);
SSN_corr_NF_dPrime=zeros(length(snrs),1);
for snrVar=1:length(snrs)
    curSNR=snrs(snrVar);
    curSNRinds=find([mr_corr_Data.SNR]==curSNR);
    PitchModIndex=dsearchn(AllmodFreqs', singleModFreq);
    
    FLNcorr_pitchRes=[mr_corr_Data(curSNRinds).FLNcorr]; %
    FLNcorr_pitchRes=FLNcorr_pitchRes(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))-2*xShift, FLNcorr_pitchRes, 'bd');
    
    FLNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).FLNcorrNF];
    FLNcorr_pitchRes_nf=FLNcorr_pitchRes_nf(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))-xShift, FLNcorr_pitchRes_nf, 'c>');
    
    SSNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).SSNcorrNF];
    SSNcorr_pitchRes_nf=SSNcorr_pitchRes_nf(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))+xShift,SSNcorr_pitchRes_nf, 'm<');
    
    SSNcorr_pitchRes=[mr_corr_Data(curSNRinds).SSNcorr]; %
    SSNcorr_pitchRes=SSNcorr_pitchRes(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))+2*xShift, SSNcorr_pitchRes, 'rd');
    
    FLN_corr_NF_dPrime(snrVar,1)=mean(FLNcorr_pitchRes);
    FLN_corr_NF_dPrime(snrVar,2)=mean(FLNcorr_pitchRes_nf);
    FLN_corr_NF_dPrime(snrVar,3)=calc_dPrime(FLNcorr_pitchRes, FLNcorr_pitchRes_nf);
    
    SSN_corr_NF_dPrime(snrVar,1)=mean(SSNcorr_pitchRes);
    SSN_corr_NF_dPrime(snrVar,2)=mean(SSNcorr_pitchRes_nf);
    SSN_corr_NF_dPrime(snrVar,3)=calc_dPrime(SSNcorr_pitchRes, SSNcorr_pitchRes_nf);
end
set(gca, 'xtick', snrs);
plot(-2*xShift+snrs, FLN_corr_NF_dPrime(:, 1), 'b-o', 'linew', lw3, 'markersize', mrkSize);
plot(2*xShift+snrs, SSN_corr_NF_dPrime(:, 1), 'r-o', 'linew', lw3, 'markersize', mrkSize);
plot(-xShift+snrs, FLN_corr_NF_dPrime(:, 2), 'oc--', 'linew', lw2, 'markersize', mrkSize);
plot(xShift+snrs, SSN_corr_NF_dPrime(:, 2), 'om--', 'linew', lw2, 'markersize', mrkSize);
grid on;

yyaxis right;
set(gca, 'ycolor', 'k');
plot(snrs, FLN_corr_NF_dPrime(:, 3), 'd-k', 'linew', 2, 'markersize', mrkSize);
plot(snrs, SSN_corr_NF_dPrime(:, 3), 'd-', 'color', .5*ones(1,3), 'linew', 2, 'markersize', mrkSize);

xlabel('SNR (dB)');
yyaxis left, ylabel('Chi (avged across time)');
yyaxis right; ylabel('d-prime');
title(sprintf('mr-sEPSM (corr) w/ one Mod Freq=%.0f Hz', singleModFreq));
legend('SSN', 'SSN-NF', 'SSN-dPrime', 'FLN', 'FLN-NF', 'FLN-dPrime', 'location', 'northwest');

set(hanOneModFreqCorr, 'units', 'inches', 'position', [1 1 8 5]);
set(gca, 'fontsize', fSize);

fName_singleModFreqCorr=sprintf('sEPSM_XpAC_modFreq%.0fHz', singleModFreq);
if saveAllFigs
    saveas(hanOneModFreqCorr, [saving_Dir fName_singleModFreqCorr], 'tiff');
    saveas(hanOneModFreqCorr, [Latex_Dir fName_singleModFreqCorr], 'epsc');
end
