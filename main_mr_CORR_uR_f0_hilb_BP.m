% Look at modulation filtered meanrates.
clear;
clc;

chinIDs= [321 322 325 338 341 343 346 347 354 355];
hanOneModFreqCorr_S_SN= 31;
hanOneModFrequnCorr_SN_N=32;

% % % % for nh=0:1
% % % %     if nh
% % % %         chinIDs=[321 322 325 338 341 343 346 347 354 355];
% % % %         hanOneModFreqCorr_S_SN=31;
% % % %         hanOneModFrequnCorr_SN_N=32;
% % % %     else
% % % %         chinIDs= 361;
% % % %         hanOneModFreqCorr_S_SN=1031;
% % % %         hanOneModFrequnCorr_SN_N=1032;
% % % %     end

% function plot_mean_rates_per_chin(chinIDs)
hanModFilter=174;
hanAllModFreqCorr=19;
fSize=14;
leg_fSize=7;
mrkSize=12;

%% Important params
saveAllFigs=0;
N_half_bp=4;
modFreq=128;
% TimeResolution=1/modFreq; % window for correlation.
TimeResolution=32e-3;
combine_chi1_mu0=0;

if combine_chi1_mu0
    postfix_fName= 'chi'; %#ok<*UNRCH>
else
    postfix_fName= 'mu';
end

warning('check time resolution 8 ms vs 20 ms| using %.1f ms now', TimeResolution*1e3);

%%
loading_Dir='/media/parida/DATAPART1/Matlab/SNRenv/n_mr_sEPSM/OUTPUT/DataAnal/';
saving_Dir='/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/mr_corr_OUTPUT/';
Latex_Dir='/home/parida/Dropbox/Study Stuff/Presentations/TorstenPurdueVisit/Figures/';

outFigDir.pdf=sprintf('/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/outFig/uR_f0_hilb_BP_%d/pdf/', modFreq);
if ~isdir(outFigDir.pdf)
    mkdir(outFigDir.pdf);
end
outFigDir.png=sprintf('/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/outFig/uR_f0_hilb_BP_%d/png/', modFreq);
if ~isdir(outFigDir.png)
    mkdir(outFigDir.png);
end

plotWindowWeights=0;
plotAllSNRrate=0;
plot_modFiltCompare=0;
fontSize=16;

if ~isdir(saving_Dir)
    mkdir(saving_Dir);
end

meanrate_binwidth=.2e-3; % good for mod freq upto 1/2/binWidth

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
        if isfield(spike_data, 'thresh') %earlier spike_data created using mr_sEPSM do not have thresh
            spike_data=rmfield(spike_data, 'thresh');
        end
        all_ChinSpikeData=[all_ChinSpikeData, spike_data];  %#ok<*AGROW>
        chin_snr_track_unit_mat=[chin_snr_track_unit_mat; [repmat(curChinID, length(spike_data),1), [spike_data.SNR]', [spike_data.track]', [spike_data.unit]']];
    end
end

unique_chin_snr_track_unit_mat=unique(chin_snr_track_unit_mat, 'rows');
mr_corr_Data=repmat(struct(...
    'SSNcorr_s_sn_pos', nan, 'SSNuncorr_sn_n_pos', nan, 'SSNcorr_s_n_pos', nan, ...
    'FLNcorr_s_sn_pos', nan, 'FLNuncorr_sn_n_pos', nan, 'FLNcorr_s_n_pos', nan, ...
    'SSNcorr_s_sn_neg', nan, 'SSNuncorr_sn_n_neg', nan, 'SSNcorr_s_n_neg', nan, ...
    'FLNcorr_s_sn_neg', nan, 'FLNuncorr_sn_n_neg', nan, 'FLNcorr_s_n_neg', nan, ...
    'CF_Hz', nan, 'SR', nan, 'SNR', nan, 'chinID', nan, 'track', nan', 'unit', nan, ...
    'FLNcorr_s_sn_pos_Final', nan, 'FLNuncorr_sn_n_pos_Final', nan, 'FLNcorr_s_n_pos_Final', nan, ...
    'SSNcorr_s_sn_pos_Final', nan, 'SSNuncorr_sn_n_pos_Final', nan, 'SSNcorr_s_n_pos_Final', nan, ...
    'FLNcorr_s_sn_neg_Final', nan, 'FLNuncorr_sn_n_neg_Final', nan, 'FLNcorr_s_n_neg_Final', nan, ...
    'SSNcorr_s_sn_neg_Final', nan, 'SSNuncorr_sn_n_neg_Final', nan, 'SSNcorr_s_n_neg_Final', nan), ...
    size(unique_chin_snr_track_unit_mat, 1), 1);

%% Main parfor loop
parfor plotVar= 1:size(unique_chin_snr_track_unit_mat,1)
    cur_inds=find(sum(repmat(unique_chin_snr_track_unit_mat(plotVar,:), size(chin_snr_track_unit_mat,1), 1)==chin_snr_track_unit_mat,2)==size(chin_snr_track_unit_mat,2));
    
    for indVar=1:length(cur_inds)
        
        curMetaData=all_ChinSpikeData(cur_inds(indVar)); %#ok<PFBNS>
        curSpikeData=curMetaData.SpikeTrains;
        fName=sprintf('data_t%d_u%02d_CF_%1.1fk_SNR%i', curMetaData.track, curMetaData.unit, curMetaData.CF_Hz/1e3, curMetaData.SNR);
        figName=sprintf('Q%d_t%d_u%02d_SNR%i', curMetaData.chinID, curMetaData.track, curMetaData.unit, curMetaData.SNR);
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
            [mr_corr_Data(plotVar).FLNcorr_s_sn_pos, mr_corr_Data(plotVar).FLNuncorr_sn_n_pos, mr_corr_Data(plotVar).FLNcorr_s_n_pos, ...
                mr_corr_Data(plotVar).FLNcorr_s_sn_neg, mr_corr_Data(plotVar).FLNuncorr_sn_n_neg, mr_corr_Data(plotVar).FLNcorr_s_n_neg]= ...
                multires_modulation_uR_f0_hilb_BP(S_rate_plus, S_rate_minus, N_rate_plus, N_rate_minus, ...
                SN_rate_plus, SN_rate_minus, 1/meanrate_binwidth, modFreq, outFigDir, TimeResolution, combine_chi1_mu0, N_half_bp, [figName '_FLN'], curMetaData.CF_Hz);
        elseif strcmp(curMetaData.noise, 'SSN')
            [mr_corr_Data(plotVar).SSNcorr_s_sn_pos, mr_corr_Data(plotVar).SSNuncorr_sn_n_pos, mr_corr_Data(plotVar).SSNcorr_s_n_pos, ...
                mr_corr_Data(plotVar).SSNcorr_s_sn_neg, mr_corr_Data(plotVar).SSNuncorr_sn_n_neg, mr_corr_Data(plotVar).SSNcorr_s_n_neg]=...
                multires_modulation_uR_f0_hilb_BP(S_rate_plus, S_rate_minus, N_rate_plus, N_rate_minus, ...
                SN_rate_plus, SN_rate_minus, 1/meanrate_binwidth, modFreq, outFigDir, TimeResolution, combine_chi1_mu0, N_half_bp, [figName '_SSN'], curMetaData.CF_Hz);
        end
        mr_corr_Data(plotVar).CF_Hz=curMetaData.CF_Hz;
        mr_corr_Data(plotVar).SR=curMetaData.SR;
        mr_corr_Data(plotVar).SNR=curMetaData.SNR;
        mr_corr_Data(plotVar).chinID=curMetaData.chinID;
        mr_corr_Data(plotVar).track=curMetaData.track;
        mr_corr_Data(plotVar).unit=curMetaData.unit;
        
        % positive uR
        mr_corr_Data(plotVar).FLNcorr_s_sn_pos_Final=nanmean([mr_corr_Data(plotVar).FLNcorr_s_sn_pos]);
        mr_corr_Data(plotVar).FLNuncorr_sn_n_pos_Final=nanmean([mr_corr_Data(plotVar).FLNuncorr_sn_n_pos]);
        mr_corr_Data(plotVar).FLNcorr_s_n_pos_Final=nanmean([mr_corr_Data(plotVar).FLNcorr_s_n_pos]);
        mr_corr_Data(plotVar).SSNcorr_s_sn_pos_Final=nanmean([mr_corr_Data(plotVar).SSNcorr_s_sn_pos]);
        mr_corr_Data(plotVar).SSNuncorr_sn_n_pos_Final=nanmean([mr_corr_Data(plotVar).SSNuncorr_sn_n_pos]);
        mr_corr_Data(plotVar).SSNcorr_s_n_pos_Final=nanmean([mr_corr_Data(plotVar).SSNcorr_s_n_pos]);
        
        % negative uR
        mr_corr_Data(plotVar).FLNcorr_s_sn_neg_Final=nanmean([mr_corr_Data(plotVar).FLNcorr_s_sn_neg]);
        mr_corr_Data(plotVar).FLNuncorr_sn_n_neg_Final=nanmean([mr_corr_Data(plotVar).FLNuncorr_sn_n_neg]);
        mr_corr_Data(plotVar).FLNcorr_s_n_neg_Final=nanmean([mr_corr_Data(plotVar).FLNcorr_s_n_neg]);
        mr_corr_Data(plotVar).SSNcorr_s_sn_neg_Final=nanmean([mr_corr_Data(plotVar).SSNcorr_s_sn_neg]);
        mr_corr_Data(plotVar).SSNuncorr_sn_n_neg_Final=nanmean([mr_corr_Data(plotVar).SSNuncorr_sn_n_neg]);
        mr_corr_Data(plotVar).SSNcorr_s_n_neg_Final=nanmean([mr_corr_Data(plotVar).SSNcorr_s_n_neg]);
    end
end

fName2Save=sprintf('%suRate_uR_f0_hilb_BP_%d.mat', saving_Dir, modFreq);
save(fName2Save, 'mr_corr_Data', 'modFreq');

%% Compute correlation between R(S) and R(SN)
figure(hanOneModFreqCorr_S_SN);
clf;

% subplot 1
subplot(121);
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

lw2=2;
lw3=3;
xShift=.1;

for snrVar=1:length(snrs)
    curSNR=snrs(snrVar);
    curSNRinds=find([mr_corr_Data.SNR]==curSNR);
    PitchModIndex=dsearchn(modFreq', modFreq);
    
    FLNcorr_pitchRes=[mr_corr_Data(curSNRinds).FLNcorr_s_sn_pos]; %
    FLNcorr_pitchRes=FLNcorr_pitchRes(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))-2*xShift, FLNcorr_pitchRes, 'bd');
    
    FLNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).FLNcorr_s_n_pos];
    FLNcorr_pitchRes_nf=FLNcorr_pitchRes_nf(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))-xShift, FLNcorr_pitchRes_nf, 'c>');
    
    SSNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).SSNcorr_s_n_pos];
    SSNcorr_pitchRes_nf=SSNcorr_pitchRes_nf(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))+xShift,SSNcorr_pitchRes_nf, 'm<');
    
    SSNcorr_pitchRes=[mr_corr_Data(curSNRinds).SSNcorr_s_sn_pos]; %
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
if combine_chi1_mu0
    yyaxis left, ylabel('$\overline{\chi ^2}$', 'interpreter', 'latex');
else
    yyaxis left, ylabel('$\overline{R_{f0,BP}(S).R_{f0,BP}(SN)}$', 'interpreter', 'latex');
    ylim([0 1]);
end

yyaxis right; ylabel('d-prime');
title(sprintf('ModFreq=%.0f Hz (corr(S,SN)), +ve uR', modFreq));
% legend('SSN', 'SSN-NF', 'SSN-dPrime', 'FLN', 'FLN-NF', 'FLN-dPrime', 'location', 'northwest', 'FontSize',10);
set(gca, 'fontsize', fSize);

% subplot 2
subplot(122);
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

lw2=2;
lw3=3;
xShift=.1;

for snrVar=1:length(snrs)
    curSNR=snrs(snrVar);
    curSNRinds=find([mr_corr_Data.SNR]==curSNR);
    PitchModIndex=dsearchn(modFreq', modFreq);
    
    FLNcorr_pitchRes=[mr_corr_Data(curSNRinds).FLNcorr_s_sn_neg]; %
    FLNcorr_pitchRes=FLNcorr_pitchRes(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))-2*xShift, FLNcorr_pitchRes, 'bd');
    
    FLNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).FLNcorr_s_n_neg];
    FLNcorr_pitchRes_nf=FLNcorr_pitchRes_nf(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))-xShift, FLNcorr_pitchRes_nf, 'c>');
    
    SSNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).SSNcorr_s_n_neg];
    SSNcorr_pitchRes_nf=SSNcorr_pitchRes_nf(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))+xShift,SSNcorr_pitchRes_nf, 'm<');
    
    SSNcorr_pitchRes=[mr_corr_Data(curSNRinds).SSNcorr_s_sn_neg]; %
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
if combine_chi1_mu0
    yyaxis left, ylabel('$\overline{\chi ^2}$', 'interpreter', 'latex');
else
    yyaxis left, ylabel('$\overline{R_{f0,BP}(S).R_{f0,BP}(SN)}$', 'interpreter', 'latex');
    ylim([0 1]);
end

yyaxis right; ylabel('d-prime');
title(sprintf('ModFreq=%.0f Hz (corr(S,SN)), -ve uR', modFreq));
LG1=legend('SSN', 'SSN-NF', 'SSN-dPrime', 'FLN', 'FLN-NF', 'FLN-dPrime', 'location', 'northwest');
LG1.FontSize=leg_fSize;
set(gca, 'fontsize', fSize);


set(hanOneModFreqCorr_S_SN, 'units', 'inches', 'position', [1 1 14 5]);
fName_singleModFreqCorr=sprintf('sEPSM_uR_f0_hilb_BP_modFreq%.0fHz_corr_%s_tRes%.0fms', modFreq, postfix_fName, TimeResolution*1e3);
if saveAllFigs
    saveas(hanOneModFreqCorr_S_SN, [saving_Dir fName_singleModFreqCorr], 'tiff');
    saveas(hanOneModFreqCorr_S_SN, [Latex_Dir fName_singleModFreqCorr], 'epsc');
end

%% Compute 1- {correlation between R(SN) and R(N)}
figure(hanOneModFrequnCorr_SN_N);
clf;

% subplot 1
subplot(121);
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
    PitchModIndex=dsearchn(modFreq', modFreq);
    
    FLNuncorr_pitchRes=[mr_corr_Data(curSNRinds).FLNuncorr_sn_n_pos]; %
    FLNuncorr_pitchRes=FLNuncorr_pitchRes(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))-2*xShift, FLNuncorr_pitchRes, 'bd');
    
    FLNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).FLNcorr_s_n_pos];
    FLNcorr_pitchRes_nf=FLNcorr_pitchRes_nf(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))-xShift, FLNcorr_pitchRes_nf, 'c>');
    
    SSNuncorr_pitchRes=[mr_corr_Data(curSNRinds).SSNuncorr_sn_n_pos]; %
    SSNuncorr_pitchRes=SSNuncorr_pitchRes(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))+2*xShift, SSNuncorr_pitchRes, 'rd');
    
    SSNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).SSNcorr_s_n_pos];
    SSNcorr_pitchRes_nf=SSNcorr_pitchRes_nf(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))+xShift,SSNcorr_pitchRes_nf, 'm<');
    
    FLN_corr_NF_dPrime(snrVar,1)=mean(FLNuncorr_pitchRes);
    FLN_corr_NF_dPrime(snrVar,2)=mean(FLNcorr_pitchRes_nf);
    FLN_corr_NF_dPrime(snrVar,3)=calc_dPrime(FLNuncorr_pitchRes, FLNcorr_pitchRes_nf);
    
    SSN_corr_NF_dPrime(snrVar,1)=mean(SSNuncorr_pitchRes);
    SSN_corr_NF_dPrime(snrVar,2)=mean(SSNcorr_pitchRes_nf);
    SSN_corr_NF_dPrime(snrVar,3)=calc_dPrime(SSNuncorr_pitchRes, SSNcorr_pitchRes_nf);
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
if combine_chi1_mu0
    yyaxis left, ylabel('$\overline{\chi ^2}$', 'interpreter', 'latex');
else
    yyaxis left, ylabel('$1-\overline{R_{f0,BP}(SN).R_{f0,BP}(N)}$', 'interpreter', 'latex');
    ylim([0 1]);
end
yyaxis right; ylabel('d-prime');
title(sprintf('ModFreq=%.0f Hz (1-corr(SN,N)), +ve uR', modFreq));
% legend('SSN', 'SSN-NF', 'SSN-dPrime', 'FLN', 'FLN-NF', 'FLN-dPrime', 'location', 'northwest', 'FontSize',10);
set(gca, 'fontsize', fSize);

% subplot 2
subplot(122);
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
    PitchModIndex=dsearchn(modFreq', modFreq);
    
    FLNuncorr_pitchRes=[mr_corr_Data(curSNRinds).FLNuncorr_sn_n_neg]; %
    FLNuncorr_pitchRes=FLNuncorr_pitchRes(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))-2*xShift, FLNuncorr_pitchRes, 'bd');
    
    FLNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).FLNcorr_s_n_neg];
    FLNcorr_pitchRes_nf=FLNcorr_pitchRes_nf(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))-xShift, FLNcorr_pitchRes_nf, 'c>');
    
    SSNuncorr_pitchRes=[mr_corr_Data(curSNRinds).SSNuncorr_sn_n_neg]; %
    SSNuncorr_pitchRes=SSNuncorr_pitchRes(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))+2*xShift, SSNuncorr_pitchRes, 'rd');
    
    SSNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).SSNcorr_s_n_neg];
    SSNcorr_pitchRes_nf=SSNcorr_pitchRes_nf(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))+xShift,SSNcorr_pitchRes_nf, 'm<');
    
    FLN_corr_NF_dPrime(snrVar,1)=mean(FLNuncorr_pitchRes);
    FLN_corr_NF_dPrime(snrVar,2)=mean(FLNcorr_pitchRes_nf);
    FLN_corr_NF_dPrime(snrVar,3)=calc_dPrime(FLNuncorr_pitchRes, FLNcorr_pitchRes_nf);
    
    SSN_corr_NF_dPrime(snrVar,1)=mean(SSNuncorr_pitchRes);
    SSN_corr_NF_dPrime(snrVar,2)=mean(SSNcorr_pitchRes_nf);
    SSN_corr_NF_dPrime(snrVar,3)=calc_dPrime(SSNuncorr_pitchRes, SSNcorr_pitchRes_nf);
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
if combine_chi1_mu0
    yyaxis left, ylabel('$\overline{\chi ^2}$', 'interpreter', 'latex');
else
    yyaxis left, ylabel('$1-\overline{R_{f0,BP}(SN).R_{f0,BP}(N)}$', 'interpreter', 'latex');
    ylim([0 1]);
end
yyaxis right; ylabel('d-prime');
title(sprintf('ModFreq=%.0f Hz (1-corr(SN,N)), -ve uR', modFreq));
LG2=legend('SSN', 'SSN-NF', 'SSN-dPrime', 'FLN', 'FLN-NF', 'FLN-dPrime', 'location', 'northwest');
LG2.FontSize=leg_fSize;
set(gca, 'fontsize', fSize);

set(hanOneModFrequnCorr_SN_N, 'units', 'inches', 'position', [1 1 14 5]);
fName_singleModFrequnCorr=sprintf('sEPSM_uR_f0_hilb_BP_modFreq%.0fHz_uncorr_%s_tRes%.0fms', modFreq, postfix_fName, TimeResolution*1e3);

if saveAllFigs
    saveas(hanOneModFrequnCorr_SN_N, [saving_Dir fName_singleModFrequnCorr], 'tiff');
    saveas(hanOneModFrequnCorr_SN_N, [Latex_Dir fName_singleModFrequnCorr], 'epsc');
end
% % % % % end