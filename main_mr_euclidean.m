% Look at modulation filtered meanrates.
clear;
clc;
hanOneModFreq_SNR= 101;

chinIDs= [321 322 325 338 341 343 346 347 354 355];
% chinIDs=[361];

% % % % for nh=0:1
% % % %     if nh
% % % %         chinIDs=[321 322 325 338 341 343 346 347 354 355];
% % % %         hanOneModFreqdist_S_SN=31;
% % % %         hanOneModFreqdist_SN_N=32;
% % % %     else
% % % %         chinIDs= 361;
% % % %         hanOneModFreqdist_S_SN=1031;
% % % %         hanOneModFreqdist_SN_N=1032;
% % % %     end

% function plot_mean_rates_per_chin(chinIDs)
hanModFilter=174;
hanAllModFreqdist=19;

%% Important params
saveAllFigs=0;
N_half_bp=4;
modFreq=128;
% TimeResolution=1/modFreq; % window for distelation.
TimeResolution=50e-3;
combine_chi1_mu0=0;

if combine_chi1_mu0
    postfix_fName= 'chi'; %#ok<*UNRCH>
else
    postfix_fName= 'mu';
end

warning('check time resolution 8 ms vs 20 ms| using %.1f ms now', TimeResolution*1e3);

%%
loading_Dir='/media/parida/DATAPART1/Matlab/SNRenv/n_mr_sEPSM/OUTPUT/DataAnal/';
saving_Dir='/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/mr_dist_OUTPUT/';
Latex_Dir='/home/parida/Dropbox/Study Stuff/Presentations/TorstenPurdueVisit/Figures/';

outFigDir.pdf=sprintf('/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/outFig/uR_f0_euclidean_%d/pdf/', modFreq);
if ~isdir(outFigDir.pdf)
    mkdir(outFigDir.pdf);
end
outFigDir.png=sprintf('/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/outFig/uR_f0_euclidean_%d/png/', modFreq);
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
mr_dist_Data=repmat(struct(...
    'SSNdist_s_sn_pos', nan, 'SSNdist_sn_n_pos', nan, 'SSNdist_s_n_pos', nan, 'SSNdist_snr_pos', nan, ...
    'FLNdist_s_sn_pos', nan, 'FLNdist_sn_n_pos', nan, 'FLNdist_s_n_pos', nan, 'FLNdist_snr_pos', nan,  ...
    'SSNdist_s_sn_neg', nan, 'SSNdist_sn_n_neg', nan, 'SSNdist_s_n_neg', nan, 'SSNdist_snr_neg', nan,  ...
    'FLNdist_s_sn_neg', nan, 'FLNdist_sn_n_neg', nan, 'FLNdist_s_n_neg', nan, 'FLNdist_snr_neg', nan,  ...
    'CF_Hz', nan, 'SR', nan, 'SNR', nan, 'chinID', nan, 'track', nan', 'unit', nan, ...
    'FLNdist_s_sn_pos_Final', nan, 'FLNdist_sn_n_pos_Final', nan, 'FLNdist_s_n_pos_Final', nan, 'SSNdist_snr_pos_Final', nan, ...
    'SSNdist_s_sn_pos_Final', nan, 'SSNdist_sn_n_pos_Final', nan, 'SSNdist_s_n_pos_Final', nan, 'FLNdist_snr_pos_Final', nan, ...
    'FLNdist_s_sn_neg_Final', nan, 'FLNdist_sn_n_neg_Final', nan, 'FLNdist_s_n_neg_Final', nan, 'SSNdist_snr_neg_Final', nan, ...
    'SSNdist_s_sn_neg_Final', nan, 'SSNdist_sn_n_neg_Final', nan, 'SSNdist_s_n_neg_Final', nan, 'FLNdist_snr_neg_Final', nan), ...
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
            [mr_dist_Data(plotVar).FLNdist_s_sn_pos, mr_dist_Data(plotVar).FLNdist_sn_n_pos, mr_dist_Data(plotVar).FLNdist_s_n_pos, ...
                mr_dist_Data(plotVar).FLNdist_s_sn_neg, mr_dist_Data(plotVar).FLNdist_sn_n_neg, mr_dist_Data(plotVar).FLNdist_s_n_neg, ...
                mr_dist_Data(plotVar).FLNdist_snr_pos, mr_dist_Data(plotVar).FLNdist_snr_neg]= ...
                multires_modulation_uR_euclidean(S_rate_plus, S_rate_minus, N_rate_plus, N_rate_minus, ...
                SN_rate_plus, SN_rate_minus, 1/meanrate_binwidth, modFreq, outFigDir, TimeResolution, combine_chi1_mu0, [figName '_FLN'], curMetaData.CF_Hz);
        elseif strcmp(curMetaData.noise, 'SSN')
            [mr_dist_Data(plotVar).SSNdist_s_sn_pos, mr_dist_Data(plotVar).SSNdist_sn_n_pos, mr_dist_Data(plotVar).SSNdist_s_n_pos, ...
                mr_dist_Data(plotVar).SSNdist_s_sn_neg, mr_dist_Data(plotVar).SSNdist_sn_n_neg, mr_dist_Data(plotVar).SSNdist_s_n_neg, ...
                mr_dist_Data(plotVar).SSNdist_snr_pos, mr_dist_Data(plotVar).SSNdist_snr_neg]=...
                multires_modulation_uR_euclidean(S_rate_plus, S_rate_minus, N_rate_plus, N_rate_minus, ...
                SN_rate_plus, SN_rate_minus, 1/meanrate_binwidth, modFreq, outFigDir, TimeResolution, combine_chi1_mu0, [figName '_SSN'], curMetaData.CF_Hz);
        end
        mr_dist_Data(plotVar).CF_Hz=curMetaData.CF_Hz;
        mr_dist_Data(plotVar).SR=curMetaData.SR;
        mr_dist_Data(plotVar).SNR=curMetaData.SNR;
        mr_dist_Data(plotVar).chinID=curMetaData.chinID;
        mr_dist_Data(plotVar).track=curMetaData.track;
        mr_dist_Data(plotVar).unit=curMetaData.unit;
        
        % positive uR
        mr_dist_Data(plotVar).FLNdist_s_sn_pos_Final=nanmean([mr_dist_Data(plotVar).FLNdist_s_sn_pos]);
        mr_dist_Data(plotVar).FLNdist_sn_n_pos_Final=nanmean([mr_dist_Data(plotVar).FLNdist_sn_n_pos]);
        mr_dist_Data(plotVar).FLNdist_s_n_pos_Final=nanmean([mr_dist_Data(plotVar).FLNdist_s_n_pos]);
        mr_dist_Data(plotVar).FLNdist_snr_pos_Final=nanmean([mr_dist_Data(plotVar).FLNdist_snr_pos].^2);
        
        mr_dist_Data(plotVar).SSNdist_s_sn_pos_Final=nanmean([mr_dist_Data(plotVar).SSNdist_s_sn_pos]);
        mr_dist_Data(plotVar).SSNdist_sn_n_pos_Final=nanmean([mr_dist_Data(plotVar).SSNdist_sn_n_pos]);
        mr_dist_Data(plotVar).SSNdist_s_n_pos_Final=nanmean([mr_dist_Data(plotVar).SSNdist_s_n_pos]);
        mr_dist_Data(plotVar).SSNdist_snr_pos_Final=nanmean([mr_dist_Data(plotVar).SSNdist_snr_pos].^2);
        
        % negative uR
        mr_dist_Data(plotVar).FLNdist_s_sn_neg_Final=nanmean([mr_dist_Data(plotVar).FLNdist_s_sn_neg]);
        mr_dist_Data(plotVar).FLNdist_sn_n_neg_Final=nanmean([mr_dist_Data(plotVar).FLNdist_sn_n_neg]);
        mr_dist_Data(plotVar).FLNdist_s_n_neg_Final=nanmean([mr_dist_Data(plotVar).FLNdist_s_n_neg]);
        mr_dist_Data(plotVar).FLNdist_snr_neg_Final=nanmean([mr_dist_Data(plotVar).FLNdist_snr_neg].^2);

        mr_dist_Data(plotVar).SSNdist_s_sn_neg_Final=nanmean([mr_dist_Data(plotVar).SSNdist_s_sn_neg]);
        mr_dist_Data(plotVar).SSNdist_sn_n_neg_Final=nanmean([mr_dist_Data(plotVar).SSNdist_sn_n_neg]);
        mr_dist_Data(plotVar).SSNdist_s_n_neg_Final=nanmean([mr_dist_Data(plotVar).SSNdist_s_n_neg]);
        mr_dist_Data(plotVar).SSNdist_snr_neg_Final=nanmean([mr_dist_Data(plotVar).SSNdist_snr_neg].^2);
    end
end

fName2Save=sprintf('%suRate_uR_f0_eucl_SNR_%d.mat', saving_Dir, modFreq);
save(fName2Save, 'mr_dist_Data', 'modFreq');

%% Compute {correlation between R(SN) and R(S)}/{correlation between R(SN) and R(N)}
lw2=2;
lw3=3;
xShift=.1;
mrkSize=12;
fSize=14;
leg_fSize=7;


figure(hanOneModFreq_SNR);
clf;

% subplot 1
subplot(121);
hold on;
subplot(122);
hold on;


snrs=unique([mr_dist_Data.SNR]);
FLN_snr_eucl_pos_mean=zeros(length(snrs),1);
SSN_snr_eucl_pos_mean=zeros(length(snrs),1);
FLN_snr_eucl_neg_mean=zeros(length(snrs),1);
SSN_snr_eucl_neg_mean=zeros(length(snrs),1);

for snrVar=1:length(snrs)
    curSNR=snrs(snrVar);
    curSNRinds=find([mr_dist_Data.SNR]==curSNR);
    
    subplot(121);
    
    FLN_sn_ratio_pos= [mr_dist_Data(curSNRinds).FLNdist_sn_n_pos_Final];
    plot(curSNR*ones(1, length(curSNRinds))-2*xShift, FLN_sn_ratio_pos, 'bd');
    
    SSN_sn_ratio_pos= [mr_dist_Data(curSNRinds).SSNdist_sn_n_pos_Final];
    plot(curSNR*ones(1, length(curSNRinds))+2*xShift, SSN_sn_ratio_pos, 'rd');
    
    FLN_snr_eucl_pos_mean(snrVar)= mean(FLN_sn_ratio_pos);
    SSN_snr_eucl_pos_mean(snrVar)= mean(SSN_sn_ratio_pos);
    
    subplot(122);
    

    FLN_sn_ratio_neg= [mr_dist_Data(curSNRinds).FLNdist_sn_n_neg_Final];
    plot(curSNR*ones(1, length(curSNRinds))-2*xShift, FLN_sn_ratio_neg, 'bd');
    
    SSN_sn_ratio_neg= [mr_dist_Data(curSNRinds).SSNdist_sn_n_neg_Final];
    plot(curSNR*ones(1, length(curSNRinds))+2*xShift, SSN_sn_ratio_neg, 'rd');
    
    FLN_snr_eucl_neg_mean(snrVar)= mean(FLN_sn_ratio_neg);
    SSN_snr_eucl_neg_mean(snrVar)= mean(SSN_sn_ratio_neg);
end


subplot(121);
set(gca, 'xtick', snrs);
plot(-2*xShift+snrs, FLN_snr_eucl_pos_mean, 'b-o', 'linew', lw3, 'markersize', mrkSize);
plot(2*xShift+snrs, SSN_snr_eucl_pos_mean, 'r-o', 'linew', lw3, 'markersize', mrkSize);
grid on;
xlabel('SNR (dB)');
title(sprintf('ModFreq=%.0f Hz (corr(SN,S)/corr(SN,N)), +ve uR', modFreq));
set(gca, 'fontsize', fSize);


subplot(122);
set(gca, 'xtick', snrs);
plot(-2*xShift+snrs, FLN_snr_eucl_neg_mean, 'b-o', 'linew', lw3, 'markersize', mrkSize);
plot(2*xShift+snrs, SSN_snr_eucl_neg_mean, 'r-o', 'linew', lw3, 'markersize', mrkSize);
grid on;
xlabel('SNR (dB)');
title(sprintf('ModFreq=%.0f Hz (corr(SN,S)/corr(SN,N)), -ve uR', modFreq));
set(gca, 'fontsize', fSize);


% LG2=legend('SSN', 'SSN-NF', 'SSN-dPrime', 'FLN', 'FLN-NF', 'FLN-dPrime', 'location', 'northwest');
% LG2.FontSize=leg_fSize;
% set(gca, 'fontsize', fSize);

set(hanOneModFreq_SNR, 'units', 'inches', 'position', [1 1 14 5]);
% fName_singleModFrequnCorr=sprintf('sEPSM_uR_f0_eucl_dist_%s_tRes%.0fms', postfix_fName, TimeResolution*1e3);

if saveAllFigs
    saveas(hanOneModFreq_SNR, [saving_Dir fName_singleModFrequnCorr], 'tiff');
    saveas(hanOneModFreq_SNR, [Latex_Dir fName_singleModFrequnCorr], 'epsc');
end


% % % % % end