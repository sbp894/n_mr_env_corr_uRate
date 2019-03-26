% Look at modulation filtered meanrates.
clear;
clc;

% chinIDs=[321 322 325 338 341 343 346 347 354 355];
% hanOneModFreqCorr_S_SN= 31;
% hanOneModFrequnCorr_SN_N=32;

ax= nan(4,1);
for nh=0:1
    if nh
        chinIDs=[321 322 325 338 341 343 346 347 354 355];
        hanOneModFreqCorr_S_SN=31;
        hanOneModFrequnCorr_SN_N=32;
        filePostFix= 'NH';
    else
        chinIDs= [361 362 358 360];
        hanOneModFreqCorr_S_SN=1031;
        hanOneModFrequnCorr_SN_N=1032;
        filePostFix= 'PTS';
    end
    
    % function plot_mean_rates_per_chin(chinIDs)
    hanModFilter=174;
    hanAllModFreqCorr=19;
    % fSize=14;
    leg_fSize=7;
    mrkSize=12;
    
    %% Important params
    saveAllFigs=0;
    N_half_bp=4;
    modFreq=64;
    fs=1e3;
    force_redo= 0;
    
    %%
    loading_Dir='/media/parida/DATAPART1/Matlab/SNRenv/n_mr_sEPSM/OUTPUT/DataAnal/';
    saving_Dir='/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/mr_corr_OUTPUT/';
    % Latex_Dir='/home/parida/Dropbox/Study Stuff/Presentations/TorstenPurdueVisit/Figures/';
    
    outFigDir.pdf=sprintf('/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/outFig/uR_info_modFreq_%d/pdf/', modFreq);
    if ~isdir(outFigDir.pdf)
        mkdir(outFigDir.pdf);
    end
    outFigDir.png=sprintf('/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/outFig/uR_info_modFreq_%d/png/', modFreq);
    if ~isdir(outFigDir.png)
        mkdir(outFigDir.png);
    end
    
    if ~isdir(saving_Dir)
        mkdir(saving_Dir);
    end
    
    meanrate_binwidth=.2e-3; % good for mod freq upto 1/2/binWidth
    
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
            
            if isfield(spike_data, 'thresh_dB') % slightly newer version of spike_data created using mr_sEPSM do not have thresh_dB
                spike_data=rmfield(spike_data, 'thresh_dB');
            end
            
            all_ChinSpikeData=[all_ChinSpikeData, spike_data];
            chin_snr_track_unit_mat=[chin_snr_track_unit_mat; [repmat(curChinID, length(spike_data),1), [spike_data.SNR]', [spike_data.track]', [spike_data.unit]']];
        end
    end
    
    unique_chin_snr_track_unit_mat=unique(chin_snr_track_unit_mat, 'rows');
    mr_corr_Data=repmat(struct('uRatePos', nan, 'uRateNeg', nan, 'CF_Hz', nan, 'SR', nan, 'SNR', nan, 'chinID', nan, 'track', nan, 'unit', nan), size(unique_chin_snr_track_unit_mat, 1), 1);
    
    %% Main parfor loop
    fName2Save=sprintf('%suRate_uR_info_modFreq_%d_%s.mat', saving_Dir, modFreq, filePostFix);
    
    if ~exist(fName2Save, 'file') || force_redo
        parfor plotVar= 1:size(unique_chin_snr_track_unit_mat,1)
            cur_inds=find(sum(repmat(unique_chin_snr_track_unit_mat(plotVar,:), size(chin_snr_track_unit_mat,1), 1)==chin_snr_track_unit_mat,2)==size(chin_snr_track_unit_mat,2));
            
            for indVar=1:length(cur_inds)
                
                curMetaData=all_ChinSpikeData(cur_inds(indVar));
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
                
                
                [mr_corr_Data(plotVar).uRatePos, mr_corr_Data(plotVar).uRateNeg]=mod_filter(S_rate_plus, S_rate_minus, 1/meanrate_binwidth, modFreq, N_half_bp, fs);
                
                mr_corr_Data(plotVar).CF_Hz=curMetaData.CF_Hz;
                mr_corr_Data(plotVar).SR=curMetaData.SR;
                mr_corr_Data(plotVar).SNR=curMetaData.SNR;
                mr_corr_Data(plotVar).chinID=curMetaData.chinID;
                mr_corr_Data(plotVar).track=curMetaData.track;
                mr_corr_Data(plotVar).unit=curMetaData.unit;
            end
        end
        save(fName2Save, 'mr_corr_Data', 'modFreq', 'fs');
    else
        temp= load(fName2Save);
        mr_corr_Data= temp.mr_corr_Data;
        modFreq=temp.modFreq;
        fs=temp.fs;
        clear temp;
    end
    
    [freq, sortInd]= unique([mr_corr_Data.CF_Hz]);
    time= (1:length(mr_corr_Data(1).uRatePos))/fs;
    % [freq, sortInd]= sort([mr_corr_Data.CF_Hz]);
    
    freqTicks= [.2 1:10]*1e3;
    freqInds= nan(size(freqTicks));
    
    for freqVar=1:length(freqTicks)
        freqInds(freqVar)= dsearchn(freq', freqTicks(freqVar));
    end
    
    [freqInds, upd_inds]= unique(freqInds);
    freqTicks= freqTicks(upd_inds);
    freqTickLabels= cellfun(@(x) num2str(x), num2cell(freqTicks), 'uniformoutput', false);
    
    tSize= length(time);
    fSize= length(freq);
    
    uRatePos= reshape([mr_corr_Data(sortInd).uRatePos], tSize, fSize)';
    uRatePos= uRatePos./repmat(rms(uRatePos,2), 1, size(uRatePos,2));
    figure(nh+1);
    clf;
    
    [stim,fs]= audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_08_15-Q361_AN_PTS/Signals/MH/SNRenv/SNR_0/FLN_Stim_S_P.wav');
    stim_time= (1:length(stim))/fs;
    
    ax(2*nh+1)=subplot(5,1,5);
    plot(stim_time ,stim);
    xlabel('time (sec)');
    ylabel('Stim Amp');
    
    ax(2*nh+2)=subplot(5,1,1:4);
    surf(time, 1:length(freq), uRatePos);
    view(2);
    grid on;
    xlabel('time (sec)');
    ylabel('CF (Hz)');
    zlabel('Filter output');
    title(sprintf('%s, mod filter CF = %.0f Hz', filePostFix, modFreq));
    set(gca, 'YTick', freqInds, 'YTickLabel', freqTickLabels);
    
    
    if saveAllFigs
        saveas(gcf, [fileparts(outFigDir.pdf(1:end-1)) filesep 'fig_for' filePostFix]);
    end
end
linkaxes(ax, 'x');