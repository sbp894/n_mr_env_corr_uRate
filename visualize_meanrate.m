% Look at modulation filtered meanrates.

chinIDs=[321 322 325 338 341 343];
% function plot_mean_rates_per_chin(chinIDs)
hanModFilter=174;

loading_Dir='/media/parida/DATAPART1/Matlab/SNRenv/n_sEPSM/OUTPUT/DataAnal/';
saving_Dir='/media/parida/DATAPART1/Matlab/SNRenv/mr_n_SNRrate/OUT_figures/meanrates/';
LatexDir='/home/parida/Dropbox/Study Stuff/Presentations/TorstenPurdueVisit/Figures/';

saveFigs=0;

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
mr_corr_Data=repmat(struct('SSNcorr', nan, 'FLNcorr', nan, 'CF_Hz', nan, 'SR', nan), size(unique_chin_snr_track_unit_mat, 1), 1);

for plotVar=5%1:length(unique_chin_snr_track_unit_mat)
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
        
        AllmodFreqs=2.^(0:8);
        %         multires_modulation_corr_anal(S_rate_env, SN_rate_env, N_rate_env, 1/meanrate_binwidth, AllmodFreqs);
        clf;
        plot_filtered_modulation(S_rate_env, 1/meanrate_binwidth)
    end
    if plotVar==5 % good example plot for Torsten's visit, showing pitch related modulation as well as syllable-rate modulation
        set(gcf, 'units', 'inches', 'position', [1 1 10 10]);
        set(gca, 'fontsize', 12);
        ylabel('arb. units');
        title(sprintf('uRate S-alone/ CF=%.0f kHz, SR=%.0f', curMetaData.CF_Hz/1e3, curMetaData.SR));
        fName2Save=strrep(sprintf('%s_modFilteredData', fName), '.', '_');
        ylim([-18.5 2]);
        xlim([-.01 1.31])
        if saveFigs
            saveas(gcf, [LatexDir fName2Save], 'epsc');
        end
    end
end