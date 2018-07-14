%%% Should be added to path
% Should work on any machine. 

function [SpikeTrains,StimsData,StimsFNames]=get_SpikeTrains_per_snr(PicData,stim_list,~,stimParams,curNoiseType,~)
% sigType_SNR_polarity_nType_Mat(isnan(sigType_SNR_polarity_nType_Mat))=inf;
% curSNR=unique(sigType_SNR_polarity_nType_Mat(~isinf(sigType_SNR_polarity_nType_Mat(:,2)),2));

%% create sigType_SNR_polarity_nType_Mat
% sigType_SNR_polarity_nType_Mat=inf(length(stimParams),4);

% sigType_SNR_polarity_nType_Mat(:,1)={stimParams.sigType}'


%%
SpikeTrains=cell(3,2);
% StimsData=cell(3,2);
all_spikes_times=PicData.spikes{1};
all_files_played=PicData.Line.file;
all_files_played=all_files_played(1:PicData.Stimuli.fully_presented_lines);

% Speech Positive 
% stim_name=stim_list(ismember(sigType_SNR_polarity_nType_Mat(:,1:3),[1 inf 1],'rows'));
SPname=stim_list{strcmp({stimParams.sigType}', 'S') & ([stimParams.polarity]==1)'};
% SpikeTrainsS_plus={};

% [~,DataFolderName,~]=fileparts(pwd);
% UserName=DataFolderName(1:2);
UserName='MH';

[path_stim,~,~]=fileparts(strrep(SPname, '\', filesep));
stim_name2load=strtrim(strrep(SPname(length(fileparts(path_stim))+2:end), '\', filesep));
StimsS_plus_fName=[pwd filesep 'Signals' filesep UserName filesep 'SNRenv' filesep stim_name2load];
SpeechStimsS_plus=audioread(StimsS_plus_fName);

fileInds=find(strcmp(all_files_played, SPname)==1);
SpikeTrainsS_plus=cell(length(fileInds),1);


for file_var=1:length(fileInds)
        SpikeTrainsS_plus{file_var}=sort(all_spikes_times((all_spikes_times(:,1)==fileInds(file_var)),2));
end

% Speech Negative 
SNname=stim_list{strcmp({stimParams.sigType}', 'S') & ([stimParams.polarity]==-1)'};
% SpikeTrainsS_minus={};

[path_stim,~,~]=fileparts(strrep(SNname, '\', filesep));
stim_name2load=strtrim(strrep(SNname(length(fileparts(path_stim))+2:end), '\', filesep));
StimsS_minus_fName=[pwd filesep 'Signals' filesep UserName filesep 'SNRenv' filesep stim_name2load];
SpeechStimsS_minus=audioread(StimsS_minus_fName);

fileInds=find(strcmp(all_files_played, SNname)==1);
SpikeTrainsS_minus=cell(length(fileInds),1);

for file_var=1:length(fileInds)
        SpikeTrainsS_minus{file_var}=sort(all_spikes_times((all_spikes_times(:,1)==fileInds(file_var)),2));
end



% Noise Positive 
NPname=stim_list{strcmp({stimParams.sigType}', 'N') & ([stimParams.polarity]==1)' & strcmp({stimParams.noiseType}' ,curNoiseType)};
% stim_name=stim_list(ismember(sigType_SNR_polarity_nType_Mat,[-1 curSNR 1 curNoiseType],'rows'));
% SpikeTrainsN_plus={};

[path_stim,~,~]=fileparts(strrep(NPname, '\', filesep));
stim_name2load=strtrim(strrep(NPname(length(fileparts(path_stim))+2:end), '\', filesep));
StimsN_plus_fName=[pwd filesep 'Signals' filesep UserName filesep 'SNRenv' filesep stim_name2load];
SpeechStimsN_plus=audioread(StimsN_plus_fName);

fileInds=find(strcmp(all_files_played, NPname)==1);
SpikeTrainsN_plus=cell(length(fileInds),1);

for file_var=1:length(fileInds)
        SpikeTrainsN_plus{file_var}=sort(all_spikes_times((all_spikes_times(:,1)==fileInds(file_var)),2));
end

% Noise dB Negative 
NNname=stim_list{strcmp({stimParams.sigType}', 'N') & ([stimParams.polarity]==-1)' & strcmp({stimParams.noiseType}' ,curNoiseType)};
% stim_name=stim_list(ismember(sigType_SNR_polarity_nType_Mat,[-1 curSNR 1 curNoiseType],'rows'));
% SpikeTrainsN_plus={};

[path_stim,~,~]=fileparts(strrep(NNname, '\', filesep));
stim_name2load=strtrim(strrep(NNname(length(fileparts(path_stim))+2:end), '\', filesep));
StimsN_minus_fName=[pwd filesep 'Signals' filesep UserName filesep 'SNRenv' filesep stim_name2load];
SpeechStimsN_minus=audioread(StimsN_minus_fName);

fileInds=find(strcmp(all_files_played, NNname)==1);
SpikeTrainsN_minus=cell(length(fileInds),1);

for file_var=1:length(fileInds)
        SpikeTrainsN_minus{file_var}=sort(all_spikes_times((all_spikes_times(:,1)==fileInds(file_var)),2));
end

% Speech+Noise SNR dB Positive 
SNPname=stim_list{strcmp({stimParams.sigType}', 'SN') & ([stimParams.polarity]==1)' & strcmp({stimParams.noiseType}' ,curNoiseType)};
% stim_name=stim_list(ismember(sigType_SNR_polarity_nType_Mat,[-1 curSNR 1 curNoiseType],'rows'));
% SpikeTrainsN_plus={};

[path_stim,~,~]=fileparts(strrep(SNPname, '\', filesep));
stim_name2load=strtrim(strrep(SNPname(length(fileparts(path_stim))+2:end), '\', filesep));
StimsSN_plus_fName=[pwd filesep 'Signals' filesep UserName filesep 'SNRenv' filesep stim_name2load];
SpeechStimsSN_plus=audioread(StimsSN_plus_fName);

fileInds=find(strcmp(all_files_played, SNPname)==1);
SpikeTrainsSN_plus=cell(length(fileInds),1);

for file_var=1:length(fileInds)
        SpikeTrainsSN_plus{file_var}=sort(all_spikes_times((all_spikes_times(:,1)==fileInds(file_var)),2));
end


% stim_name=stim_list(ismember(sigType_SNR_polarity_nType_Mat,[0 curSNR 1 curNoiseType],'rows'));
% SpikeTrainsSN_plus={};
% 
% [path_stim,~,~]=fileparts(strrep(stim_name, '\', filesep));
% stim_name2load=strtrim(strrep(stim_name(length(fileparts(path_stim))+2:end), '\', filesep));
% StimsSN_plus_fName=[pwd filesep 'Signals' filesep UserName filesep 'SNRenv' filesep stim_name2load];
% SpeechStimsSN_plus=audioread(StimsSN_plus_fName);

% for file_var=1:length(all_files_played)
%     if strcmp(all_files_played{file_var},stim_name)
%         SpikeTrainsSN_plus{end+1}=(sort(all_spikes_times((all_spikes_times(:,1)==file_var),2)));
%     end
%     
% end

% Speech+Noise SNR dB Negative 
SNNname=stim_list{strcmp({stimParams.sigType}', 'SN') & ([stimParams.polarity]==1)' & strcmp({stimParams.noiseType}' ,curNoiseType)};
% stim_name=stim_list(ismember(sigType_SNR_polarity_nType_Mat,[-1 curSNR 1 curNoiseType],'rows'));
% SpikeTrainsN_plus={};

[path_stim,~,~]=fileparts(strrep(SNNname, '\', filesep));
stim_name2load=strtrim(strrep(SNNname(length(fileparts(path_stim))+2:end), '\', filesep));
StimsSN_minus_fName=[pwd filesep 'Signals' filesep UserName filesep 'SNRenv' filesep stim_name2load];
SpeechStimsSN_minus=audioread(StimsSN_minus_fName);

fileInds=find(strcmp(all_files_played, SNNname)==1);
SpikeTrainsSN_minus=cell(length(fileInds),1);

for file_var=1:length(fileInds)
        SpikeTrainsSN_minus{file_var}=sort(all_spikes_times((all_spikes_times(:,1)==fileInds(file_var)),2));
end

% stim_name=stim_list(ismember(sigType_SNR_polarity_nType_Mat,[0 curSNR -1 curNoiseType],'rows'));
% SpikeTrainsSN_minus={};
% 
% [path_stim,~,~]=fileparts(strrep(stim_name, '\', filesep));
% stim_name2load=strtrim(strrep(stim_name(length(fileparts(path_stim))+2:end), '\', filesep));
% StimsSN_minus_fName=[pwd filesep 'Signals' filesep UserName filesep 'SNRenv' filesep stim_name2load];
% SpeechStimsSN_minus=audioread(StimsSN_minus_fName);
% 
% for file_var=1:length(all_files_played)
%     if strcmp(all_files_played{file_var},stim_name)
%         SpikeTrainsSN_minus{end+1}=(sort(all_spikes_times((all_spikes_times(:,1)==file_var),2)));
%     end
%     
% end

% SpikeTrains=[SpikeTrainsS_plus,SpikeTrainsS_minus;SpikeTrainsN_plus,SpikeTrainsN_minus;SpikeTrainsSN_plus,SpikeTrainsSN_minus];

SpikeTrains{1,1}=SpikeTrainsS_plus;
SpikeTrains{1,2}=SpikeTrainsS_minus;
SpikeTrains{2,1}=SpikeTrainsN_plus;
SpikeTrains{2,2}=SpikeTrainsN_minus;
SpikeTrains{3,1}=SpikeTrainsSN_plus;
SpikeTrains{3,2}=SpikeTrainsSN_minus;

StimsData={{SpeechStimsS_plus},{SpeechStimsS_minus};{SpeechStimsN_plus},{SpeechStimsN_minus};{SpeechStimsSN_plus},{SpeechStimsSN_minus}};
StimsFNames={{StimsS_plus_fName},{StimsS_minus_fName};{StimsN_plus_fName},{StimsN_minus_fName};{StimsSN_plus_fName},{StimsSN_minus_fName}};