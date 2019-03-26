clear;
clc;


% create a directory to save spike_stim_data

dirStruct.WorkingDir='/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/mr_corr_OUTPUT/';
dirStruct.MATDataDir='/media/parida/DATAPART1/Matlab/ExpData/MatData/';
dirStruct.CodesDir= '/media/parida/DATAPART1/Matlab/SNRenv/n_mr_sEPSM/Codes/';
dirStruct.OutDir=[dirStruct.WorkingDir 'InData/DanishData/'];
if ~isdir(dirStruct.OutDir)
    mkdir(dirStruct.OutDir);
end

addpath(dirStruct.CodesDir);
chinIDs=[321 322 325 338 341 343 346 347 354 355 358 360 361 362];
% chinIDs=[355 361 362];

parfor chinVar=1:length(chinIDs)
    curChinID= chinIDs(chinVar);
    fName2Save=sprintf('%sQ%d_spikestimulusData', dirStruct.OutDir, curChinID);
    loopChinID=chinIDs(chinVar);
    loopDataDir=dir([dirStruct.MATDataDir '*Q' num2str(loopChinID) '*AN*']);
    %     if length(loopDataDir)>1
    %         loopDataDir= loopDataDir(contains({loopDataDir.name}', '_AN_'));
    %     end
    [spike_data, ~, ~]= DataAnal.load_data_per_snr(DataAnal.get_ExpControlParams, [dirStruct.MATDataDir loopDataDir.name]);
    
    dummy_save([fName2Save '.mat'], spike_data);
end

rmpath(dirStruct.CodesDir);


function dummy_save(fName, spike_data) %#ok<INUSD>
save(fName, 'spike_data');
end