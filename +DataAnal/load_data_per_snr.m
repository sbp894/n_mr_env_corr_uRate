% Works for the templates SNRenv_0, SNRenv_m5 and SNRenv_m10.

function [spike_data,StimData,StimsFNames]=load_data_per_snr(ExpControlParams,DataDir,~)

%% To return: These are the things that will be used later.
% Sampling Frequency, presentation level, condition_table, spike data corresponding
% to each condition in the condition table, paramsIN.durA_msec=dur_sec*1000;
% paramsIN.durB_msec=dur_sec*1000;
% paramsIN.durC_msec=dur_sec*1000;
% paramsIN.CF_A_Hz=CF_A_kHz*1000;
% paramsIN.CF_B_Hz=CF_B_kHz*1000;
% paramsIN.CF_C_Hz=CF_C_kHz*1000;
% % Need to include CF_A, CF_B, CF_C for more generality
% paramsIN.MAXspikes=2500;
% paramsIN.SNR2use_dB=SNR2use_dB;
% %% *MH July 7 2015 - change to longer delays to allow for lower modulations
% paramsIN.MAXdelay_sec=1;  % sentence duration is 2.7 sec, so 1 sec delays are about 1/3 of stim duration, which is pushing what we can estimate

CurDir=pwd;
addpath(CurDir);

%% load data
if ~strcmp(DataDir(end), filesep)
   DataDir= strcat(DataDir, filesep); 
end
cd (DataDir);
allfiles=dir('*SNRenv*');
% NoiseOrder=ExpControlParams.NoiseOrder; % 0 and 1 resp

condStruct=struct([]);
spike_data=struct([]);
StimData={};
StimsFNames={};
Ncases=0;

allCalibfiles=dir('*calib*');
fprintf('Using %s as calibration file\n', allCalibfiles(end).name);
CalibData=load(allCalibfiles(end).name);
CalibData=CalibData.data.CalibData;
TFiltWidthTC=5;

% rng('default'); % Important to reproduce same deviations (cfDevn2separateFiles) so that update_progress function can work properly

for file_var=1:length(allfiles)
    SNRenvFname=allfiles(file_var).name;
    
    fprintf('Working on %s\n', SNRenvFname);
    PicData=load(SNRenvFname);
    PicData=PicData.data;
    
    if isfield(PicData, 'screening')
        if PicData.screening.refract_violate_percent < ExpControlParams.THRESHOLD_percent_less_than_refractory
           proceed_flag=1; 
        else % Bad/ noisy triggering.
            proceed_flag=0;
            fprintf('%s file is excluded. %.2f%% of the ISI was less than absolute refractory period.\n', SNRenvFname, PicData.screening.refract_violate_percent);
        end
    else
       proceed_flag=1;
       warning('ISI has not been checked against absolute refractory period for file \n \t ---->%s', SNRenvFname);
       PicData.screening.refract_violate_percent=nan;
    end
    
    if proceed_flag==1
        unitNUM=PicData.General.unit;
        %     paramsUnit.unitNUM=PicData.General.track;
        trackNUM=PicData.General.track;
        %     paramsUnit.trackNUM=PicData.General.track;
        
        SRfName=[DataDir  sprintf('*_u%1d_%02d_SR*',trackNUM,unitNUM)];
        RLFfName=[DataDir  sprintf('*_u%1d_%02d_RLV*',trackNUM,unitNUM)];
        TCfName=[DataDir  sprintf('*_u%1d_%02d_tc*',trackNUM,unitNUM)];
        
        
        %% General parameters
        stim_list=PicData.Stimuli.list';
        paramsAnal.trigger=PicData.General.trigger;
        
        %% Rest
        %     eval(sprintf('x=Unit_%d_%02d;',trackNUM,unitNUM));
        Unit_Data=load(sprintf('Unit_%d_%02d.mat',trackNUM,unitNUM));
        Unit_Data=Unit_Data.data;
        
        paramsAnal.dur_msec=PicData.Hardware.Trigger.StmOn;
        paramsAnal.dur_total=PicData.Hardware.Trigger.StmOn+PicData.Hardware.Trigger.StmOff;
        paramsAnal.MAXspikes=inf;
        
        
        %     error('change to structure');
        sigType_SNR_polarity_nType_Mat=nan(length(stim_list),4); % Why not nan?
        stimParams=repmat(struct('sigType','','noiseType','','snr',nan, 'polarity',nan), length(stim_list), 1);
        
        
        for list_var=1:length(stim_list)
            cur_stim_name=stim_list{list_var};
            
            if contains(cur_stim_name,'SSN')
                sigType_SNR_polarity_nType_Mat(list_var,4)=0;
                stimParams(list_var).noiseType='SSN';
            elseif contains(cur_stim_name,'FLN')
                sigType_SNR_polarity_nType_Mat(list_var,4)=1;
                stimParams(list_var).noiseType='FLN';
            end
            
            cur_stim_name=cur_stim_name(strfind(cur_stim_name,'Stim'):end);
            
            if contains(cur_stim_name,'_N_')
                sigType_SNR_polarity_nType_Mat(list_var,1)=-1;
                stimParams(list_var).sigType='N';
            elseif contains(cur_stim_name,'_SN_')
                sigType_SNR_polarity_nType_Mat(list_var,1)=0;
                stimParams(list_var).sigType='SN';
            elseif contains(cur_stim_name,'_S_')
                sigType_SNR_polarity_nType_Mat(list_var,1)=1;
                stimParams(list_var).sigType='S';
            end
            
            SNRval=sscanf(cur_stim_name,'Stim%f*dB_%*s_%*s.wav');
            if ~isempty(SNRval)
                sigType_SNR_polarity_nType_Mat(list_var,2)=SNRval;
                stimParams(list_var).snr=SNRval;
            end
            
            if strcmp(cur_stim_name(end-5),'P')
                sigType_SNR_polarity_nType_Mat(list_var,3)=1;
                stimParams(list_var).polarity=1;
            else
                sigType_SNR_polarity_nType_Mat(list_var,3)=-1;
                stimParams(list_var).polarity=-1;
            end
        end
        
        %% How many csCells to add (case 1: multiple SNRs for same noise type. Case-2: single SNR for two different noises)
        
        %     noiseTypes_for_current_unit=unique(sigType_SNR_polarity_nType_Mat(:,4));
        noiseTypes_for_current_unit=unique({stimParams.noiseType}');
        if length(noiseTypes_for_current_unit)~=1
            condStruct(end+length(noiseTypes_for_current_unit)).CF=nan; %#ok<*AGROW>
            if isfield(Unit_Data, 'BFmod')
                fiberBF = 1e3*Unit_Data.BFmod; %% Update
            elseif isfield(Unit_Data, 'BF')
                fiberBF = 1e3*Unit_Data.BF; %% Update
            end
            %         csCell(end-length(noiseTypes_for_current_unit)+1:end,1)=fiberBF+cfDevn2separateFiles;
            [condStruct(end-length(noiseTypes_for_current_unit)+1:end).CF]=deal(fiberBF);
        end
        
        %% Calculate spontaneous rate
        if ~isempty(dir(SRfName))
            temp=dir(SRfName);
            data=load([DataDir temp(1).name]);
            data=data.data;
            num_of_lines= data.Stimuli.fully_presented_lines - length(data.Stimuli.bad_lines);
            if isempty(data.spikes{1,1})
                paramsUnit.SR=0;
                fprintf('The spont rate should be 0 for %s\n',SRFname);
            else
                paramsUnit.SR=size(data.spikes{1,1},1)/num_of_lines/(data.Hardware.Trigger.StmOn + data.Hardware.Trigger.StmOff);
            end
            
        elseif ~isempty(dir(SNRenvFname))
            temp=dir(SNRenvFname);
            data=load([DataDir temp(1).name]);
            data=data.data;
            spkData=data.spikes{1,1};
            tIgnore=.05;
            tStimOffStart=(data.Hardware.Trigger.StmOn+tIgnore*1e3)/1e3;
            spkData=spkData(spkData(:,1)<=data.Stimuli.fully_presented_stimuli,2);
            num_of_lines= data.Stimuli.fully_presented_lines - length(data.Stimuli.bad_lines);
            
            paramsUnit.SR=sum(spkData>tStimOffStart)/num_of_lines/(data.Hardware.Trigger.StmOff/1e3-tIgnore);
            warning('%s not present, need to be updated to exclude badlines (SNRenv file=%s)', SRfName, SNRenvFname);
        else
            paramsUnit.SR=nan;
        end
        
        %% Calculate Rate-Level Function
        if ~isempty(dir(RLFfName))
            temp=dir(RLFfName);
            data=load([DataDir temp(1).name]);
            data=data.data;
            spkData=data.spikes{1,1};
            dBs2run=nan(data.Stimuli.fully_presented_lines,1);
            rates=nan(data.Stimuli.fully_presented_lines,1);
            
            if ~isempty(data.Stimuli.bad_lines)
               disp('Debugging-- found ya!'); 
            end
            %%
            PlotRLV=0;
            plotFittedRLV=0;
            %%
            
            for line_var=1:length(rates)
                dBs2run(line_var)=line_var; % actual value doesn't matter. We are interested in sat rate
                spkCurLine=spkData(spkData(:,1)==line_var,2);
                if isempty(spkCurLine)
                    rates(line_var)=0;
                else
                    spkCurLine=spkCurLine(spkCurLine<data.Hardware.Trigger.StmOn/1e3);
                    rates(line_var)=length(spkCurLine)/(data.Hardware.Trigger.StmOn/1e3);
                end
            end
            [RLVparams,~,~]=NELfuns.fitRLfun(dBs2run,rates,PlotRLV,plotFittedRLV);
            %             CF_SR_SAT_Q10_3(unit_var)=RLVparams.R_Sat;
            paramsUnit.SatR=RLVparams.R_Sat;
            %         pause;
            %         clf;
            
        else
            %             CF_SR_SAT_Q10_3(unit_var)=nan;
            paramsUnit.SatR=nan;
        end
        
        %% Find Q10
        cd(CurDir);
        if isfield(Unit_Data, 'Q10_mod')
            paramsUnit.Q10=Unit_Data.Q10_mod;
        else
            if ~isempty(dir(TCfName))
                temp=dir(TCfName);
                TCdata=load([DataDir temp(1).name]);
                TCdata=TCdata.data.TcData;
                TCdata=TCdata(TCdata(:,1)~=0,:); % Get rid of all 0 freqs
                for i=1:size(TCdata,1)
                    TCdata(i,3)=NELfuns.CalibInterp(TCdata(i,1),CalibData)-TCdata(i,2);
                end
                TCdata(:,4)=NELfuns.trifilt(TCdata(:,3)',TFiltWidthTC)';
                [tempQ10,~,~,~] = NELfuns.findQ10(TCdata(:,1),TCdata(:,4),fiberBF/1e3);
                %             CF_SR_SAT_Q10_4(unit_var)=tempQ10;
                paramsUnit.Q10=tempQ10;
                
            else
                %             CF_SR_SAT_Q10_4(unit_var)=nan;
                paramsUnit.Q10=nan;
            end
        end
        cd (DataDir);
        
        %%
        
        
        %%
        
        %     paramsIN.level=CalibInterp(csCell(1,1)/1e3,CalibData)-PicData.Stimuli.attens;
        %     %%% If we are interested in energy at CF, else the value below
        %     should be alright.
        paramsAnal.level=65; % Hardcoded for NH. Also, also hardcoded in parse_saved_data_for_all_conds_SNRenv  %120-PicData.Stimuli.attens;
        
        
        %%
        for noise_var=1:length(noiseTypes_for_current_unit)
            Ncases=Ncases+1;
            %         curSNR=unique(sigType_SNR_polarity_nType_Mat(~isnan(sigType_SNR_polarity_nType_Mat(:,2)),2));
            curSNR=unique([stimParams.snr]);
            curSNR(isnan(curSNR))=[];
            
            cur_spike_data_ind=length(spike_data)+1;
            
            [spike_data(cur_spike_data_ind).SpikeTrains,StimData{cur_spike_data_ind},StimsFNames{cur_spike_data_ind}]=DataAnal.get_SpikeTrains_per_snr(PicData,stim_list,sigType_SNR_polarity_nType_Mat,stimParams,noiseTypes_for_current_unit(noise_var),curSNR);
            spike_data(cur_spike_data_ind).SNR=curSNR;
            spike_data(cur_spike_data_ind).noise=noiseTypes_for_current_unit{noise_var};
            spike_data(cur_spike_data_ind).track=trackNUM;
            spike_data(cur_spike_data_ind).unit=unitNUM;
            spike_data(cur_spike_data_ind).Q10=paramsUnit.Q10;
            spike_data(cur_spike_data_ind).SatR =paramsUnit.SatR;
            spike_data(cur_spike_data_ind).SR=paramsUnit.SR;
            spike_data(cur_spike_data_ind).CF_Hz=condStruct(Ncases).CF;
            spike_data(cur_spike_data_ind).StimsFNames=StimsFNames{end}; % StimsFNames is unnecessary.
            spike_data(cur_spike_data_ind).percent_less_than_refractory=PicData.screening.refract_violate_percent;
            
            paramsAnal.CF_Hz=spike_data(cur_spike_data_ind).CF_Hz;
            paramsAnal.Fs=PicData.Stimuli.updateRate_Hz;
            paramsAnal.SNR2use_dB=curSNR;
            paramsAnal.MAXdelay_sec=ExpControlParams.MAXdelay_fractionStimDur*paramsAnal.dur_msec/1e3;
            paramsAnal.SCC_onsetIGNORE_sec=ExpControlParams.SCC_onsetIGNORE_sec;
            paramsAnal.DELAYbinwidth_sec=ExpControlParams.DELAYbinwidth_sec;
            paramsAnal.ModSumExponent=ExpControlParams.ModSumExponent;
            fTypeNames={'LSR', 'MSR', 'HSR'};
            srTypeInd=([spike_data(cur_spike_data_ind).SR<.5 (spike_data(cur_spike_data_ind).SR>.5 & spike_data(cur_spike_data_ind).SR<18) spike_data(cur_spike_data_ind).SR>18]);
            paramsAnal.FiberType=fTypeNames{srTypeInd};
            paramsAnal.nReps=mean(mean(cellfun(@(x) numel(x),spike_data(cur_spike_data_ind).SpikeTrains)));
            paramsAnal.UseSlepian=ExpControlParams.UseSlepian;
            paramsAnal.BOOTSTRAP_fractiondata=ExpControlParams.BOOTSTRAP_fractiondata;
            
            
            spike_data(cur_spike_data_ind).paramsAnal=paramsAnal;
            spike_data(cur_spike_data_ind).SPL=paramsAnal.level;
            spike_data(cur_spike_data_ind).ChinID=cell2mat(cellfun(@(x) sscanf(char(x{1}), '-Q%d*'), regexp(DataDir,'(-Q\d+_)','tokens'), 'UniformOutput', 0));

            %         paramsIN.plt=1;
            %         paramsAnal.fTypes={'LSR','MSR','HSR'};
            %         spike_data(cur_spike_data_ind).FiberType=paramsAnal.fTypes{[paramsUnit.SR<.5 (.5<paramsUnit.SR & paramsUnit.SR<18) paramsUnit.SR>18]};
            %         spike_data(cur_spike_data_ind).CF=condStruct(Ncases).CF;
            %         spike_data(cur_spike_data_ind).paramsUnit=paramsUnit;
            
            
        end
    end
end

%%
cd (CurDir);
%csCell(:,3)=csCell(:,2);
[condStruct.sent]=deal(1);

for i=1:length(condStruct)
    condStruct(i).track=spike_data(i).track;
    condStruct(i).unit=spike_data(i).unit;
    condStruct(i).SNR=spike_data(i).SNR;
    condStruct(i).noise=spike_data(i).noise;
    condStruct(i).level=spike_data(i).paramsAnal.level;
end

% save([resultsDir 'conditions.mat'],'condStruct'); 
% earlier: function [spike_data,StimData,StimsFNames]=load_data_per_snr(ExpControlParams,DataDir,resultsDir)
% Now removed resultsDir. 
% Not sure if there will be any problem in running n_mr_sEPSM