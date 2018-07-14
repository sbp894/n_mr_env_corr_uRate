function [corrValsModFreq, corrValsModFreqNF]=multires_modulation_hilb_tfs_f0(S_rate_plus_Org, S_rate_minus_Org, N_rate_plus_Org, N_rate_minus_Org, SN_rate_plus_Org, SN_rate_minus_Org, fsOrg, modFreq)

debugModeOn=0;
fs=5e3;

corrValsModFreq=nan(length(modFreq),1);
corrValsModFreqNF=nan(length(modFreq),1);

S_rate_plus=resample(S_rate_plus_Org, 1, fsOrg/fs);
S_rate_minus=resample(S_rate_minus_Org, 1, fsOrg/fs);
S_rate_tfs=(S_rate_plus-S_rate_minus)/2;
S_rate_tfs=abs(hilbert(S_rate_tfs));

N_rate_plus=resample(N_rate_plus_Org, 1, fsOrg/fs);
N_rate_minus=resample(N_rate_minus_Org, 1, fsOrg/fs);
N_rate_tfs=(N_rate_plus-N_rate_minus)/2;
N_rate_tfs=abs(hilbert(N_rate_tfs));


SN_rate_plus=resample(SN_rate_plus_Org, 1, fsOrg/fs);
SN_rate_minus=resample(SN_rate_minus_Org, 1, fsOrg/fs);
SN_rate_tfs=(SN_rate_plus-SN_rate_minus)/2;
SN_rate_tfs=abs(hilbert(SN_rate_tfs));


stimDur=length(S_rate_plus)/fs;
stimTime=(1:length(S_rate_plus))/fs;

for modFreqVar=1:length(modFreq)
    curModFreq=modFreq(modFreqVar);
    %     curTimeResolution=1/curModFreq;
    warning('Change');
    curTimeResolution=20e-3;
    N_bp=8;
    HalfPowerFrequency1=modFreq/sqrt(2);
    HalfPowerFrequency2=modFreq*sqrt(2);
    curFilt= designfilt('bandpassiir','FilterOrder',N_bp, ...get(p,props)
        'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
        'SampleRate',fs);
    
    if debugModeOn
        fvtool(curFilt);
        xlim([curModFreq/2 2*curModFreq]/1e3);
        ylim([-200 0])
    end
    
    
    if rem(stimDur, curTimeResolution) < .5*curTimeResolution
        curCorrVals=nan(floor(stimDur/curTimeResolution), 1);
        curCorrValsNF=nan(floor(stimDur/curTimeResolution), 1);
    else
        curCorrVals=nan(ceil(stimDur/curTimeResolution), 1);
        curCorrValsNF=nan(ceil(stimDur/curTimeResolution), 1);
    end
    
    S_rate_tfs_filt=filtfilt(curFilt, S_rate_tfs);
    SN_rate_tfs_filt=filtfilt(curFilt, SN_rate_tfs);
    N_rate_tfs_filt=filtfilt(curFilt, N_rate_tfs);
    
    
    
    for windowVar=1:length(curCorrVals)
        tStart=(windowVar-1)*curTimeResolution;
        indStart=max(1, round(tStart*fs));
        tEnd=windowVar*curTimeResolution;
        indEnd=min(length(S_rate_plus), round(tEnd*fs));
        validINDs=indStart:indEnd;
        
        curCorrVals(windowVar)=corr2(S_rate_tfs_filt(validINDs), SN_rate_tfs_filt(validINDs));
        curCorrValsNF(windowVar)=corr2(S_rate_tfs_filt(validINDs), N_rate_tfs_filt(validINDs));
        
        if debugModeOn
            figure(1012);
            %             if windowVar==1
            clf;
            %             end
            hold on;
            lw=.8;
            Ashift=max(S_rate_plus);
            plot(stimTime, Ashift+S_rate_plus, 'b-.', 'linew', lw);
            plot(stimTime, Ashift+S_rate_minus, 'm:', 'linew', lw);
            plot(stimTime, SN_rate_plus, 'b-.', 'linew', lw);
            plot(stimTime, SN_rate_minus, 'm:', 'linew', lw);
            plot(stimTime, -Ashift+N_rate_plus, 'b-.', 'linew', lw);
            plot(stimTime, -Ashift+N_rate_minus, 'm:', 'linew', lw);
            
            yyaxis right;
            lw=1.5;
            Ashift=max(abs(S_rate_tfs_filt));
            plot(stimTime, Ashift+S_rate_tfs_filt, 'k-', 'linew', lw);
            plot(stimTime, SN_rate_tfs_filt, 'k-', 'linew', lw);
            plot(stimTime, -Ashift+N_rate_tfs_filt, 'k--', 'linew', lw);
        end
    end
    curCorrVals(curCorrVals<0)=0;
    curCorrValsNF(curCorrValsNF<0)=0;
    %     corrValsModFreq(modFreqVar)=sqrt(nansum(curCorrVals.^2));
    %     corrValsModFreqNF(modFreqVar)=sqrt(nansum(curCorrValsNF.^2));
    warning('chjagnekam');
    corrValsModFreq(modFreqVar)=(nansum(curCorrVals.^1));
    corrValsModFreqNF(modFreqVar)=(nansum(curCorrValsNF.^1));
end