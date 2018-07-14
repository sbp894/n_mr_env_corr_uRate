function [corrValsModFreq, corrValsModFreqNF]=multires_modulation_XpAC_anal(S_rate_plus_Org, S_rate_minus_Org, N_rate_plus_Org, N_rate_minus_Org, SN_rate_plus_Org, SN_rate_minus_Org, fsOrg, AllmodFreqs)

debugModeOn=0;
fs=5e3;

corrValsModFreq=nan(length(AllmodFreqs),1);
corrValsModFreqNF=nan(length(AllmodFreqs),1);

S_rate_plus=resample(S_rate_plus_Org, 1, fsOrg/fs);
S_rate_minus=resample(S_rate_minus_Org, 1, fsOrg/fs);
N_rate_plus=resample(N_rate_plus_Org, 1, fsOrg/fs);
N_rate_minus=resample(N_rate_minus_Org, 1, fsOrg/fs);
SN_rate_plus=resample(SN_rate_plus_Org, 1, fsOrg/fs);
SN_rate_minus=resample(SN_rate_minus_Org, 1, fsOrg/fs);
stimDur=length(S_rate_plus)/fs;
stimTime=(1:length(S_rate_plus))/fs;

for modFreqVar=1:length(AllmodFreqs)
    curModFreq=AllmodFreqs(modFreqVar);
    curTimeResolution=1/curModFreq;
    if modFreqVar==1
        N_lp=8;
        Fstop_lp=sqrt(AllmodFreqs(1)*AllmodFreqs(2));
        curFilt= designfilt('lowpassiir','FilterOrder',N_lp, ...
            'PassbandFrequency',Fstop_lp,'PassbandRipple',0.2, ...
            'SampleRate',fs);
    else
        N_bp=8;
        HalfPowerFrequency1=sqrt(AllmodFreqs(modFreqVar)*AllmodFreqs(modFreqVar-1));
        if modFreqVar~=length(AllmodFreqs)
            HalfPowerFrequency2=sqrt(AllmodFreqs(modFreqVar)*AllmodFreqs(modFreqVar+1));
        else
            HalfPowerFrequency2=sqrt(2)*AllmodFreqs(modFreqVar);
        end
        
        curFilt= designfilt('bandpassiir','FilterOrder',N_bp, ...get(p,props)
            'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
            'SampleRate',fs);
    end
    
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
    
    S_rate_plus_filt=filtfilt(curFilt, S_rate_plus);
    S_rate_minus_filt=filtfilt(curFilt, S_rate_minus);
    N_rate_plus_filt=filtfilt(curFilt, N_rate_plus);
    N_rate_minus_filt=filtfilt(curFilt, N_rate_minus);
    SN_rate_plus_filt=filtfilt(curFilt, SN_rate_plus);
    SN_rate_minus_filt=filtfilt(curFilt, SN_rate_minus);
    
    
    for windowVar=1:length(curCorrVals)
        tStart=(windowVar-1)*curTimeResolution;
        indStart=max(1, round(tStart*fs));
        tEnd=windowVar*curTimeResolution;
        indEnd=min(length(S_rate_plus), round(tEnd*fs));
        validINDs=indStart:indEnd;
        
        s_sn_1=corr2(S_rate_plus_filt(validINDs), SN_rate_minus_filt(validINDs));
        s_sn_2=corr2(S_rate_minus_filt(validINDs), SN_rate_plus_filt(validINDs));
        curCorrVals(windowVar)=(s_sn_1+s_sn_2)/2;
        
        s_n_1=corr2(S_rate_plus_filt(validINDs), N_rate_minus_filt(validINDs));
        s_n_2=corr2(S_rate_minus_filt(validINDs), N_rate_plus_filt(validINDs));
        curCorrValsNF(windowVar)=(s_n_1+s_n_2)/2;
        
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
            Ashift=max(abs(S_rate_plus_filt));
            plot(stimTime, Ashift+S_rate_plus_filt, 'k-', 'linew', lw);
            plot(stimTime, Ashift+S_rate_minus_filt, 'g--', 'linew', lw);
            plot(stimTime, SN_rate_plus_filt, 'k-', 'linew', lw);
            plot(stimTime, SN_rate_minus_filt, 'g--', 'linew', lw);
            plot(stimTime, -Ashift+N_rate_plus_filt, 'k-', 'linew', lw);
            plot(stimTime, -Ashift+N_rate_minus_filt, 'g--', 'linew', lw);
        end
    end
    curCorrVals(curCorrVals<0)=0;
    curCorrValsNF(curCorrValsNF<0)=0;
    corrValsModFreq(modFreqVar)=sqrt(nansum(curCorrVals.^2));
    corrValsModFreqNF(modFreqVar)=sqrt(nansum(curCorrValsNF.^2));
end