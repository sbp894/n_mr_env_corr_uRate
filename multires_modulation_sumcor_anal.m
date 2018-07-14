function [corrValsModFreq, corrValsModFreqNF]=multires_modulation_sumcor_anal(S_rate_plus_Org, S_rate_minus_Org, N_rate_plus_Org, N_rate_minus_Org, SN_rate_plus_Org, SN_rate_minus_Org, fsOrg, AllmodFreqs)

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
    
    for windowVar=1:length(curCorrVals)
        tStart=(windowVar-1)*curTimeResolution;
        indStart=max(1, round(tStart*fs));
        tEnd=windowVar*curTimeResolution;
        indEnd=min(length(S_rate_plus), round(tEnd*fs));
        validINDs=indStart:indEnd;
        
        s_acf=(xcorr(S_rate_plus(validINDs))+xcorr(S_rate_minus(validINDs)))/2;
        s_xcf=xcorr(S_rate_plus(validINDs), S_rate_minus(validINDs));
        s_env=(s_acf+s_xcf)/2;
        s_env_filt=filter(curFilt, s_env);
        
        n_acf=(xcorr(N_rate_plus(validINDs))+xcorr(N_rate_minus(validINDs)))/2;
        n_xcf=xcorr(N_rate_plus(validINDs), N_rate_minus(validINDs));
        n_env=(n_acf+n_xcf)/2;
        n_env_filt=filter(curFilt, n_env);
        
        sn_acf=(xcorr(SN_rate_plus(validINDs))+xcorr(SN_rate_minus(validINDs)))/2;
        [sn_xcf, delay_ind]=xcorr(SN_rate_plus(validINDs), SN_rate_minus(validINDs));
        sn_env=(sn_acf+sn_xcf)/2;
        sn_env_filt=filter(curFilt, sn_env);

        curCorrVals(windowVar)=corr2(s_env_filt, sn_env_filt);
        curCorrValsNF(windowVar)=corr2(s_env_filt, n_env_filt);
        
        if debugModeOn
            figure(1012);
            if windowVar==1
                clf;
            end
            hold on;
            delay=delay_ind/fs;
            plot(2*curTimeResolution*(windowVar-1)+delay, s_env, 'k');
            plot(2*curTimeResolution*(windowVar-1)+delay, s_env_filt, 'b');
            plot(2*curTimeResolution*(windowVar-1)+delay, -max(s_env_filt)+sn_env_filt, 'r');
            plot(2*curTimeResolution*(windowVar-1)+delay, -2*max(s_env_filt)+n_env_filt, 'g');
        end
    end
    curCorrVals(curCorrVals<0)=0;
    curCorrValsNF(curCorrValsNF<0)=0;
    corrValsModFreq(modFreqVar)=sqrt(nansum(curCorrVals.^2));
    corrValsModFreqNF(modFreqVar)=sqrt(nansum(curCorrValsNF.^2));
end