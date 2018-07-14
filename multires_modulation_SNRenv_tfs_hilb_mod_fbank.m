function [SNRenv_modFreq]=...
    multires_modulation_SNRenv_tfs_hilb_mod_fbank(S_rate_plus_Org, S_rate_minus_Org, N_rate_plus_Org, N_rate_minus_Org, SN_rate_plus_Org, SN_rate_minus_Org, fsOrg, modFreq)

debugModeOn1=0;
debugModeOn2=0;
fs=2e3;

SNRenv_modFreq=nan(length(modFreq),1); 

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
    curTimeResolution=1/curModFreq;
    if curModFreq==1
        N_lp=4;
        Fstop_lp=curModFreq*sqrt(2);
        curFilt= designfilt('lowpassiir','FilterOrder',N_lp, ...
            'PassbandFrequency',Fstop_lp,'PassbandRipple',0.2, ...
            'SampleRate',fs);
    else
        N_bp=8;
        HalfPowerFrequency1=curModFreq/sqrt(2);
        HalfPowerFrequency2=curModFreq*sqrt(2);
        curFilt= designfilt('bandpassiir','FilterOrder',N_bp, ...get(p,props)
            'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
            'SampleRate',fs);
    end
    
    if debugModeOn1
        fvtool(curFilt);
        xlim([curModFreq/2 2*curModFreq]/1e3);
        ylim([-200 0])
    end
    
    
    if rem(stimDur, curTimeResolution) < .5*curTimeResolution
        cur_SNRenv_vals=nan(floor(stimDur/curTimeResolution), 1);
    else
        cur_SNRenv_vals=nan(ceil(stimDur/curTimeResolution), 1);
    end
    
    S_rate_tfs_filt=filtfilt(curFilt, S_rate_tfs);
    SN_rate_tfs_filt=filtfilt(curFilt, SN_rate_tfs);
    N_rate_tfs_filt=filtfilt(curFilt, N_rate_tfs);
    
    
    for windowVar=1:length(cur_SNRenv_vals)
        tStart=(windowVar-1)*curTimeResolution;
        indStart=max(1, round(tStart*fs));
        tEnd=windowVar*curTimeResolution;
        indEnd=min(length(S_rate_plus), round(tEnd*fs));
        validINDs=indStart:indEnd;
        
        pSN=var(SN_rate_tfs_filt(validINDs));
        pN=var(N_rate_tfs_filt(validINDs));
        pSest=max(0, pSN-pN);
                
        cur_SNRenv_vals(windowVar)=pSest/pN;
        
        if debugModeOn2
            figure(1012);
            %             if windowVar==1
            clf;
            %             end
            hold on;
            lw=.8;
            Ashift=max(S_rate_plus);
            plot(stimTime, Ashift+S_rate_plus, 'b', 'linew', lw);
            plot(stimTime, Ashift+S_rate_minus, 'm', 'linew', lw);
            plot(stimTime, SN_rate_plus, 'b', 'linew', lw);
            plot(stimTime, SN_rate_minus, 'm', 'linew', lw);
            plot(stimTime, -Ashift+N_rate_plus, 'b', 'linew', lw);
            plot(stimTime, -Ashift+N_rate_minus, 'm', 'linew', lw);
            
            yyaxis right;
            lw=1.5;
            Ashift=max(abs(S_rate_tfs_filt));
            plot(stimTime, Ashift+S_rate_tfs_filt-mean(S_rate_tfs_filt), 'k-', 'linew', lw);
            plot(stimTime, SN_rate_tfs_filt-mean(SN_rate_tfs_filt), 'k-', 'linew', lw);
            plot(stimTime, -Ashift+N_rate_tfs_filt-mean(N_rate_tfs_filt), 'k-', 'linew', lw);
        end
    end
    
    cur_SNRenv_vals(cur_SNRenv_vals<0)=0;

    SNRenv_modFreq(modFreqVar)=nanmean(cur_SNRenv_vals);
end