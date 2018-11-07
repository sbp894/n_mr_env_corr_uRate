function [corr_s_sn_ValsModFreq, uncorr_sn_n_ValsModFreq, corr_s_n_ValsModFreq]=...
    multires_modulation_tfs_hilb_f0(S_rate_plus_Org, S_rate_minus_Org, N_rate_plus_Org, N_rate_minus_Org, SN_rate_plus_Org, SN_rate_minus_Org, fsOrg, modFreq, outFigDir, fName)

if ~exist('outFigDir', 'var')
    outFigDir=[];
    fName=[];
end


debugModeOn1=0;
debugModeOn2=1;
fs=5e3;
fSize=14;
textSize=20;
SetOnsetZeroDuration=.01;
OnsetInds2SetZero=round(fs*SetOnsetZeroDuration);

corr_s_sn_ValsModFreq=nan(length(modFreq),1);
uncorr_sn_n_ValsModFreq=nan(length(modFreq),1);
corr_s_n_ValsModFreq=nan(length(modFreq),1);

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
    N_bp=8;
    HalfPowerFrequency1=modFreq/sqrt(2);
    HalfPowerFrequency2=modFreq*sqrt(2);
    curFilt= designfilt('bandpassiir','FilterOrder',N_bp, ... 
        'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
        'SampleRate',fs);
    
    if debugModeOn1
        fvtool(curFilt);
        xlim([curModFreq/2 2*curModFreq]/1e3);
        ylim([-200 0])
    end
    
    
    if rem(stimDur, curTimeResolution) < .5*curTimeResolution
        cur_s_sn_CorrVals=nan(floor(stimDur/curTimeResolution), 1);
        cur_sn_n_CorrVals=nan(floor(stimDur/curTimeResolution), 1);
        cur_s_n_CorrVals=nan(floor(stimDur/curTimeResolution), 1);
    else
        cur_s_sn_CorrVals=nan(ceil(stimDur/curTimeResolution), 1);
        cur_sn_n_CorrVals=nan(ceil(stimDur/curTimeResolution), 1);
        cur_s_n_CorrVals=nan(ceil(stimDur/curTimeResolution), 1);
    end
    
    S_rate_tfs_filt=filtfilt(curFilt, S_rate_tfs);
    S_rate_tfs_filt(1:OnsetInds2SetZero)=0;
    SN_rate_tfs_filt=filtfilt(curFilt, SN_rate_tfs);
    SN_rate_tfs_filt(1:OnsetInds2SetZero)=0;
    N_rate_tfs_filt=filtfilt(curFilt, N_rate_tfs);
    N_rate_tfs_filt(1:OnsetInds2SetZero)=0;
    
    
    for windowVar=1:length(cur_s_sn_CorrVals)
        tStart=(windowVar-1)*curTimeResolution;
        indStart=max(1, round(tStart*fs));
        tEnd=windowVar*curTimeResolution;
        indEnd=min(length(S_rate_plus), round(tEnd*fs));
        validINDs=indStart:indEnd;
        
        cur_s_sn_CorrVals(windowVar)=corr2(S_rate_tfs_filt(validINDs), SN_rate_tfs_filt(validINDs));
        cur_sn_n_CorrVals(windowVar)=corr2(SN_rate_tfs_filt(validINDs), N_rate_tfs_filt(validINDs));
        cur_s_n_CorrVals(windowVar)=corr2(S_rate_tfs_filt(validINDs), N_rate_tfs_filt(validINDs));
    end
    
    cur_s_sn_CorrVals(cur_s_sn_CorrVals<0)=0;
    cur_sn_n_CorrVals(cur_sn_n_CorrVals<0)=0;
    cur_sn_n_unCorrVals=1-cur_sn_n_CorrVals;
    cur_s_n_CorrVals(cur_s_n_CorrVals<0)=0;
    
    corr_s_sn_ValsModFreq(modFreqVar)=sqrt(nansum(cur_s_sn_CorrVals.^2));
    uncorr_sn_n_ValsModFreq(modFreqVar)=sqrt(nansum(cur_sn_n_unCorrVals.^2));
    corr_s_n_ValsModFreq(modFreqVar)=sqrt(nansum(cur_s_n_CorrVals.^2));
    
    if ~isempty(outFigDir) && debugModeOn2
        %             figure(1012);
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
        Ashift=1.75*max(abs([S_rate_tfs_filt SN_rate_tfs_filt N_rate_tfs_filt]));
        plot(stimTime, Ashift+S_rate_tfs_filt, 'k-', 'linew', lw);
        plot(stimTime, SN_rate_tfs_filt, 'k-', 'linew', lw);
        plot(stimTime, -Ashift+N_rate_tfs_filt, 'k-', 'linew', lw);
        
        yyaxis right;
        text(1.005*stimDur, Ashift, 'S', 'fontsize', textSize, 'color', 'b');
        text(1.005*stimDur, 0, 'SN', 'fontsize', textSize, 'color', 'b');
        text(1.005*stimDur, -Ashift, 'N', 'fontsize', textSize, 'color', 'b');
        
        xlabel('time (sec)');
        yyaxis left,
        ylabel(['mean Rate (' fName(end-2:end) ')']);
        set(gca, 'ycolor', 'r');
        
        yyaxis right,
        ylabel('filtered Rate');
        set(gca, 'ycolor', 'k');
        
        title(strrep(sprintf('%s, %d Hz filtered hilb(tfs)', fName, modFreq), '_', '-'));
        set(gca, 'fontsize', fSize);
        
        set(gcf, 'units', 'inches', 'position', [1 1 9 6]);
        print([outFigDir.pdf fName], '-dpdf', '-bestfit');
        print([outFigDir.png fName], '-dpng');
    end
end