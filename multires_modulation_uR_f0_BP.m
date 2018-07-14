% in this function, plus and pos are synonymous. Same with minus and neg
function [corr_s_sn_ValsModFreq_pos, uncorr_sn_n_ValsModFreq_pos, corr_s_n_ValsModFreq_pos, corr_s_sn_ValsModFreq_neg, uncorr_sn_n_ValsModFreq_neg, corr_s_n_ValsModFreq_neg]=...
    multires_modulation_uR_f0_BP(S_rate_plus_Org, S_rate_minus_Org, N_rate_plus_Org, N_rate_minus_Org, SN_rate_plus_Org, SN_rate_minus_Org, fsOrg, modFreq, outFigDir, curTimeResolution, combine_chi1_mu0, N_bp, fName)

if ~exist('outFigDir', 'var')
    outFigDir=[];
    fName=[];
end

% combine_mean0_rms1=0;

debugModeOn1=0;
debugModeOn2=0;
fs=5e3;
fSize=14;
textSize=20;
SetOnsetZeroDuration=.01;
OnsetInds2SetZero=round(fs*SetOnsetZeroDuration);

corr_s_sn_ValsModFreq_pos=nan(length(modFreq),1);
uncorr_sn_n_ValsModFreq_pos=nan(length(modFreq),1);
corr_s_n_ValsModFreq_pos=nan(length(modFreq),1);
corr_s_sn_ValsModFreq_neg=nan(length(modFreq),1);
uncorr_sn_n_ValsModFreq_neg=nan(length(modFreq),1);
corr_s_n_ValsModFreq_neg=nan(length(modFreq),1);

S_rate_plus=resample(S_rate_plus_Org, 1, fsOrg/fs);
S_rate_minus=resample(S_rate_minus_Org, 1, fsOrg/fs);
% S_rate_tfs=S_rate_plus; %(S_rate_plus-S_rate_minus)/2;

N_rate_plus=resample(N_rate_plus_Org, 1, fsOrg/fs);
N_rate_minus=resample(N_rate_minus_Org, 1, fsOrg/fs);
% N_rate_tfs=N_rate_plus; %(N_rate_plus-N_rate_minus)/2;

SN_rate_plus=resample(SN_rate_plus_Org, 1, fsOrg/fs);
SN_rate_minus=resample(SN_rate_minus_Org, 1, fsOrg/fs);
% SN_rate_tfs=SN_rate_plus; %(SN_rate_plus-SN_rate_minus)/2;

stimDur=length(S_rate_plus)/fs;
stimTime=(1:length(S_rate_plus))/fs;

for modFreqVar=1:length(modFreq)
    curModFreq=modFreq(modFreqVar);
    %     curTimeResolution=20e-3; %1/curModFreq;
    %     warning('check time resolution 8 ms vs 20 ms');
    %     N_bp=2;
    HalfPowerFrequency1=modFreq/sqrt(2);
    HalfPowerFrequency2=modFreq*sqrt(2);
    curFilt= designfilt('bandpassiir','FilterOrder',N_bp, ...get(p,props)
        'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
        'SampleRate',fs);
    
    if debugModeOn1
        fvtool(curFilt);
        xlim([curModFreq/2 2*curModFreq]/1e3);
        ylim([-200 0])
    end
    
    
    if rem(stimDur, curTimeResolution) < .5*curTimeResolution
        cur_s_sn_CorrVals_pos=nan(floor(stimDur/curTimeResolution), 1);
        cur_sn_n_CorrVals_pos=nan(floor(stimDur/curTimeResolution), 1);
        cur_s_n_CorrVals_pos=nan(floor(stimDur/curTimeResolution), 1);
        cur_s_sn_CorrVals_neg=nan(floor(stimDur/curTimeResolution), 1);
        cur_sn_n_CorrVals_neg=nan(floor(stimDur/curTimeResolution), 1);
        cur_s_n_CorrVals_neg=nan(floor(stimDur/curTimeResolution), 1);
    else
        cur_s_sn_CorrVals_pos=nan(ceil(stimDur/curTimeResolution), 1);
        cur_sn_n_CorrVals_pos=nan(ceil(stimDur/curTimeResolution), 1);
        cur_s_n_CorrVals_pos=nan(ceil(stimDur/curTimeResolution), 1);
        cur_s_sn_CorrVals_neg=nan(ceil(stimDur/curTimeResolution), 1);
        cur_sn_n_CorrVals_neg=nan(ceil(stimDur/curTimeResolution), 1);
        cur_s_n_CorrVals_neg=nan(ceil(stimDur/curTimeResolution), 1);
    end
    
    % for +ve polatiry
    S_rate_pos_filt=filtfilt(curFilt, S_rate_plus);
    S_rate_pos_filt(1:OnsetInds2SetZero)=0;
    SN_rate_pos_filt=filtfilt(curFilt, SN_rate_plus);
    SN_rate_pos_filt(1:OnsetInds2SetZero)=0;
    N_rate_pos_filt=filtfilt(curFilt, N_rate_plus);
    N_rate_pos_filt(1:OnsetInds2SetZero)=0;
    
    % for -ve polatiry
    S_rate_neg_filt=filtfilt(curFilt, S_rate_minus);
    S_rate_neg_filt(1:OnsetInds2SetZero)=0;
    SN_rate_neg_filt=filtfilt(curFilt, SN_rate_minus);
    SN_rate_neg_filt(1:OnsetInds2SetZero)=0;
    N_rate_neg_filt=filtfilt(curFilt, N_rate_minus);
    N_rate_neg_filt(1:OnsetInds2SetZero)=0;
    
    for windowVar=1:length(cur_s_sn_CorrVals_pos)
        tStart=(windowVar-1)*curTimeResolution;
        indStart=max(1, round(tStart*fs));
        tEnd=windowVar*curTimeResolution;
        indEnd=min(length(S_rate_plus), round(tEnd*fs));
        validINDs=indStart:indEnd;
        
        cur_s_sn_CorrVals_pos(windowVar)=corr2(S_rate_pos_filt(validINDs), SN_rate_pos_filt(validINDs));
        cur_sn_n_CorrVals_pos(windowVar)=corr2(SN_rate_pos_filt(validINDs), N_rate_pos_filt(validINDs));
        cur_s_n_CorrVals_pos(windowVar)=corr2(S_rate_pos_filt(validINDs), N_rate_pos_filt(validINDs));
        
        cur_s_sn_CorrVals_neg(windowVar)=corr2(S_rate_neg_filt(validINDs), SN_rate_neg_filt(validINDs));
        cur_sn_n_CorrVals_neg(windowVar)=corr2(SN_rate_neg_filt(validINDs), N_rate_neg_filt(validINDs));
        cur_s_n_CorrVals_neg(windowVar)=corr2(S_rate_neg_filt(validINDs), N_rate_neg_filt(validINDs));
    end
    
    cur_s_sn_CorrVals_pos(cur_s_sn_CorrVals_pos<0)=0;
    cur_sn_n_CorrVals_pos(cur_sn_n_CorrVals_pos<0)=0;
    cur_sn_n_unCorrVals_pos=1-cur_sn_n_CorrVals_pos;
    cur_s_n_CorrVals_pos(cur_s_n_CorrVals_pos<0)=0;
    
    cur_s_sn_CorrVals_neg(cur_s_sn_CorrVals_neg<0)=0;
    cur_sn_n_CorrVals_neg(cur_sn_n_CorrVals_neg<0)=0;
    cur_sn_n_unCorrVals_neg=1-cur_sn_n_CorrVals_neg;
    cur_s_n_CorrVals_neg(cur_s_n_CorrVals_neg<0)=0;
    
    if combine_chi1_mu0
        corr_s_sn_ValsModFreq_pos(modFreqVar)=sqrt(nansum(cur_s_sn_CorrVals_pos.^2));
        uncorr_sn_n_ValsModFreq_pos(modFreqVar)=sqrt(nansum(cur_sn_n_unCorrVals_pos.^2));
        corr_s_n_ValsModFreq_pos(modFreqVar)=sqrt(nansum(cur_s_n_CorrVals_pos.^2));
        corr_s_sn_ValsModFreq_neg(modFreqVar)=sqrt(nansum(cur_s_sn_CorrVals_neg.^2));
        uncorr_sn_n_ValsModFreq_neg(modFreqVar)=sqrt(nansum(cur_sn_n_unCorrVals_neg.^2));
        corr_s_n_ValsModFreq_neg(modFreqVar)=sqrt(nansum(cur_s_n_CorrVals_neg.^2));
    else % nanmean
        corr_s_sn_ValsModFreq_pos(modFreqVar)=nanmean(cur_s_sn_CorrVals_pos);
        uncorr_sn_n_ValsModFreq_pos(modFreqVar)=nanmean(cur_sn_n_unCorrVals_pos);
        corr_s_n_ValsModFreq_pos(modFreqVar)=nanmean(cur_s_n_CorrVals_pos);
        corr_s_sn_ValsModFreq_neg(modFreqVar)=nanmean(cur_s_sn_CorrVals_neg);
        uncorr_sn_n_ValsModFreq_neg(modFreqVar)=nanmean(cur_sn_n_unCorrVals_neg);
        corr_s_n_ValsModFreq_neg(modFreqVar)=nanmean(cur_s_n_CorrVals_neg);
    end
    
    
    if ~isempty(outFigDir) && debugModeOn2
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
        Ashift=1.75*max(abs([S_rate_pos_filt SN_rate_pos_filt N_rate_pos_filt]));
        plot(stimTime, Ashift+S_rate_pos_filt, 'k-', 'linew', lw);
        plot(stimTime, SN_rate_pos_filt, 'k-', 'linew', lw);
        plot(stimTime, -Ashift+N_rate_pos_filt, 'k-', 'linew', lw);
        
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
        
        title(strrep(sprintf('%s, %d Hz filtered tfs', fName, modFreq), '_', '-'));
        set(gca, 'fontsize', fSize);
        
        set(gcf, 'units', 'inches', 'position', [1 1 9 6]);
        print([outFigDir.pdf fName], '-dpdf', '-bestfit');
        print([outFigDir.png fName], '-dpng');
    end
end