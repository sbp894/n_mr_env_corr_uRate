% in this function, plus and pos are synonymous. Same with minus and neg
function ...
    [corr_s_sn_ValsModFreq_pos, uncorr_sn_n_ValsModFreq_pos, corr_s_n_ValsModFreq_pos, corr_s_sn_ValsModFreq_neg, uncorr_sn_n_ValsModFreq_neg, corr_s_n_ValsModFreq_neg]= ... %out
    multires_modulation_uR_f0_BP... %fun
    (S_rate_plus_Org, S_rate_minus_Org, N_rate_plus_Org, N_rate_minus_Org, SN_rate_plus_Org, SN_rate_minus_Org, fsOrg, modFreq, outFigDir, curTimeResolution, combine_chi1_mu0, N_bp_half, fName, curMetaData) %in

if ~exist('outFigDir', 'var')
    outFigDir=[];
    fName=[];
end

CF_Hz= curMetaData.CF_Hz;
SNR= curMetaData.SNR;

% combine_mean0_rms1=0;

debugModeOn1=0;
debugModeOn2=1;
saveDebugFiles2= true;
plotNegPol= 0;

fs=fsOrg;
fSize=20;
textSize=20;
SetOnsetZeroDuration=1/2/modFreq;
OnsetInds2SetZero=max(1, round(fs*SetOnsetZeroDuration));

corr_s_sn_ValsModFreq_pos=nan(length(modFreq),1);
uncorr_sn_n_ValsModFreq_pos=nan(length(modFreq),1);
corr_s_n_ValsModFreq_pos=nan(length(modFreq),1);
corr_s_sn_ValsModFreq_neg=nan(length(modFreq),1);
uncorr_sn_n_ValsModFreq_neg=nan(length(modFreq),1);
corr_s_n_ValsModFreq_neg=nan(length(modFreq),1);

S_rate_plus=gen_resample(S_rate_plus_Org, fsOrg, fs);
S_rate_minus=gen_resample(S_rate_minus_Org, fsOrg, fs);
% S_rate_tfs=S_rate_plus; %(S_rate_plus-S_rate_minus)/2;

N_rate_plus=gen_resample(N_rate_plus_Org, fsOrg, fs);
N_rate_minus=gen_resample(N_rate_minus_Org, fsOrg, fs);
% N_rate_tfs=N_rate_plus; %(N_rate_plus-N_rate_minus)/2;

SN_rate_plus=gen_resample(SN_rate_plus_Org, fsOrg, fs);
SN_rate_minus=gen_resample(SN_rate_minus_Org, fsOrg, fs);
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
    curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
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
    windowVec= ones(1, numel(S_rate_plus));
    %     windowVec= hamming(numel(S_rate_plus))';
    S_rate_pos_filt=filtfilt(curFilt, S_rate_plus.*windowVec);
    S_rate_pos_filt(1:OnsetInds2SetZero)=linspace(0, 1, OnsetInds2SetZero).*S_rate_pos_filt(1:OnsetInds2SetZero);
    SN_rate_pos_filt=filtfilt(curFilt, SN_rate_plus.*windowVec);
    SN_rate_pos_filt(1:OnsetInds2SetZero)=linspace(0, 1, OnsetInds2SetZero).*SN_rate_pos_filt(1:OnsetInds2SetZero);
    N_rate_pos_filt=filtfilt(curFilt, N_rate_plus.*windowVec);
    N_rate_pos_filt(1:OnsetInds2SetZero)=linspace(0, 1, OnsetInds2SetZero).*N_rate_pos_filt(1:OnsetInds2SetZero);
    
    % for -ve polatiry
    S_rate_neg_filt=filtfilt(curFilt, S_rate_minus.*windowVec);
    S_rate_neg_filt(1:OnsetInds2SetZero)=linspace(0, 1, OnsetInds2SetZero).*S_rate_neg_filt(1:OnsetInds2SetZero);
    SN_rate_neg_filt=filtfilt(curFilt, SN_rate_minus.*windowVec);
    SN_rate_neg_filt(1:OnsetInds2SetZero)=linspace(0, 1, OnsetInds2SetZero).*SN_rate_neg_filt(1:OnsetInds2SetZero);
    N_rate_neg_filt=filtfilt(curFilt, N_rate_minus.*windowVec);
    N_rate_neg_filt(1:OnsetInds2SetZero)=linspace(0, 1, OnsetInds2SetZero).*N_rate_neg_filt(1:OnsetInds2SetZero);
    
    for windowVar=1:length(cur_s_sn_CorrVals_pos)
        tStart=(windowVar-1)*curTimeResolution;
        indStart=max(1, round(tStart*fs));
        tEnd=windowVar*curTimeResolution;
        indEnd=min(length(S_rate_plus), round(tEnd*fs));
        validINDs=indStart:indEnd;
        
        cur_s_sn_CorrVals_pos(windowVar)=corr_norm(S_rate_pos_filt(validINDs), SN_rate_pos_filt(validINDs));
        cur_sn_n_CorrVals_pos(windowVar)=corr_norm(SN_rate_pos_filt(validINDs), N_rate_pos_filt(validINDs));
        cur_s_n_CorrVals_pos(windowVar)=corr_norm(S_rate_pos_filt(validINDs), N_rate_pos_filt(validINDs));
        
        cur_s_sn_CorrVals_neg(windowVar)=corr_norm(S_rate_neg_filt(validINDs), SN_rate_neg_filt(validINDs));
        cur_sn_n_CorrVals_neg(windowVar)=corr_norm(SN_rate_neg_filt(validINDs), N_rate_neg_filt(validINDs));
        cur_s_n_CorrVals_neg(windowVar)=corr_norm(S_rate_neg_filt(validINDs), N_rate_neg_filt(validINDs));
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
    
    col.S= [.1 .1 .8];
    col.SN= [.1 .6 .1];
    col.N= [.9 .2 .2];
    
    if (~isempty(outFigDir) && ~exist([outFigDir.png fName '.png'], 'file')) || debugModeOn2
        figure(1012);
        %             if windowVar==1
        clf;
        %             end
        hold on;
        lw=1;
        lw2=2.5;
        Ashift=max(S_rate_plus(min(3*OnsetInds2SetZero, round(length(S_rate_plus)/2)):end));
        plot(stimTime, Ashift+S_rate_plus, '-', 'color', col.S, 'linew', lw);
        plot(stimTime, SN_rate_plus, '-', 'color', col.SN, 'linew', lw);
        plot(stimTime, -Ashift+N_rate_plus, '-', 'color', col.N, 'linew', lw);
        
        if plotNegPol
            plot(stimTime, Ashift+S_rate_minus, 'm:', 'linew', lw);
            plot(stimTime, -Ashift+N_rate_minus, 'm:', 'linew', lw);
            plot(stimTime, SN_rate_minus, 'm:', 'linew', lw);
        end
        axis tight;
        ylL= ylim();
        ylL_new= max(abs(ylL));
        ylim([-ylL_new*.67 ylL_new]);
        
        yyaxis right;
        lw=1.5;
        inds2checkMax_avoid_onset= min(2*OnsetInds2SetZero, round(length(S_rate_plus)/2)):length(S_rate_plus);
%         Ashift_mod=1.75*max(abs([S_rate_pos_filt(inds2checkMax_avoid_onset) SN_rate_pos_filt(inds2checkMax_avoid_onset) N_rate_pos_filt(inds2checkMax_avoid_onset)]));
        Ashift_mod=2.25*max(abs(S_rate_pos_filt(inds2checkMax_avoid_onset)));
        plot(stimTime, Ashift_mod+S_rate_pos_filt, 'k-', 'linew', lw2);
        plot(stimTime, SN_rate_pos_filt, 'k-', 'linew', lw2);
        plot(stimTime, -Ashift_mod+N_rate_pos_filt, 'k-', 'linew', lw2);
        ylR= ylim();
        ylR_new= max(abs(ylR));
        ylim([-1.75*Ashift_mod 1.75*Ashift_mod]);
        xlim([-.01 1.4]);
        
        yyaxis right;
        text(1.005*stimDur, Ashift_mod, 'S', 'fontsize', textSize, 'color', col.S);
        text(1.005*stimDur, 0, 'SN', 'fontsize', textSize, 'color', col.SN);
        text(1.005*stimDur, -Ashift_mod, 'N', 'fontsize', textSize, 'color', col.N);
        
        xlabel('time (sec)');
        yyaxis left,
        ylabel(sprintf('$PSTH_{+ve}$ (CF= %.1f kHz)', CF_Hz/1e3), 'interpreter', 'latex');
        set(gca, 'ycolor', 'r');
        
        yyaxis right,
        ylabel(sprintf('Output: $f_{mod}=%.0f$ Hz', modFreq), 'interpreter', 'latex');
        set(gca, 'ycolor', 'k');
        
        %         title(strrep(sprintf('%s, %d Hz filtered uR', fName, modFreq), '_', '-'));
        if ismember(str2double(fName(2:4)), [358 360 361 362])
            title(['HI: SNR=' num2str(SNR) ' dB']);
        else
            title(['NH: SNR=' num2str(SNR) ' dB']);
        end
        
        set(gca, 'fontsize', fSize);
        set(gcf, 'units', 'inches', 'position', [1 1 10 6]);
        
        if saveDebugFiles2
            %             print([outFigDir.pdf fName], '-dpdf', '-bestfit');
            print([outFigDir.png fName], '-dpng');
        end
        if contains(fName, {'Q362_t1_u04_SNR0_SSN', 'Q347_t1_u03_SNR0_SSN'})
            saveas(gcf, [outFigDir.latex fName], 'epsc');
        end
    end
end