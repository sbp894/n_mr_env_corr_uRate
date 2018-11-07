% in this function, plus and pos are synonymous. Same with minus and neg
function [corr_s_sn_ValsModFreq_pos, uncorr_sn_n_ValsModFreq_pos, corr_s_n_ValsModFreq_pos, corr_s_sn_ValsModFreq_neg, uncorr_sn_n_ValsModFreq_neg, corr_s_n_ValsModFreq_neg]=...
    multires_modulation_uR_f0_LP(S_rate_plus_Org, S_rate_minus_Org, N_rate_plus_Org, N_rate_minus_Org, SN_rate_plus_Org, SN_rate_minus_Org, fsOrg, modFreq, outFigDir, curTimeResolution, combine_chi1_mu0, N_lp, fName, CF_Hz)

if ~exist('outFigDir', 'var')
    outFigDir=[];
    fName=[];
end

% combine_mean0_rms1=0;

debugModeOn1=0;
debugModeOn2=1;
saveDebugFiles2= 1;
fs=2e3;
fSize=14;
textSize=20;
SetOnsetZeroDuration=.01;
OnsetInds2SetZero=round(fs*SetOnsetZeroDuration);

corrCoef0Corr1= 1;

%% find p and q for resampling
temp= rats(fs/fsOrg);
temp(temp==' ')= [];
if ~isempty(find(temp=='/', 1))
    pqEndPartition= find(temp=='/');
    p= str2double(temp(1:pqEndPartition-1));
    q= str2double(temp(pqEndPartition+1:end));
else
    p= str2double(temp);
    q= 1;
end

S_rate_plus=resample(S_rate_plus_Org, p, q);
S_rate_minus=resample(S_rate_minus_Org, p, q);
% S_rate_tfs=S_rate_plus; %(S_rate_plus-S_rate_minus)/2;

N_rate_plus=resample(N_rate_plus_Org, p, q);
N_rate_minus=resample(N_rate_minus_Org, p, q);
% N_rate_tfs=N_rate_plus; %(N_rate_plus-N_rate_minus)/2;

SN_rate_plus=resample(SN_rate_plus_Org, p, q);
SN_rate_minus=resample(SN_rate_minus_Org, p, q);
% SN_rate_tfs=SN_rate_plus; %(SN_rate_plus-SN_rate_minus)/2;

stimDur=length(S_rate_plus)/fs;
stimTime=(1:length(S_rate_plus))/fs;

%     curTimeResolution=20e-3; %1/curModFreq;
%     warning('check time resolution 8 ms vs 20 ms');
%     N_bp=2;

curFilt= designfilt('lowpassiir','FilterOrder',N_lp, ...
    'PassbandFrequency',modFreq,'PassbandRipple',0.2, ...
    'SampleRate',fs);

if debugModeOn1
    fvtool(curFilt);
    xlim([0 2*modFreq]);
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
    
    if corrCoef0Corr1
        cur_s_sn_CorrVals_pos(windowVar)=corr2(S_rate_pos_filt(validINDs), SN_rate_pos_filt(validINDs));
        cur_sn_n_CorrVals_pos(windowVar)=corr2(SN_rate_pos_filt(validINDs), N_rate_pos_filt(validINDs));
        cur_s_n_CorrVals_pos(windowVar)=corr2(S_rate_pos_filt(validINDs), N_rate_pos_filt(validINDs));
        
        cur_s_sn_CorrVals_neg(windowVar)=corr2(S_rate_neg_filt(validINDs), SN_rate_neg_filt(validINDs));
        cur_sn_n_CorrVals_neg(windowVar)=corr2(SN_rate_neg_filt(validINDs), N_rate_neg_filt(validINDs));
        cur_s_n_CorrVals_neg(windowVar)=corr2(S_rate_neg_filt(validINDs), N_rate_neg_filt(validINDs));
        
    else
        cur_s_sn_CorrVals_pos(windowVar)=corr2(S_rate_pos_filt(validINDs), SN_rate_pos_filt(validINDs));
        cur_sn_n_CorrVals_pos(windowVar)=corr2(SN_rate_pos_filt(validINDs), N_rate_pos_filt(validINDs));
        cur_s_n_CorrVals_pos(windowVar)=corr2(S_rate_pos_filt(validINDs), N_rate_pos_filt(validINDs));
        
        cur_s_sn_CorrVals_neg(windowVar)=corr2(S_rate_neg_filt(validINDs), SN_rate_neg_filt(validINDs));
        cur_sn_n_CorrVals_neg(windowVar)=corr2(SN_rate_neg_filt(validINDs), N_rate_neg_filt(validINDs));
        cur_s_n_CorrVals_neg(windowVar)=corr2(S_rate_neg_filt(validINDs), N_rate_neg_filt(validINDs));
    end
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
    corr_s_sn_ValsModFreq_pos=sqrt(nansum(cur_s_sn_CorrVals_pos.^2));
    uncorr_sn_n_ValsModFreq_pos=sqrt(nansum(cur_sn_n_unCorrVals_pos.^2));
    corr_s_n_ValsModFreq_pos=sqrt(nansum(cur_s_n_CorrVals_pos.^2));
    corr_s_sn_ValsModFreq_neg=sqrt(nansum(cur_s_sn_CorrVals_neg.^2));
    uncorr_sn_n_ValsModFreq_neg=sqrt(nansum(cur_sn_n_unCorrVals_neg.^2));
    corr_s_n_ValsModFreq_neg=sqrt(nansum(cur_s_n_CorrVals_neg.^2));
else % nanmean
    corr_s_sn_ValsModFreq_pos=nanmean(cur_s_sn_CorrVals_pos);
    uncorr_sn_n_ValsModFreq_pos=nanmean(cur_sn_n_unCorrVals_pos);
    corr_s_n_ValsModFreq_pos=nanmean(cur_s_n_CorrVals_pos);
    corr_s_sn_ValsModFreq_neg=nanmean(cur_s_sn_CorrVals_neg);
    uncorr_sn_n_ValsModFreq_neg=nanmean(cur_sn_n_unCorrVals_neg);
    corr_s_n_ValsModFreq_neg=nanmean(cur_s_n_CorrVals_neg);
end


if ~isempty(outFigDir) && debugModeOn2 && ~exist([outFigDir.png fName '.png'], 'file')
    %             figure(1012);
    %             if windowVar==1
    clf;
    %             end
    hold on;
    lw=1.5;
    Ashift=1.25*max(S_rate_plus);
    plot(stimTime, Ashift+S_rate_plus, 'b-.', 'linew', lw);
    plot(stimTime, Ashift+S_rate_minus, 'm:', 'linew', lw);
    plot(stimTime, SN_rate_plus, 'b-.', 'linew', lw);
    plot(stimTime, SN_rate_minus, 'm:', 'linew', lw);
    plot(stimTime, -Ashift+N_rate_plus, 'b-.', 'linew', lw);
    plot(stimTime, -Ashift+N_rate_minus, 'm:', 'linew', lw);
    
    yyaxis right;
    lw=2.5;
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
    ylabel(sprintf('mean Rate (CF= %.1f kHz)', CF_Hz/1e3));
    set(gca, 'ycolor', 'r');
    
    yyaxis right,
    ylabel('filtered Rate');
    set(gca, 'ycolor', 'k');
    
    title(strrep(sprintf('%s, %d Hz filtered +PSTH', fName, modFreq), '_', '-'));
    set(gca, 'fontsize', fSize);
    grid on;
   
    if saveDebugFiles2
        set(gcf, 'units', 'inches', 'position', [1 1 9 6]);
        print([outFigDir.pdf fName], '-dpdf', '-bestfit');
        print([outFigDir.png fName], '-dpng');
    end
end