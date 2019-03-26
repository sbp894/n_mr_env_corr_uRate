function [S_rate_pos_filt, S_rate_neg_filt, time]=mod_filter(S_rate_plus_Org, S_rate_minus_Org, fsOrg, modFreq, N_bp_half, fs)

debugModeOn1=0;
SetOnsetZeroDuration=.01;
OnsetInds2SetZero=max(1, round(fs*SetOnsetZeroDuration));

S_rate_plus=resample(S_rate_plus_Org, 1, fsOrg/fs);
S_rate_minus=resample(S_rate_minus_Org, 1, fsOrg/fs);
% S_rate_tfs=S_rate_plus; %(S_rate_plus-S_rate_minus)/2;

warning('Debugging! Uncomment the following line');

HalfPowerFrequency1=1 ;% modFreq/sqrt(2);
HalfPowerFrequency2=modFreq*sqrt(2);
curFilt= designfilt('bandpassiir','FilterOrder',N_bp_half, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs);

if debugModeOn1
    fvtool(curFilt);
    xlim([modFreq/2 2*modFreq]/1e3);
    ylim([-200 0])
end

% for +ve polatiry
windowVec= ones(1, numel(S_rate_plus));
%     windowVec= hamming(numel(S_rate_plus))';
S_rate_pos_filt=filtfilt(curFilt, S_rate_plus.*windowVec);
S_rate_pos_filt(1:OnsetInds2SetZero)=0;

% for -ve polatiry
S_rate_neg_filt=filtfilt(curFilt, S_rate_minus);
S_rate_neg_filt(1:OnsetInds2SetZero)=0;

time= (1:length(S_rate_pos_filt))/fs;