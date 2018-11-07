% in this function, plus and pos are synonymous. Same with minus and neg
function [dist_s_sn_vals_pos, dist_sn_n_vals_pos, dist_s_n_vals_pos, dist_s_sn_vals_neg, dist_sn_n_vals_neg, corr_s_n_ValsModFreq_neg, snr_dist_pos, snr_dist_neg]=...
    multires_modulation_uR_euclidean(S_rate_plus_Org, S_rate_minus_Org, N_rate_plus_Org, N_rate_minus_Org, SN_rate_plus_Org, SN_rate_minus_Org, fsOrg, modFreq, outFigDir, curTimeResolution, combine_chi1_mu0, fName, CF_Hz)

if ~exist('outFigDir', 'var')
    outFigDir=[];
    fName=[];
end

debugModeOn2=0;
saveDebugFiles2= 0;
fs=2e3;
fSize=14;
textSize=20;
mrkSize= 12;

S_rate_pos=gen_resample(S_rate_plus_Org, fsOrg, fs);
S_rate_neg=gen_resample(S_rate_minus_Org, fsOrg, fs);

N_rate_pos=gen_resample(N_rate_plus_Org, fsOrg, fs);
N_rate_neg=gen_resample(N_rate_minus_Org, fsOrg, fs);

SN_rate_pos=gen_resample(SN_rate_plus_Org, fsOrg, fs);
SN_rate_neg=gen_resample(SN_rate_minus_Org, fsOrg, fs);

stimDur=length(S_rate_pos)/fs;
stimTime=(1:length(S_rate_pos))/fs;

if rem(stimDur, curTimeResolution) < .5*curTimeResolution
    s_rate_weight= nan(floor(stimDur/curTimeResolution), 1);
    s_sn_dist_pos=nan(floor(stimDur/curTimeResolution), 1);
    sn_n_dist_pos=nan(floor(stimDur/curTimeResolution), 1);
    s_n_dist_pos=nan(floor(stimDur/curTimeResolution), 1);
    s_sn_dist_neg=nan(floor(stimDur/curTimeResolution), 1);
    sn_n_dist_neg=nan(floor(stimDur/curTimeResolution), 1);
    s_n_dist_neg=nan(floor(stimDur/curTimeResolution), 1);
else
    s_rate_weight=nan(ceil(stimDur/curTimeResolution), 1); 
    s_sn_dist_pos=nan(ceil(stimDur/curTimeResolution), 1);
    sn_n_dist_pos=nan(ceil(stimDur/curTimeResolution), 1);
    s_n_dist_pos=nan(ceil(stimDur/curTimeResolution), 1);
    s_sn_dist_neg=nan(ceil(stimDur/curTimeResolution), 1);
    sn_n_dist_neg=nan(ceil(stimDur/curTimeResolution), 1);
    s_n_dist_neg=nan(ceil(stimDur/curTimeResolution), 1);
end

for windowVar=1:length(s_sn_dist_pos)
    tStart=(windowVar-1)*curTimeResolution;
    indStart=max(1, round(tStart*fs));
    tEnd=windowVar*curTimeResolution;
    indEnd=min(length(S_rate_pos), round(tEnd*fs));
    validINDs=indStart:indEnd;
    
    s_rate_weight(windowVar)= sum(S_rate_pos(validINDs));
    
    s_sn_dist_pos(windowVar)=norm(S_rate_pos(validINDs)-SN_rate_pos(validINDs));
    sn_n_dist_pos(windowVar)=norm(SN_rate_pos(validINDs)- N_rate_pos(validINDs));
    s_n_dist_pos(windowVar)=norm(S_rate_pos(validINDs)-N_rate_pos(validINDs));
    
    s_sn_dist_neg(windowVar)=norm(S_rate_neg(validINDs)-SN_rate_neg(validINDs));
    sn_n_dist_neg(windowVar)=norm(SN_rate_neg(validINDs)-N_rate_neg(validINDs));
    s_n_dist_neg(windowVar)=norm(S_rate_neg(validINDs)-N_rate_neg(validINDs));
end

maxVal= 10;
s_rate_weight=s_rate_weight/max(s_rate_weight);
snr_dist_pos= s_rate_weight.*sn_n_dist_pos./s_sn_dist_pos;
snr_dist_pos(isnan(snr_dist_pos))=0;
snr_dist_pos(isinf(snr_dist_pos))=0;
snr_dist_pos(snr_dist_pos>maxVal)=maxVal;

snr_dist_neg=s_rate_weight.*sn_n_dist_neg./s_sn_dist_neg;
snr_dist_neg(isnan(snr_dist_neg))=0;
snr_dist_neg(isinf(snr_dist_neg))=0;
snr_dist_neg(snr_dist_neg>maxVal)= maxVal;

if combine_chi1_mu0
    dist_s_sn_vals_pos=sqrt(nansum(s_sn_dist_pos.^2));
    dist_sn_n_vals_pos=sqrt(nansum(sn_n_dist_pos.^2));
    dist_s_n_vals_pos=sqrt(nansum(s_n_dist_pos.^2));
    dist_s_sn_vals_neg=sqrt(nansum(s_sn_dist_neg.^2));
    dist_sn_n_vals_neg=sqrt(nansum(sn_n_dist_neg.^2));
    corr_s_n_ValsModFreq_neg=sqrt(nansum(s_n_dist_neg.^2));
else % nanmean
    dist_s_sn_vals_pos=nanmean(s_sn_dist_pos);
    dist_sn_n_vals_pos=nanmean(sn_n_dist_pos);
    dist_s_n_vals_pos=nanmean(s_n_dist_pos);
    dist_s_sn_vals_neg=nanmean(s_sn_dist_neg);
    dist_sn_n_vals_neg=nanmean(sn_n_dist_neg);
    corr_s_n_ValsModFreq_neg=nanmean(s_n_dist_neg);
end

if ~isempty(outFigDir) && debugModeOn2 && ~exist([outFigDir.png fName '.png'], 'file')
    figure(1012);
    %             if windowVar==1
    clf;
    dist_time= (1:length(s_sn_dist_pos))*curTimeResolution-curTimeResolution/2;
    %             end
    hold on;
    lw=.8;
    Ashift=max(S_rate_pos);
    plot(stimTime, Ashift+S_rate_pos, 'b-', 'linew', lw);
    plot(stimTime, Ashift+S_rate_neg, 'c-', 'linew', lw);
    plot(stimTime, SN_rate_pos, 'b-', 'linew', lw);
    plot(stimTime, SN_rate_neg, 'c-', 'linew', lw);
    plot(stimTime, -Ashift+N_rate_pos, 'b-', 'linew', lw);
    plot(stimTime, -Ashift+N_rate_neg, 'c-', 'linew', lw);
    
    plot(dist_time, Ashift+demean_scale(s_sn_dist_pos, Ashift), 'kd-', 'linew', lw);
    plot(dist_time, demean_scale(sn_n_dist_pos, Ashift), 'kd-', 'linew', lw);
    plot(dist_time, -Ashift+demean_scale(s_n_dist_pos, Ashift), 'kd-', 'linew', lw);
    
    plot(dist_time, Ashift+demean_scale(s_sn_dist_neg, Ashift), 'kd-', 'linew', lw);
    plot(dist_time, demean_scale(sn_n_dist_neg, Ashift), 'kd-', 'linew', lw);
    plot(dist_time, -Ashift+demean_scale(s_n_dist_neg, Ashift), 'kd-', 'linew', lw);
    
    text(1.005*stimDur, Ashift, 'S vs SN', 'fontsize', textSize, 'color', 'b');
    text(1.005*stimDur, 0, 'N vs SN', 'fontsize', textSize, 'color', 'b');
    text(1.005*stimDur, -Ashift, 'S vs N', 'fontsize', textSize, 'color', 'b');
    
    
    yyaxis right;
    lw2=2;
    plot(dist_time, snr_dist_pos, 'rd-', 'linew', lw2, 'markersize', mrkSize);
    plot(dist_time, snr_dist_neg, 'md-', 'linew', lw2, 'markersize', mrkSize);
    
    xlabel('time (sec)');
    yyaxis left,
    ylabel(sprintf('mean Rate (CF= %.1f kHz)', CF_Hz/1e3));
    set(gca, 'ycolor', 'r');
    
    yyaxis right,
    ylabel('SNR (eucl)');
    set(gca, 'ycolor', 'k');
    
    title(strrep(sprintf('%s, %d Hz filtered tfs', fName, modFreq), '_', '-'));
    set(gca, 'fontsize', fSize);
    
    if saveDebugFiles2
        set(gcf, 'units', 'inches', 'position', [1 1 9 6]);
        print([outFigDir.pdf fName], '-dpdf', '-bestfit');
        print([outFigDir.png fName], '-dpng');
    end
end