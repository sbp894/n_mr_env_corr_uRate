function plot_population_plots(mr_corr_Data, modFreq, combine_chi1_mu0, postfix_fName, TimeResolution, saveAllFigs)


mrkSize=12;
hanOneModFreqCorr_S_SN= 31;
hanOneModFrequnCorr_SN_N=32;
hanOneModFreq_Corr_ratio=33;
fSize=14;
leg_fSize=7;

lw2=2;
lw3=3;
xShift=.1;

plot_S_SN_corr= 0;

%% Compute correlation between R(S) and R(SN)
if plot_S_SN_corr
    figure(hanOneModFreqCorr_S_SN);
    clf;
    
    % subplot 1
    subplot(121);
    hold on;
    plot(nan, nan, '-dr', 'linew',2)
    plot(nan, nan, 'm<--', 'linew',2)
    plot(nan, nan, '-', 'color', .5*ones(1,3), 'linew',2)
    plot(nan, nan, '-bd', 'linew',2)
    plot(nan, nan, 'c>--', 'linew',2)
    plot(nan, nan, 'k-', 'linew',2)
    
    snrs=unique([mr_corr_Data.SNR]);
    FLN_corr_NF_dPrime=zeros(length(snrs),1);
    SSN_corr_NF_dPrime=zeros(length(snrs),1);
    
    
    for snrVar=1:length(snrs)
        curSNR=snrs(snrVar);
        curSNRinds=find([mr_corr_Data.SNR]==curSNR);
        PitchModIndex=dsearchn(modFreq', modFreq);
        
        FLNcorr_pitchRes=[mr_corr_Data(curSNRinds).FLNcorr_s_sn_pos]; %
        FLNcorr_pitchRes=FLNcorr_pitchRes(PitchModIndex, :);
        plot(curSNR*ones(1, length(curSNRinds))-2*xShift, FLNcorr_pitchRes, 'bd');
        
        FLNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).FLNcorr_s_n_pos];
        FLNcorr_pitchRes_nf=FLNcorr_pitchRes_nf(PitchModIndex, :);
        plot(curSNR*ones(1, length(curSNRinds))-xShift, FLNcorr_pitchRes_nf, 'c>');
        
        SSNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).SSNcorr_s_n_pos];
        SSNcorr_pitchRes_nf=SSNcorr_pitchRes_nf(PitchModIndex, :);
        plot(curSNR*ones(1, length(curSNRinds))+xShift,SSNcorr_pitchRes_nf, 'm<');
        
        SSNcorr_pitchRes=[mr_corr_Data(curSNRinds).SSNcorr_s_sn_pos]; %
        SSNcorr_pitchRes=SSNcorr_pitchRes(PitchModIndex, :);
        plot(curSNR*ones(1, length(curSNRinds))+2*xShift, SSNcorr_pitchRes, 'rd');
        
        FLN_corr_NF_dPrime(snrVar,1)=mean(FLNcorr_pitchRes);
        FLN_corr_NF_dPrime(snrVar,2)=mean(FLNcorr_pitchRes_nf);
        FLN_corr_NF_dPrime(snrVar,3)=calc_dPrime(FLNcorr_pitchRes, FLNcorr_pitchRes_nf);
        
        SSN_corr_NF_dPrime(snrVar,1)=mean(SSNcorr_pitchRes);
        SSN_corr_NF_dPrime(snrVar,2)=mean(SSNcorr_pitchRes_nf);
        SSN_corr_NF_dPrime(snrVar,3)=calc_dPrime(SSNcorr_pitchRes, SSNcorr_pitchRes_nf);
    end
    set(gca, 'xtick', snrs);
    plot(-2*xShift+snrs, FLN_corr_NF_dPrime(:, 1), 'b-o', 'linew', lw3, 'markersize', mrkSize);
    plot(2*xShift+snrs, SSN_corr_NF_dPrime(:, 1), 'r-o', 'linew', lw3, 'markersize', mrkSize);
    plot(-xShift+snrs, FLN_corr_NF_dPrime(:, 2), 'oc--', 'linew', lw2, 'markersize', mrkSize);
    plot(xShift+snrs, SSN_corr_NF_dPrime(:, 2), 'om--', 'linew', lw2, 'markersize', mrkSize);
    grid on;
    
    yyaxis right;
    set(gca, 'ycolor', 'k');
    plot(snrs, FLN_corr_NF_dPrime(:, 3), 'd-k', 'linew', 2, 'markersize', mrkSize);
    plot(snrs, SSN_corr_NF_dPrime(:, 3), 'd-', 'color', .5*ones(1,3), 'linew', 2, 'markersize', mrkSize);
    
    xlabel('SNR (dB)');
    if combine_chi1_mu0
        yyaxis left, ylabel('$\overline{\chi ^2}$', 'interpreter', 'latex');
    else
        yyaxis left, ylabel('$\overline{R_{f0,BP}(S).R_{f0,BP}(SN)}$', 'interpreter', 'latex');
        ylim([0 1]);
    end
    
    yyaxis right; ylabel('d-prime');
    title(sprintf('ModFreq=%.0f Hz (corr(S,SN)), +ve uR', modFreq));
    % legend('SSN', 'SSN-NF', 'SSN-dPrime', 'FLN', 'FLN-NF', 'FLN-dPrime', 'location', 'northwest', 'FontSize',10);
    set(gca, 'fontsize', fSize);
    
    % subplot 2
    subplot(122);
    hold on;
    plot(nan, nan, '-dr', 'linew',2)
    plot(nan, nan, 'm<--', 'linew',2)
    plot(nan, nan, '-', 'color', .5*ones(1,3), 'linew',2)
    plot(nan, nan, '-bd', 'linew',2)
    plot(nan, nan, 'c>--', 'linew',2)
    plot(nan, nan, 'k-', 'linew',2)
    
    snrs=unique([mr_corr_Data.SNR]);
    FLN_corr_NF_dPrime=zeros(length(snrs),1);
    SSN_corr_NF_dPrime=zeros(length(snrs),1);
    
    lw2=2;
    lw3=3;
    xShift=.1;
    
    for snrVar=1:length(snrs)
        curSNR=snrs(snrVar);
        curSNRinds=find([mr_corr_Data.SNR]==curSNR);
        PitchModIndex=dsearchn(modFreq', modFreq);
        
        FLNcorr_pitchRes=[mr_corr_Data(curSNRinds).FLNcorr_s_sn_neg]; %
        FLNcorr_pitchRes=FLNcorr_pitchRes(PitchModIndex, :);
        plot(curSNR*ones(1, length(curSNRinds))-2*xShift, FLNcorr_pitchRes, 'bd');
        
        FLNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).FLNcorr_s_n_neg];
        FLNcorr_pitchRes_nf=FLNcorr_pitchRes_nf(PitchModIndex, :);
        plot(curSNR*ones(1, length(curSNRinds))-xShift, FLNcorr_pitchRes_nf, 'c>');
        
        SSNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).SSNcorr_s_n_neg];
        SSNcorr_pitchRes_nf=SSNcorr_pitchRes_nf(PitchModIndex, :);
        plot(curSNR*ones(1, length(curSNRinds))+xShift,SSNcorr_pitchRes_nf, 'm<');
        
        SSNcorr_pitchRes=[mr_corr_Data(curSNRinds).SSNcorr_s_sn_neg]; %
        SSNcorr_pitchRes=SSNcorr_pitchRes(PitchModIndex, :);
        plot(curSNR*ones(1, length(curSNRinds))+2*xShift, SSNcorr_pitchRes, 'rd');
        
        FLN_corr_NF_dPrime(snrVar,1)=mean(FLNcorr_pitchRes);
        FLN_corr_NF_dPrime(snrVar,2)=mean(FLNcorr_pitchRes_nf);
        FLN_corr_NF_dPrime(snrVar,3)=calc_dPrime(FLNcorr_pitchRes, FLNcorr_pitchRes_nf);
        
        SSN_corr_NF_dPrime(snrVar,1)=mean(SSNcorr_pitchRes);
        SSN_corr_NF_dPrime(snrVar,2)=mean(SSNcorr_pitchRes_nf);
        SSN_corr_NF_dPrime(snrVar,3)=calc_dPrime(SSNcorr_pitchRes, SSNcorr_pitchRes_nf);
    end
    set(gca, 'xtick', snrs);
    plot(-2*xShift+snrs, FLN_corr_NF_dPrime(:, 1), 'b-o', 'linew', lw3, 'markersize', mrkSize);
    plot(2*xShift+snrs, SSN_corr_NF_dPrime(:, 1), 'r-o', 'linew', lw3, 'markersize', mrkSize);
    plot(-xShift+snrs, FLN_corr_NF_dPrime(:, 2), 'oc--', 'linew', lw2, 'markersize', mrkSize);
    plot(xShift+snrs, SSN_corr_NF_dPrime(:, 2), 'om--', 'linew', lw2, 'markersize', mrkSize);
    grid on;
    
    yyaxis right;
    set(gca, 'ycolor', 'k');
    plot(snrs, FLN_corr_NF_dPrime(:, 3), 'd-k', 'linew', 2, 'markersize', mrkSize);
    plot(snrs, SSN_corr_NF_dPrime(:, 3), 'd-', 'color', .5*ones(1,3), 'linew', 2, 'markersize', mrkSize);
    
    xlabel('SNR (dB)');
    if combine_chi1_mu0
        yyaxis left, ylabel('$\overline{\chi ^2}$', 'interpreter', 'latex');
    else
        yyaxis left, ylabel('$\overline{R_{f0,BP}(S).R_{f0,BP}(SN)}$', 'interpreter', 'latex');
        ylim([0 1]);
    end
    
    yyaxis right; ylabel('d-prime');
    title(sprintf('ModFreq=%.0f Hz (corr(S,SN)), -ve uR', modFreq));
    LG1=legend('SSN', 'SSN-NF', 'SSN-dPrime', 'FLN', 'FLN-NF', 'FLN-dPrime', 'location', 'northwest');
    LG1.FontSize=leg_fSize;
    set(gca, 'fontsize', fSize);
    
    
    set(hanOneModFreqCorr_S_SN, 'units', 'inches', 'position', [1 1 14 5]);
    fName_singleModFreqCorr=sprintf('sEPSM_uR_f0_BP_modFreq%.0fHz_corr_%s_tRes%.0fms', modFreq, postfix_fName, TimeResolution*1e3);
    if saveAllFigs
        saveas(hanOneModFreqCorr_S_SN, [saving_Dir fName_singleModFreqCorr], 'tiff');
        saveas(hanOneModFreqCorr_S_SN, [Latex_Dir fName_singleModFreqCorr], 'epsc');
    end
    
end


%% check single modualtion frequency
figure(hanOneModFrequnCorr_SN_N);
clf;

% subplot 1
subplot(121);
hold on;
plot(nan, nan, '-dr', 'linew',2)
plot(nan, nan, 'm<--', 'linew',2)
plot(nan, nan, '-', 'color', .5*ones(1,3), 'linew',2)
plot(nan, nan, '-bd', 'linew',2)
plot(nan, nan, 'c>--', 'linew',2)
plot(nan, nan, 'k-', 'linew',2)

snrs=unique([mr_corr_Data.SNR]);
FLN_corr_NF_dPrime=zeros(length(snrs),1);
SSN_corr_NF_dPrime=zeros(length(snrs),1);
for snrVar=1:length(snrs)
    curSNR=snrs(snrVar);
    curSNRinds=find([mr_corr_Data.SNR]==curSNR);
    PitchModIndex=dsearchn(modFreq', modFreq);
    
    FLNuncorr_pitchRes=[mr_corr_Data(curSNRinds).FLNuncorr_sn_n_pos]; %
    FLNuncorr_pitchRes=FLNuncorr_pitchRes(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))-2*xShift, FLNuncorr_pitchRes, 'bd');
    
    FLNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).FLNcorr_s_n_pos];
    FLNcorr_pitchRes_nf=FLNcorr_pitchRes_nf(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))-xShift, FLNcorr_pitchRes_nf, 'c>');
    
    SSNuncorr_pitchRes=[mr_corr_Data(curSNRinds).SSNuncorr_sn_n_pos]; %
    SSNuncorr_pitchRes=SSNuncorr_pitchRes(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))+2*xShift, SSNuncorr_pitchRes, 'rd');
    
    SSNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).SSNcorr_s_n_pos];
    SSNcorr_pitchRes_nf=SSNcorr_pitchRes_nf(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))+xShift,SSNcorr_pitchRes_nf, 'm<');
    
    FLN_corr_NF_dPrime(snrVar,1)=mean(FLNuncorr_pitchRes);
    FLN_corr_NF_dPrime(snrVar,2)=mean(FLNcorr_pitchRes_nf);
    FLN_corr_NF_dPrime(snrVar,3)=calc_dPrime(FLNuncorr_pitchRes, FLNcorr_pitchRes_nf);
    
    SSN_corr_NF_dPrime(snrVar,1)=mean(SSNuncorr_pitchRes);
    SSN_corr_NF_dPrime(snrVar,2)=mean(SSNcorr_pitchRes_nf);
    SSN_corr_NF_dPrime(snrVar,3)=calc_dPrime(SSNuncorr_pitchRes, SSNcorr_pitchRes_nf);
end
set(gca, 'xtick', snrs);
plot(-2*xShift+snrs, FLN_corr_NF_dPrime(:, 1), 'b-o', 'linew', lw3, 'markersize', mrkSize);
plot(2*xShift+snrs, SSN_corr_NF_dPrime(:, 1), 'r-o', 'linew', lw3, 'markersize', mrkSize);
plot(-xShift+snrs, FLN_corr_NF_dPrime(:, 2), 'oc--', 'linew', lw2, 'markersize', mrkSize);
plot(xShift+snrs, SSN_corr_NF_dPrime(:, 2), 'om--', 'linew', lw2, 'markersize', mrkSize);
grid on;

yyaxis right;
set(gca, 'ycolor', 'k');
plot(snrs, FLN_corr_NF_dPrime(:, 3), 'd-k', 'linew', 2, 'markersize', mrkSize);
plot(snrs, SSN_corr_NF_dPrime(:, 3), 'd-', 'color', .5*ones(1,3), 'linew', 2, 'markersize', mrkSize);

xlabel('SNR (dB)');
if combine_chi1_mu0
    yyaxis left, ylabel('$\overline{\chi ^2}$', 'interpreter', 'latex');
else
    yyaxis left, ylabel('$1-\overline{R_{f0,LP}(SN).R_{f0,LP}(N)}$', 'interpreter', 'latex');
    ylim([0 1]);
end
yyaxis right; ylabel('d-prime');
title(sprintf('ModFreq=%.0f Hz (1-corr(SN,N)), +ve uR', modFreq));
% legend('SSN', 'SSN-NF', 'SSN-dPrime', 'FLN', 'FLN-NF', 'FLN-dPrime', 'location', 'northwest', 'FontSize',10);
set(gca, 'fontsize', fSize);

% subplot 2
subplot(122);
hold on;
plot(nan, nan, '-dr', 'linew',2)
plot(nan, nan, 'm<--', 'linew',2)
plot(nan, nan, '-', 'color', .5*ones(1,3), 'linew',2)
plot(nan, nan, '-bd', 'linew',2)
plot(nan, nan, 'c>--', 'linew',2)
plot(nan, nan, 'k-', 'linew',2)

snrs=unique([mr_corr_Data.SNR]);
FLN_corr_NF_dPrime=zeros(length(snrs),1);
SSN_corr_NF_dPrime=zeros(length(snrs),1);
for snrVar=1:length(snrs)
    curSNR=snrs(snrVar);
    curSNRinds=find([mr_corr_Data.SNR]==curSNR);
    PitchModIndex=dsearchn(modFreq', modFreq);
    
    FLNuncorr_pitchRes=[mr_corr_Data(curSNRinds).FLNuncorr_sn_n_neg]; %
    FLNuncorr_pitchRes=FLNuncorr_pitchRes(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))-2*xShift, FLNuncorr_pitchRes, 'bd');
    
    FLNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).FLNcorr_s_n_neg];
    FLNcorr_pitchRes_nf=FLNcorr_pitchRes_nf(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))-xShift, FLNcorr_pitchRes_nf, 'c>');
    
    SSNuncorr_pitchRes=[mr_corr_Data(curSNRinds).SSNuncorr_sn_n_neg]; %
    SSNuncorr_pitchRes=SSNuncorr_pitchRes(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))+2*xShift, SSNuncorr_pitchRes, 'rd');
    
    SSNcorr_pitchRes_nf=[mr_corr_Data(curSNRinds).SSNcorr_s_n_neg];
    SSNcorr_pitchRes_nf=SSNcorr_pitchRes_nf(PitchModIndex, :);
    plot(curSNR*ones(1, length(curSNRinds))+xShift,SSNcorr_pitchRes_nf, 'm<');
    
    FLN_corr_NF_dPrime(snrVar,1)=mean(FLNuncorr_pitchRes);
    FLN_corr_NF_dPrime(snrVar,2)=mean(FLNcorr_pitchRes_nf);
    FLN_corr_NF_dPrime(snrVar,3)=calc_dPrime(FLNuncorr_pitchRes, FLNcorr_pitchRes_nf);
    
    SSN_corr_NF_dPrime(snrVar,1)=mean(SSNuncorr_pitchRes);
    SSN_corr_NF_dPrime(snrVar,2)=mean(SSNcorr_pitchRes_nf);
    SSN_corr_NF_dPrime(snrVar,3)=calc_dPrime(SSNuncorr_pitchRes, SSNcorr_pitchRes_nf);
end
set(gca, 'xtick', snrs);
plot(-2*xShift+snrs, FLN_corr_NF_dPrime(:, 1), 'b-o', 'linew', lw3, 'markersize', mrkSize);
plot(2*xShift+snrs, SSN_corr_NF_dPrime(:, 1), 'r-o', 'linew', lw3, 'markersize', mrkSize);
plot(-xShift+snrs, FLN_corr_NF_dPrime(:, 2), 'oc--', 'linew', lw2, 'markersize', mrkSize);
plot(xShift+snrs, SSN_corr_NF_dPrime(:, 2), 'om--', 'linew', lw2, 'markersize', mrkSize);
grid on;

yyaxis right;
set(gca, 'ycolor', 'k');
plot(snrs, FLN_corr_NF_dPrime(:, 3), 'd-k', 'linew', 2, 'markersize', mrkSize);
plot(snrs, SSN_corr_NF_dPrime(:, 3), 'd-', 'color', .5*ones(1,3), 'linew', 2, 'markersize', mrkSize);

xlabel('SNR (dB)');
if combine_chi1_mu0
    yyaxis left, ylabel('$\overline{\chi ^2}$', 'interpreter', 'latex');
else
    yyaxis left, ylabel('$1-\overline{R_{f0,LP}(SN).R_{f0,LP}(N)}$', 'interpreter', 'latex');
    ylim([0 1]);
end
yyaxis right; ylabel('d-prime');
title(sprintf('ModFreq=%.0f Hz (1-corr(SN,N)), -ve uR', modFreq));
LG2=legend('SSN', 'SSN-NF', 'SSN-dPrime', 'FLN', 'FLN-NF', 'FLN-dPrime', 'location', 'northwest');
LG2.FontSize=leg_fSize;
set(gca, 'fontsize', fSize);


%% Compute {correlation between R(SN) and R(S)}/{correlation between R(SN) and R(N)}
figure(hanOneModFreq_Corr_ratio);
clf;

% subplot 1
subplot(121);
hold on;
subplot(122);
hold on;


snrs=unique([mr_corr_Data.SNR]);
FLN_sn_ratio_pos_mean=zeros(length(snrs),1);
SSN_sn_ratio_pos_mean=zeros(length(snrs),1);
FLN_sn_ratio_neg_mean=zeros(length(snrs),1);
SSN_sn_ratio_neg_mean=zeros(length(snrs),1);

for snrVar=1:length(snrs)
    curSNR=snrs(snrVar);
    curSNRinds=find([mr_corr_Data.SNR]==curSNR);
    PitchModIndex=dsearchn(modFreq', modFreq);
    
    subplot(121);
    
    FLN_s_sn_corr_pos= [mr_corr_Data(curSNRinds).FLNcorr_s_sn_pos]; %
    FLN_s_sn_corr_pos=FLN_s_sn_corr_pos(PitchModIndex, :);
    FLN_n_sn_corr_pos= 1-[mr_corr_Data(curSNRinds).FLNuncorr_sn_n_pos]; %
    FLN_n_sn_corr_pos=FLN_n_sn_corr_pos(PitchModIndex, :);
    FLN_sn_ratio_pos= FLN_s_sn_corr_pos./FLN_n_sn_corr_pos;
    plot(curSNR*ones(1, length(curSNRinds))-2*xShift, FLN_sn_ratio_pos, 'bd');
    
    SSN_s_sn_corr_pos= [mr_corr_Data(curSNRinds).SSNcorr_s_sn_pos]; %
    SSN_s_sn_corr_pos=SSN_s_sn_corr_pos(PitchModIndex, :);
    SSN_n_sn_corr_pos= 1-[mr_corr_Data(curSNRinds).SSNuncorr_sn_n_pos]; %
    SSN_n_sn_corr_pos=SSN_n_sn_corr_pos(PitchModIndex, :);
    SSN_sn_ratio_pos= SSN_s_sn_corr_pos./SSN_n_sn_corr_pos;
    plot(curSNR*ones(1, length(curSNRinds))+2*xShift, SSN_sn_ratio_pos, 'rd');
    
    FLN_sn_ratio_pos_mean(snrVar)= mean(FLN_sn_ratio_pos);
    SSN_sn_ratio_pos_mean(snrVar)= mean(SSN_sn_ratio_pos);
    
    subplot(122);
    
    FLN_s_sn_corr_neg= [mr_corr_Data(curSNRinds).FLNcorr_s_sn_neg]; %
    FLN_s_sn_corr_neg=FLN_s_sn_corr_neg(PitchModIndex, :);
    FLN_n_sn_corr_neg= 1-[mr_corr_Data(curSNRinds).FLNuncorr_sn_n_neg]; %
    FLN_n_sn_corr_neg=FLN_n_sn_corr_neg(PitchModIndex, :);
    FLN_sn_ratio_neg= FLN_s_sn_corr_neg./FLN_n_sn_corr_neg;
    plot(curSNR*ones(1, length(curSNRinds))-2*xShift, FLN_sn_ratio_neg, 'bd');
    
    SSN_s_sn_corr_neg= [mr_corr_Data(curSNRinds).SSNcorr_s_sn_neg]; %
    SSN_s_sn_corr_neg=SSN_s_sn_corr_neg(PitchModIndex, :);
    SSN_n_sn_corr_neg= 1-[mr_corr_Data(curSNRinds).SSNuncorr_sn_n_neg]; %
    SSN_n_sn_corr_neg=SSN_n_sn_corr_neg(PitchModIndex, :);
    SSN_sn_ratio_neg= SSN_s_sn_corr_neg./SSN_n_sn_corr_neg;
    plot(curSNR*ones(1, length(curSNRinds))+2*xShift, SSN_sn_ratio_neg, 'rd');
    
    FLN_sn_ratio_neg_mean(snrVar)= mean(FLN_sn_ratio_neg);
    SSN_sn_ratio_neg_mean(snrVar)= mean(SSN_sn_ratio_neg);
end


subplot(121);
set(gca, 'xtick', snrs);
plot(-2*xShift+snrs, FLN_sn_ratio_pos_mean, 'b-o', 'linew', lw3, 'markersize', mrkSize);
plot(2*xShift+snrs, SSN_sn_ratio_pos_mean, 'r-o', 'linew', lw3, 'markersize', mrkSize);
grid on;
xlabel('SNR (dB)');
ylabel('$\overline{R_{f0,BP}(S).R_{f0,BP}(SN)/R_{f0,BP}(S).R_{f0,BP}(N)}$', 'interpreter', 'latex');
title(sprintf('ModFreq=%.0f Hz (corr(SN,S)/corr(SN,N)), +ve uR', modFreq));
set(gca, 'fontsize', fSize);


subplot(122);
set(gca, 'xtick', snrs);
plot(-2*xShift+snrs, FLN_sn_ratio_neg_mean, 'b-o', 'linew', lw3, 'markersize', mrkSize);
plot(2*xShift+snrs, SSN_sn_ratio_neg_mean, 'r-o', 'linew', lw3, 'markersize', mrkSize);
grid on;
xlabel('SNR (dB)');
ylabel('$\overline{R_{f0,BP}(S).R_{f0,BP}(SN)/R_{f0,BP}(S).R_{f0,BP}(N)}$', 'interpreter', 'latex');
title(sprintf('ModFreq=%.0f Hz (corr(SN,S)/corr(SN,N)), -ve uR', modFreq));
set(gca, 'fontsize', fSize);


% LG2=legend('SSN', 'SSN-NF', 'SSN-dPrime', 'FLN', 'FLN-NF', 'FLN-dPrime', 'location', 'northwest');
% LG2.FontSize=leg_fSize;
% set(gca, 'fontsize', fSize);

set(hanOneModFreq_Corr_ratio, 'units', 'inches', 'position', [1 1 14 5]);
fName_singleModFrequnCorr=sprintf('sEPSM_uR_f0_BP_modFreq%.0fHz_ratio_corr_%s_tRes%.0fms', modFreq, postfix_fName, TimeResolution*1e3);

if saveAllFigs
    saveas(hanOneModFreq_Corr_ratio, [saving_Dir fName_singleModFrequnCorr], 'tiff');
    saveas(hanOneModFreq_Corr_ratio, [Latex_Dir fName_singleModFrequnCorr], 'epsc');
end