% Only for steady-state noise for now
% for various modulation frequency and time window durations, plots
% corr(SN,S), corr(SN, N) and corr(S, N). Has two subplots: one for normal
% hearing and one for hearing impaired.
% Have to run generate_bandpass_data.m to generate first stage of data. Bad
% idea. Make it independnent?

clear;
clc;

chinIDs.NH= [321 322 325 338 341 343 346 347 354 355];
chinIDs.HI= [361 362];

saveFigs= 1;
cache_figHandle=1111;
std_dev1_err0=1;
if std_dev1_err0
    err_post_fix= 'std';
else
    err_post_fix= 'se';
end

hypo_result_syms= '-*';

loading_Dir='/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/mr_corr_OUTPUT/BP_data/';
outfigDir.parent= '/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/outFig/intrinsic_plots/corr_ofs_ssnn/';
outfigDir.tif= [outfigDir.parent 'tifs' filesep];
if ~isdir(outfigDir.tif)
    mkdir(outfigDir.tif);
end
outfigDir.fig= [outfigDir.parent 'figs' filesep];
if ~isdir(outfigDir.fig)
    mkdir(outfigDir.fig);
end

modFitlerType= {'BP'};

allModFreq=[32 64 128];
allWindows= [32 64]*1e-3;
fSize= 12;
lw=2;
SNRs2use= [-10 -5 0];

figure(cache_figHandle);
clf;
co= get(gca, 'colororder');
colScheme.s_sn_corr= co(1,:); % blue
colScheme.n_sn_corr= co(2,:); % red
colScheme.s_n_corr= co(5,:); % green

condMat= [...
    16 128 1; ...
    32 64 1; ...
    64 32 1; ...
    128 32 1;]; % freqMod | timeWin | filtType

% for winVar= 1:length(allWindows)
%     TimeResolution=allWindows(winVar);
%
%     for filtVar=1:length(modFitlerType)
%         curModFiltType= modFitlerType{filtVar};
%
%
%         for modVar=1:length(allModFreq)
%             modFreq=allModFreq(modVar);
%

for condVar=1:size(condMat, 1)
    modFreq=condMat(condVar, 1);
    TimeResolution=condMat(condVar, 2)/1e3;
    curModFiltType=modFitlerType{condMat(condVar, 3)};
    
    
    fName2Load=sprintf('%suRate_uR_f0_%s_%d_win_%.0fms.mat', loading_Dir, curModFiltType, modFreq, TimeResolution*1e3);
    figName2Save.tif=sprintf('%scorr_ofs_ssnn_all_snrs_%s_%d_win_%.0fms_%s', outfigDir.tif, curModFiltType, modFreq, TimeResolution*1e3, err_post_fix);
    figName2Save.fig=sprintf('%scorr_ofs_ssnn_all_snrs_%s_%d_win_%.0fms_%s', outfigDir.fig, curModFiltType, modFreq, TimeResolution*1e3, err_post_fix);
    
    
    temp= load(fName2Load);
    mr_corr_Data= temp.mr_corr_Data;
    modFreq= temp.modFreq;
    clear temp;
    
    corr_data_for_cur_win_mod_filt.s_sn= repmat(struct('nh', nan, 'hi', nan), 1, length(SNRs2use));
    corr_data_for_cur_win_mod_filt.s_n= repmat(struct('nh', nan, 'hi', nan), 1, length(SNRs2use));
    corr_data_for_cur_win_mod_filt.n_sn= repmat(struct('nh', nan, 'hi', nan), 1, length(SNRs2use));
    innerStruct= struct('hypoVal', nan', 'pVal', nan);
    stat.corr_sn_s= repmat(struct('nh', innerStruct, 'hi', innerStruct), 1, length(SNRs2use));
    stat.corr_sn_n= repmat(struct('nh', innerStruct, 'hi', innerStruct), 1, length(SNRs2use));
    
    for snrVar=1:length(SNRs2use)
        curSNR= SNRs2use(snrVar);
        
        nhINDs= ismember([mr_corr_Data.chinID], chinIDs.NH);
        hiINDs= ismember([mr_corr_Data.chinID], chinIDs.HI);
        snrINDs= [mr_corr_Data.SNR]== curSNR;
        validINDs.NH= nhINDs & snrINDs;
        validINDs.HI= hiINDs & snrINDs;
        
        
        corr_data_for_cur_win_mod_filt.n_sn(snrVar).nh= 1- [[mr_corr_Data(validINDs.NH).SSNuncorr_sn_n_pos_Final], [mr_corr_Data(validINDs.NH).SSNuncorr_sn_n_neg_Final]];
        corr_data_for_cur_win_mod_filt.n_sn(snrVar).hi= 1- [[mr_corr_Data(validINDs.HI).SSNuncorr_sn_n_pos_Final], [mr_corr_Data(validINDs.HI).SSNuncorr_sn_n_neg_Final]];
        
        corr_data_for_cur_win_mod_filt.s_sn(snrVar).nh= [[mr_corr_Data(validINDs.NH).SSNcorr_s_sn_pos_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_neg_Final]];
        corr_data_for_cur_win_mod_filt.s_sn(snrVar).hi= [[mr_corr_Data(validINDs.HI).SSNcorr_s_sn_pos_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_sn_neg_Final]];
        
        corr_data_for_cur_win_mod_filt.s_n(snrVar).nh= [[mr_corr_Data(validINDs.NH).SSNcorr_s_n_pos_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_n_neg_Final]];
        corr_data_for_cur_win_mod_filt.s_n(snrVar).hi= [[mr_corr_Data(validINDs.HI).SSNcorr_s_n_pos_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_n_neg_Final]];
        
        % Add [**] vs [-*] vs [*-] for NF [first * for NH, second * for PTS]
        [stat.corr_sn_s(snrVar).nh.hypoVal, stat.corr_sn_s(snrVar).nh.pVal]= ttest2(corr_data_for_cur_win_mod_filt.s_sn(snrVar).nh, corr_data_for_cur_win_mod_filt.s_n(snrVar).nh);
        [stat.corr_sn_s(snrVar).hi.hypoVal, stat.corr_sn_s(snrVar).hi.pVal]= ttest2(corr_data_for_cur_win_mod_filt.s_sn(snrVar).hi, corr_data_for_cur_win_mod_filt.s_n(snrVar).hi);
        [stat.corr_sn_n(snrVar).nh.hypoVal, stat.corr_sn_n(snrVar).nh.pVal]= ttest2(corr_data_for_cur_win_mod_filt.n_sn(snrVar).nh, corr_data_for_cur_win_mod_filt.s_n(snrVar).nh);
        [stat.corr_sn_n(snrVar).hi.hypoVal, stat.corr_sn_n(snrVar).hi.pVal]= ttest2(corr_data_for_cur_win_mod_filt.n_sn(snrVar).hi, corr_data_for_cur_win_mod_filt.s_n(snrVar).hi);
        
        fprintf('tRes=%.0f ms, filtType= %s, fMod=%.0f, snr=%.0f dB, SN[%s%s]N d''= %.2f, SN[%s%s]S d''=%.2f, where stat[NH PTS]\n', TimeResolution*1e3, curModFiltType, ...
            modFreq, curSNR, hypo_result_syms(stat.corr_sn_n(snrVar).nh.hypoVal+1),  hypo_result_syms(stat.corr_sn_n(snrVar).hi.hypoVal+1), ...
            calc_dPrime(corr_data_for_cur_win_mod_filt.n_sn(snrVar).hi, corr_data_for_cur_win_mod_filt.n_sn(snrVar).nh), ...
            hypo_result_syms(stat.corr_sn_s(snrVar).nh.hypoVal+1),  hypo_result_syms(stat.corr_sn_s(snrVar).hi.hypoVal+1), ...
            calc_dPrime(corr_data_for_cur_win_mod_filt.s_sn(snrVar).hi, corr_data_for_cur_win_mod_filt.s_sn(snrVar).nh));
    end
    
    fprintf('--------------------------------------------------------------------------------\n--------------------------------------------------------------------------------\n');
    figure(cache_figHandle);
    clf;
    spAx= nan(2,1);
    % Normal hearing
    spAx(1)=subplot(1, 2, 1);
    hold on;
    ax(1)=errorbar_plot_fun(SNRs2use, {corr_data_for_cur_win_mod_filt.n_sn.nh}, std_dev1_err0, colScheme.n_sn_corr,lw);
    ax(2)=errorbar_plot_fun(SNRs2use, {corr_data_for_cur_win_mod_filt.s_sn.nh}, std_dev1_err0, colScheme.s_sn_corr,lw);
    ax(3)=errorbar_plot_fun(SNRs2use, {corr_data_for_cur_win_mod_filt.s_n.nh}, std_dev1_err0, colScheme.s_n_corr,lw);
    
    set(gca, 'xtick', SNRs2use, 'fontsize', fSize);
    legend(ax, 'N & SN','S & SN','S & N');
    grid on;
    title('Normal Hearing');
    xlim([-11 1]);
    ylabel('Corr (Only SSN)');
    xlabel('SNR (dB)');
    
    % Hearing Impaired
    spAx(2)=subplot(1, 2, 2);
    hold on;
    ax(1)=errorbar_plot_fun(SNRs2use, {corr_data_for_cur_win_mod_filt.n_sn.hi}, std_dev1_err0, colScheme.n_sn_corr,lw);
    ax(2)=errorbar_plot_fun(SNRs2use, {corr_data_for_cur_win_mod_filt.s_sn.hi}, std_dev1_err0, colScheme.s_sn_corr,lw);
    ax(3)=errorbar_plot_fun(SNRs2use, {corr_data_for_cur_win_mod_filt.s_n.hi}, std_dev1_err0, colScheme.s_n_corr,lw);
    
    set(gca, 'xtick', SNRs2use, 'fontsize', fSize);
    legend(ax, 'N & SN','S & SN','S & N');
    grid on;
    title('Hearing Impaired');
    xlim([-11 1]);
    tx= text(0.01, .01, sprintf('$%s |  F_{co}=%.0f Hz | tRes=%.0fms$', curModFiltType, modFreq, TimeResolution*1e3), 'interpreter', 'latex', 'units', 'normalized');
    ylabel('Corr (Only SSN)');
    xlabel('SNR (dB)');
    
    linkaxes(spAx, 'y');
    
    if saveFigs
        figure(cache_figHandle);
        set(cache_figHandle, 'units', 'inches', 'position', [1 1 14 4]);
        saveas(gcf, figName2Save.tif, 'tiff');
        saveas(gcf, figName2Save.fig);
    end
end
%     end
% end
% end

function ax= errorbar_plot_fun(x_data, y_data, std_dev1_err0, plot_col, lw)
if std_dev1_err0
    ax=errorbar(x_data, cellfun(@(x) nanmean(x), y_data), cellfun(@(x) nanstd(x), y_data), 'color', plot_col, 'linew', lw);
else
    ax=errorbar(x_data, cellfun(@(x) nanmean(x), y_data), cellfun(@(x) nanstd(x)/sqrt(numel(x)), y_data), 'color', plot_col, 'linew', lw);
end
end