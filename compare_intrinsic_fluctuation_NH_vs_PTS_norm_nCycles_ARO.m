clear;
clc;

chinIDs.NH= [321 322 325 338 341 343 346 347 354 355];
chinIDs.HI= [361 362];

saveFigs= 1;
cache_figHandle=1111;

outfigDir= '/home/parida/Dropbox/Conferences/ARO-2019/latex/Figures/eps/';
Root_loading_Dir='/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/mr_corr_OUTPUT/BP_data_nCycles/';
modFitlerType= {'BP'};

if ~isdir(outfigDir)
    mkdir(outfigDir);
end

allModFreq= 8; %[4 8 16 32 64 128];
all_numCycles= 3; %[2 3 4 5 10];

SNRs2use= -5; %[-10 -5 0];

figure(cache_figHandle);
clf;
co= get(gca, 'colororder');
mrkCol.nh= co(1,:);
mrkCol.hi= co(2,:);
mrkSize= 12;
fSize= 20;
lw= 3;
colLims= [-.05 .75];
tickvals= [0 .2 .4 .6];
corrRes= .05;

% condition_all= get_conditions(allWindows, modFitlerType, allModFreq);


for modVar=1% :4 %1:length(allModFreq)mr_corr_Data
    modFreq=allModFreq(modVar);
    figure(cache_figHandle);
    clf;
    
    for filtVar=1:length(modFitlerType)
        curModFiltType= modFitlerType{filtVar};
        
        cur_win_filttype_data= nan(2, length(SNRs2use), length(allModFreq));
        
        for cycVar= 1:length(all_numCycles)
            curNCycles=all_numCycles(cycVar);
            
            allChins= [chinIDs.NH chinIDs.HI];
            mr_corr_Data= [];
            for loadVar= 1:length(allChins)
                curChinID=  allChins(loadVar);
                loading_Dir= [Root_loading_Dir 'Q' num2str(curChinID) filesep];
                
                fName2Load=sprintf('%sQ%d_uRate_uR_f0_%s_%d_nCyc_%.0f.mat', loading_Dir, curChinID, curModFiltType, modFreq, curNCycles);
                if exist(fName2Load, 'file')
                    temp= load(fName2Load);
                    mr_corr_Data= [mr_corr_Data; temp.mr_corr_Data];
                end
            end
            
            validInds= [mr_corr_Data.CF_Hz]>5e3;
            mr_corr_Data(validInds)= [];
            
            if ~isempty(mr_corr_Data)
                modFreq= temp.modFreq;
                clear temp;
                
                
                
                for snrVar=1:length(SNRs2use)
                    curSNR= SNRs2use(snrVar);
                    
                    nhINDs= ismember([mr_corr_Data.chinID], chinIDs.NH);
                    hiINDs= ismember([mr_corr_Data.chinID], chinIDs.HI);
                    snrINDs= [mr_corr_Data.SNR]== curSNR;
                    validINDs.NH= nhINDs & snrINDs;
                    validINDs.HI= hiINDs & snrINDs;
                    
                    
                    %%
                    figure(cache_figHandle);
                    %                     subplot(length(SNRs2use), 2, (snrVar-1)*2+1);
                    %                     hold on;
                    %                     plot([mr_corr_Data(validINDs.NH).SSNcorr_s_n_pos_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_pos_Final], 'ob');
                    %                     plot([mr_corr_Data(validINDs.NH).SSNcorr_s_n_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_neg_Final], 'ob');
                    %                     plot([mr_corr_Data(validINDs.HI).SSNcorr_s_n_pos_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_sn_pos_Final], 'dr');
                    %                     plot([mr_corr_Data(validINDs.HI).SSNcorr_s_n_neg_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_sn_neg_Final], 'dr');
                    %                     [~,pVal]=ttest( ...
                    %                         [ [mr_corr_Data(validINDs.NH).SSNcorr_s_n_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_n_pos_Final]], ...
                    %                         [ [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_pos_Final] ]);
                    %
                    %                     tx= text(0.1, .6, sprintf('$F_{co}=%.0f Hz | nCyc=%.0f | pVal=%.4f$', modFreq, curNCycles, pVal), 'interpreter', 'latex');
                    %
                    %                     grid on;
                    %                     xlabel('corr (S, N)');
                    %                     ylabel('corr (S, SN)');
                    %                     set(gca, 'fontsize',  fSize);
                    %                     xlim([0 .7]);
                    %                     ylim([0 .7]);
                    %                     title(sprintf('SNR = %d dB', curSNR));
                    %                     plot([0 1], [0 1], 'k', 'linew', 2);
                    
                    ax(1)= subplot(4, 4, [2 3 4 6 7 8 10 11 12]);
                    hold on;
                    
                    % NH
                    cur_sn_n_nh= [1-[mr_corr_Data(validINDs.NH).SSNuncorr_sn_n_pos_Final], 1-[mr_corr_Data(validINDs.NH).SSNuncorr_sn_n_neg_Final]];
                    cur_sn_s_nh= [[mr_corr_Data(validINDs.NH).SSNcorr_s_sn_pos_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_neg_Final]];
                    
                    % HI
                    cur_sn_n_hi= [1-[mr_corr_Data(validINDs.HI).SSNuncorr_sn_n_pos_Final], [mr_corr_Data(validINDs.HI).SSNuncorr_sn_n_neg_Final]];
                    cur_sn_s_hi= [[mr_corr_Data(validINDs.HI).SSNcorr_s_sn_pos_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_sn_neg_Final]];
                    
                    plHan(1)= plot(cur_sn_n_nh, cur_sn_s_nh, 's', 'color', mrkCol.nh, 'markersize', mrkSize, 'linew', lw);
                    plHan(2)= plot(cur_sn_n_hi, cur_sn_s_hi, 'p', 'color', mrkCol.hi, 'markersize', mrkSize, 'linew', lw);
                    
                    grid on;
                    set(gca, 'fontsize',  fSize, 'xticklabel', '', 'yticklabel', '');
                    xlim([-0.05 .7]);
                    ylim([-0.05 .7]);
                    plot([0 1], [0 1], 'k', 'linew', 2);
                    
                    tx= text(0.1, .6, sprintf('$f_{mod}=%.0f Hz | T_{win}=%.0fms$', modFreq, curNCycles/modFreq*1e3), 'interpreter', 'latex', 'fontsize', fSize);
                    
                    xlim(colLims);
                    title(sprintf('SNR = %d dB', curSNR));
                    grid on;
                    
                    corrBinEdges= min(colLims):corrRes:max(colLims);
                    corrBinCent=  (corrBinEdges(1:end-1)+corrBinEdges(2:end))/2;
                    
                    ax(2)= subplot(4, 4, [1 5 9]);
                    hold on;
                    [hist_sn_s_nh]= histcounts(cur_sn_s_nh, corrBinEdges, 'Normalization', 'pdf');
                    [hist_sn_s_hi]= histcounts(cur_sn_s_hi, corrBinEdges, 'Normalization', 'pdf');
                    bl= barh(corrBinCent, hist_sn_s_nh, 'b', 'facealpha', .5);
                    %                     set(bl, 'color', mrkCol.nh);
                    % , 'color', mrkCol.hi
                    barh(corrBinCent, hist_sn_s_hi, 'r', 'facealpha', .5);
                    ylabel('corr (SN, S)');
                    grid on;
                    set(gca, 'fontsize',  fSize, 'ytick', tickvals);
                    axis tight;
                    xlabel('PDF');
                    
                    ax(3)= subplot(4, 4, [14 15 16]);
                    hold on;
                    [hist_sn_n_nh]= histcounts(cur_sn_n_nh, corrBinEdges, 'Normalization', 'pdf');
                    [hist_sn_n_hi]= histcounts(cur_sn_n_hi, corrBinEdges, 'Normalization', 'pdf');
                    bar(corrBinCent, hist_sn_n_nh, 'b', 'facealpha', .5);
                    bar(corrBinCent, hist_sn_n_hi, 'r', 'facealpha', .5);
                    xlabel('corr (SN, N)');
                    grid on;
                    set(gca, 'fontsize',  fSize, 'xtick', tickvals);
                    axis tight;
                    ylabel('PDF');
                    
                    legend(plHan, 'NH', 'HI', 'location', 'northwest');
                    linkaxes(ax(1, 2), 'y');
                    linkaxes(ax(1, 3), 'x');
                    xlim(colLims);
                end
                
                %%
                
                set(gcf, 'units', 'inches', 'position', [1 1 10 8]);
                if saveFigs
                    figure(cache_figHandle);
                    figName2Save=sprintf('%suRate_uR_f0_%s_%d_win_nCyc%.0f', outfigDir, curModFiltType, modFreq, curNCycles);
                    %                     saveas(gcf, figName2Save, 'tiff');
                    saveas(gcf, figName2Save, 'epsc');
                end
                %                 clf;
            end
        end
    end
end