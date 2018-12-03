clear;
clc;

chinIDs.NH= [321 322 325 338 341 343 346 347 354 355];
chinIDs.HI= [361 362];

saveFigs= 0;
cache_figHandle=1111;

outfigDir= '/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/outFig/intrinsic_plots/';
Root_loading_Dir='/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/mr_corr_OUTPUT/BP_data_nCycles/';
modFitlerType= {'BP'};

if ~isdir(outfigDir)
    mkdir(outfigDir);
end

allModFreq=[4 8 16 32 64 128];
all_numCycles= [2 3 4 5 10];

fSize= 12;
SNRs2use= [-10 -5 0];


% condition_all= get_conditions(allWindows, modFitlerType, allModFreq);


for modVar=3:4 %1:length(allModFreq)
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
                    subplot(length(SNRs2use), 2, (snrVar-1)*2+1);
                    hold on;
                    plot([mr_corr_Data(validINDs.NH).SSNcorr_s_n_pos_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_pos_Final], 'ob');
                    plot([mr_corr_Data(validINDs.NH).SSNcorr_s_n_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_neg_Final], 'ob');
                    plot([mr_corr_Data(validINDs.HI).SSNcorr_s_n_pos_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_sn_pos_Final], 'dr');
                    plot([mr_corr_Data(validINDs.HI).SSNcorr_s_n_neg_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_sn_neg_Final], 'dr');
                    [~,pVal]=ttest( ...
                        [ [mr_corr_Data(validINDs.NH).SSNcorr_s_n_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_n_pos_Final]], ...
                        [ [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_pos_Final] ]);
                    
                    tx= text(0.1, .6, sprintf('$F_{co}=%.0f Hz | nCyc=%.0f | pVal=%.4f$', modFreq, curNCycles, pVal), 'interpreter', 'latex');
                    
                    grid on;
                    xlabel('corr (S, N)');
                    ylabel('corr (S, SN)');
                    set(gca, 'fontsize',  fSize);
                    xlim([0 .7]);
                    ylim([0 .7]);
                    title(sprintf('SNR = %d dB', curSNR));
                    plot([0 1], [0 1], 'k', 'linew', 2);
                    
                    subplot(length(SNRs2use), 2, 2*snrVar);
                    hold on;
                    plot(1-[mr_corr_Data(validINDs.NH).SSNuncorr_sn_n_pos_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_pos_Final], 'ob');
                    plot(1-[mr_corr_Data(validINDs.NH).SSNuncorr_sn_n_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_neg_Final], 'ob');
                    plot(1-[mr_corr_Data(validINDs.HI).SSNuncorr_sn_n_pos_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_sn_pos_Final], 'dr');
                    plot(1-[mr_corr_Data(validINDs.HI).SSNuncorr_sn_n_neg_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_sn_neg_Final], 'dr');
                    
                    grid on;
                    xlabel('corr (SN, N)');
                    ylabel('corr (SN, S)');
                    set(gca, 'fontsize',  fSize);
                    xlim([0 .7]);
                    ylim([0 .7]);
                    plot([0 1], [0 1], 'k', 'linew', 2);
                    
                    tx= text(0.1, .6, sprintf('$F_{co}=%.0f Hz | ncyc=%.0f$', modFreq, curNCycles), 'interpreter', 'latex');
                    
                    title(sprintf('SNR = %d dB', curSNR));
                    grid on;
                end
                
                %%
                
                if saveFigs
                    figure(cache_figHandle);
                    set(gcf, 'units', 'inches', 'position', [1 1 14 4]);
                    figName2Save=sprintf('%suRate_uR_f0_%s_%d_win_%.0fms_2d', outfigDir, curModFiltType, modFreq, curNCycles*1e3);
                    saveas(gcf, figName2Save, 'tiff');
                    saveas(gcf, figName2Save);
                end
                clf;
            end
        end
    end
end