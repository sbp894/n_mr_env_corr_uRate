clear;
clc;

chinIDs.NH= [321 322 325 338 341 343 346 347 354 355];
chinIDs.HI= [358 360 361 362]; %[358 360 361 362];

saveFigs= 1;

outfigDir= '/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/outFig/dPrime_intrinsic_plots/';
Root_loading_Dir= '/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/mr_corr_OUTPUT/BP_data_nCycles/';
modFitlerType= {'BP'};

if ~isdir(outfigDir)
    mkdir(outfigDir);
end

allModFreq=[4 8 16 32 64 128];
all_numCycles= [2 3 4 5 10];
fSize= 12;
lw= 2;
mrkSize= 12;
SNRs2use= [-10 -5 0];

all_numCycles_str= cellfun(@(x) num2str(x),num2cell(all_numCycles), 'uniformoutput', false);

for modVar=1:length(allModFreq)
    modFreq=allModFreq(modVar);
    figure(modVar);
    clf;
    
    for filtVar=1:length(modFitlerType)
        curModFiltType= modFitlerType{filtVar};
        
        cur_win_filttype_data= nan(2, length(SNRs2use), length(all_numCycles)); % 2 because scatter plot of SN&S vs SN&N corr.
        check_significance= nan(length(SNRs2use), length(all_numCycles)); 
        figure(modVar);
        subplot(length(modFitlerType), length(SNRs2use), (filtVar-1)*length(SNRs2use)+1);
        ylabel(['$' curModFiltType ':(\mu _{PTS} - \mu _{NH})/var_{norm}$'], 'interpreter', 'latex');
        
        for snrVar=1:length(SNRs2use)
            curSNR= SNRs2use(snrVar);
            
            for cycVar= 1:length(all_numCycles)
                curNCycles=all_numCycles(cycVar);
                cur_tRes= 1/modFreq*curNCycles;
                
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
                    
                    
                    nhINDs= ismember([mr_corr_Data.chinID], chinIDs.NH);
                    hiINDs= ismember([mr_corr_Data.chinID], chinIDs.HI);
                    snrINDs= [mr_corr_Data.SNR]== curSNR;
                    validINDs.NH= nhINDs & snrINDs;
                    validINDs.HI= hiINDs & snrINDs;
                    
                    %%
                    % calc d prime
                    dPrime.corr_S_SN.pos= calc_dPrime([mr_corr_Data(validINDs.HI).SSNcorr_s_sn_pos_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_pos_Final]);
                    dPrime.corr_S_SN.neg= calc_dPrime([mr_corr_Data(validINDs.HI).SSNcorr_s_sn_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_neg_Final]);
                    dPrime.corr_N_SN.pos= calc_dPrime(1-[mr_corr_Data(validINDs.HI).SSNuncorr_sn_n_pos_Final], 1-[mr_corr_Data(validINDs.NH).SSNuncorr_sn_n_pos_Final]); % subtracting from 1 is unnecessary
                    dPrime.corr_N_SN.neg= calc_dPrime(1-[mr_corr_Data(validINDs.HI).SSNuncorr_sn_n_neg_Final], 1-[mr_corr_Data(validINDs.NH).SSNuncorr_sn_n_neg_Final]); % subtracting from 1 is unnecessary
                    
                    dPrime.corr_S_SN.avg= (dPrime.corr_S_SN.pos+dPrime.corr_S_SN.neg)/2;
                    dPrime.corr_N_SN.avg= (dPrime.corr_N_SN.pos+dPrime.corr_N_SN.neg)/2;
                    
                    
                    [~,pVal]=ttest( ...
                        [ [mr_corr_Data(validINDs.NH).SSNcorr_s_n_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_n_pos_Final]], ...
                        [ [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_pos_Final] ]);
                    
                    cur_win_filttype_data(1, snrVar, cycVar)= dPrime.corr_S_SN.avg;
                    cur_win_filttype_data(2, snrVar, cycVar)= dPrime.corr_N_SN.avg;
                    check_significance(snrVar, cycVar)= pVal;
                    
                end
                
            end
            figure(modVar);
            subplot(length(modFitlerType), length(SNRs2use), (filtVar-1)*length(SNRs2use)+snrVar);
            grid on;
            if filtVar==1
                title(sprintf('SNR=%.0f', SNRs2use(snrVar)));
            end
            
            inds2use_this_modFreq= 1/modFreq*all_numCycles <=1;
            
            s_sn_ydata= squeeze(cur_win_filttype_data(1, snrVar, :)); % xdata = allModFreq
            n_sn_ydata= squeeze(cur_win_filttype_data(2, snrVar, :));
            sign_inds_pVals= check_significance(snrVar, :)<.05;
            hold on;
            plot(1:length(all_numCycles), s_sn_ydata, '-o', 'linew', lw, 'markersize' , mrkSize); % xdata = allModFreq
            plot(1:length(all_numCycles), n_sn_ydata, '-d', 'linew', lw, 'markersize' , mrkSize);
            plot(find(sign_inds_pVals), 0*find(sign_inds_pVals), '*k', 'linew', lw, 'markersize' , mrkSize+4);
            set(gca, 'xtick', 1:length(allModFreq(inds2use_this_modFreq)), 'xticklabel', all_numCycles_str(inds2use_this_modFreq), 'fontsize', fSize);
        end
    end
    subplot(length(modFitlerType), length(SNRs2use), length(modFitlerType)*length(SNRs2use)-1);
    xlabel('# cycles');
    
    subplot(length(modFitlerType), length(SNRs2use), 1);
    text(.05, .9, sprintf('fMod=%.0f Hz', modFreq), 'units', 'normalized');
    
    subplot(length(modFitlerType), length(SNRs2use), length(modFitlerType)*length(SNRs2use));
    legend('S&SN', 'N&SN', 'location', 'northwest');
    
    if saveFigs
        set(gcf, 'units', 'inches', 'position', [1 1 12 8]);
        figName2Save=sprintf('%sdPrimePlot_fMod%.0fHz', outfigDir, modFreq);
        saveas(gcf, figName2Save, 'tiff');
    end
end