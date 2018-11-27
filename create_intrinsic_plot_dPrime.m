clear;
clc;

chinIDs.NH= [321 322 325 338 341 343 346 347 354 355];
chinIDs.HI= [358 360 361 362]; %[358 360 361 362];

saveFigs= 1;

outfigDir= '/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/outFig/dPrime_intrinsic_plots/';
loading_Dir='/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/mr_corr_OUTPUT/BP_data/';
% modFitlerType= {'LP', 'BP'};
modFitlerType= {'BP'};

if ~isdir(outfigDir)
    mkdir(outfigDir);
end

allModFreq=[4 8 16 32 64 128];
allWindows= [32 64 128]*1e-3;
fSize= 12;
lw= 2;
mrkSize= 12;
SNRs2use= [-10 -5 0];

allModFreq_str= cellfun(@(x) num2str(x),num2cell(allModFreq), 'uniformoutput', false);

for winVar= 1:length(allWindows)
    TimeResolution=allWindows(winVar);
    
    figure(winVar);
    clf;
    
    for filtVar=1:length(modFitlerType)
        curModFiltType= modFitlerType{filtVar};
        
        cur_win_filttype_data= nan(2, length(SNRs2use), length(allModFreq)); % 2 because scatter plot of SN&S vs SN&N corr.
        figure(winVar);
        subplot(length(modFitlerType), length(SNRs2use), (filtVar-1)*length(SNRs2use)+1);
        ylabel(['$' curModFiltType ':(\mu _{PTS} - \mu _{NH})/var_{norm}$'], 'interpreter', 'latex');
        
        for snrVar=1:length(SNRs2use)
            curSNR= SNRs2use(snrVar);
            
            for modVar=1:length(allModFreq)
                modFreq=allModFreq(modVar);
                
                fName2Load=sprintf('%suRate_uR_f0_%s_%d_win_%.0fms.mat', loading_Dir, curModFiltType, modFreq, TimeResolution*1e3);
                
                %                 ax= nan(length(SNRs2use), 1);
                
                temp= load(fName2Load);
                mr_corr_Data= temp.mr_corr_Data;
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
                
                cur_win_filttype_data(1, snrVar, modVar)= dPrime.corr_S_SN.avg;
                cur_win_filttype_data(2, snrVar, modVar)= dPrime.corr_N_SN.avg;
                
            end
            figure(winVar);
            subplot(length(modFitlerType), length(SNRs2use), (filtVar-1)*length(SNRs2use)+snrVar);
            grid on;
            if filtVar==1
                title(sprintf('SNR=%.0f', SNRs2use(snrVar)));
            end
            
            s_sn_ydata= squeeze(cur_win_filttype_data(1, snrVar, :)); % xdata = allModFreq
            n_sn_ydata= squeeze(cur_win_filttype_data(2, snrVar, :));
            hold on;
            plot(1:length(allModFreq), s_sn_ydata, '-o', 'linew', lw, 'markersize' , mrkSize); % xdata = allModFreq
            plot(1:length(allModFreq), n_sn_ydata, '-d', 'linew', lw, 'markersize' , mrkSize);
            set(gca, 'xtick', 1:length(allModFreq), 'xticklabel', allModFreq_str, 'fontsize', fSize);
        end
    end
    subplot(length(modFitlerType), length(SNRs2use), length(modFitlerType)*length(SNRs2use)-1);
    xlabel('Mod Freq (cutoff in Hz)');
    
    subplot(length(modFitlerType), length(SNRs2use), 1);
    text(.05, .9, sprintf('tRes=%.0f ms', TimeResolution*1e3), 'units', 'normalized');
    
    subplot(length(modFitlerType), length(SNRs2use), length(modFitlerType)*length(SNRs2use));
    legend('S&SN', 'N&SN', 'location', 'northwest');
    
    if saveFigs
        set(gcf, 'units', 'inches', 'position', [1 1 12 8]);
        figName2Save=sprintf('%sdPrimePlot_%.0f', outfigDir, TimeResolution*1e3);
        saveas(gcf, figName2Save, 'tiff');
    end
end