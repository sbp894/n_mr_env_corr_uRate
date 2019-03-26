clear;
clc;

FreqCutOff= 4e3;
% FreqCutOff= 5e3;
% FreqCutOff= 12e3;

chinIDs.NH= [321 322 325 338 341 343 346 347 354 355];
chinIDs.HI= [358 360 361 362]; %[358 360 361 362];

saveFigs= 1;
debug_Significance= 1;
yy_var_lim= [.3 .2 .2 .1 .1 .1];


LatexDir= '/home/parida/Dropbox/Articles/Speech_in_noise_AN_NH_HI/figures/';
outfigDir= '/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/outFig/dPrime_intrinsic_plots/';
Root_loading_Dir= '/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/mr_corr_OUTPUT/BP_data_nCycles/';
modFitlerType= {'BP'};

nSProws= length(modFitlerType);

if ~isfolder(outfigDir)
    mkdir(outfigDir);
end

allModFreq=[4 8 16 32 64 128];
all_numCycles= [2 3 4 5 10];
fSize= 12;
lw= 2;
mrkSize= 12;
SNRs2use= [-10 -5 0];

all_numCycles_str= cellfun(@(x) num2str(x),num2cell(all_numCycles), 'uniformoutput', false);
allModFreq_str= cellfun(@(x) num2str(x),num2cell(allModFreq), 'uniformoutput', false);

final_dPrime_data= struct('corr_sn_s', nan(length(allModFreq), length(SNRs2use)), 'corr_sn_n', nan(length(allModFreq), length(SNRs2use)));

for modVar=1:length(allModFreq)
    
    modFreq=allModFreq(modVar);
    figure(modVar);
    clf;
    ax= nan(2*length(SNRs2use), 1);
    
    for filtVar=1:length(modFitlerType)
        curModFiltType= modFitlerType{filtVar};
        
        cur_win_filttype_data= nan(2, length(SNRs2use), length(all_numCycles)); % 2 because scatter plot of SN&S vs SN&N corr.
        check_significance= nan(length(SNRs2use), length(all_numCycles));
        if debug_Significance
            mu_s_sn_significance= nan(length(SNRs2use), length(all_numCycles));
            mu_s_n_significance= nan(length(SNRs2use), length(all_numCycles));
            std_significance= nan(length(SNRs2use), length(all_numCycles));
        end
        
        figure(modVar);
        subplot(nSProws, length(SNRs2use), 1);
        ylabel(['$' curModFiltType ':(\mu _{HI} - \mu _{NH})/std_{norm}$'], 'interpreter', 'latex');
        
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
                    
                    highCFinds= [mr_corr_Data.CF_Hz]>FreqCutOff;
                    mr_corr_Data= mr_corr_Data(~highCFinds);
                    
                    modFreq= temp.modFreq;
                    clear temp;
                    
                    
                    nhINDs= ismember([mr_corr_Data.chinID], chinIDs.NH);
                    hiINDs= ismember([mr_corr_Data.chinID], chinIDs.HI);
                    snrINDs= [mr_corr_Data.SNR]== curSNR;
                    validINDs.NH= nhINDs & snrINDs;
                    validINDs.HI= hiINDs & snrINDs;
                    
                    %%
                    % calc d prime
                    dPrime.corr_S_SN.pos= calc_dPrime([mr_corr_Data(validINDs.HI).FLNcorr_s_sn_pos_Final], [mr_corr_Data(validINDs.NH).FLNcorr_s_sn_pos_Final]);
                    dPrime.corr_S_SN.neg= calc_dPrime([mr_corr_Data(validINDs.HI).FLNcorr_s_sn_neg_Final], [mr_corr_Data(validINDs.NH).FLNcorr_s_sn_neg_Final]);
                    dPrime.corr_N_SN.pos= calc_dPrime(1-[mr_corr_Data(validINDs.HI).FLNuncorr_sn_n_pos_Final], 1-[mr_corr_Data(validINDs.NH).FLNuncorr_sn_n_pos_Final]); % subtracting from 1 is unnecessary
                    dPrime.corr_N_SN.neg= calc_dPrime(1-[mr_corr_Data(validINDs.HI).FLNuncorr_sn_n_neg_Final], 1-[mr_corr_Data(validINDs.NH).FLNuncorr_sn_n_neg_Final]); % subtracting from 1 is unnecessary
                    
                    dPrime.corr_S_SN.avg= (dPrime.corr_S_SN.pos+dPrime.corr_S_SN.neg)/2;
                    dPrime.corr_N_SN.avg= (dPrime.corr_N_SN.pos+dPrime.corr_N_SN.neg)/2;
                    
                    
                    [~,pVal]=ttest( ...
                        [ [mr_corr_Data(validINDs.NH).FLNcorr_s_n_neg_Final], [mr_corr_Data(validINDs.NH).FLNcorr_s_n_pos_Final] ], ...
                        [ [mr_corr_Data(validINDs.NH).FLNcorr_s_sn_neg_Final], [mr_corr_Data(validINDs.NH).FLNcorr_s_sn_pos_Final] ]);
                    
                    if debug_Significance
                        mu_s_sn_significance(snrVar, cycVar)= nanmean([ [mr_corr_Data(validINDs.NH).FLNcorr_s_sn_neg_Final], [mr_corr_Data(validINDs.NH).FLNcorr_s_sn_pos_Final] ]);
                        mu_s_n_significance(snrVar, cycVar)= nanmean([ [mr_corr_Data(validINDs.NH).FLNcorr_s_n_neg_Final], [mr_corr_Data(validINDs.NH).FLNcorr_s_n_pos_Final] ]);
                        std_significance(snrVar, cycVar)= sqrt(nanstd([ [mr_corr_Data(validINDs.NH).FLNcorr_s_sn_neg_Final], [mr_corr_Data(validINDs.NH).FLNcorr_s_sn_pos_Final] ])...
                            *nanstd([ [mr_corr_Data(validINDs.NH).FLNcorr_s_sn_neg_Final], [mr_corr_Data(validINDs.NH).FLNcorr_s_sn_pos_Final] ]));
                    end
                    
                    cur_win_filttype_data(1, snrVar, cycVar)= dPrime.corr_S_SN.avg;
                    cur_win_filttype_data(2, snrVar, cycVar)= dPrime.corr_N_SN.avg;
                    check_significance(snrVar, cycVar)= pVal;
                    
                end
                if curNCycles==3
                    fprintf('Entering \n');
                    final_dPrime_data.corr_sn_s(modVar, snrVar)= dPrime.corr_S_SN.avg;
                    final_dPrime_data.corr_sn_n(modVar, snrVar)= dPrime.corr_N_SN.avg;
                end
                
            end
            figure(modVar);
            subplot(nSProws, length(SNRs2use), snrVar);
            yyaxis right;
            plot(1:length(all_numCycles), check_significance(snrVar, :), '-ok', 'linew', lw, 'markersize' , mrkSize);
            ylim([0 .2]);
            
            yyaxis left;
            ax(snrVar)= gca;
            grid on;
            if filtVar==1
                title(sprintf('FLN: SNR=%.0f', SNRs2use(snrVar)));
            end
            
            inds2use_this_modFreq= 1/modFreq*all_numCycles <=1;
            
            s_sn_ydata= squeeze(cur_win_filttype_data(1, snrVar, :)); % xdata = allModFreq
            n_sn_ydata= squeeze(cur_win_filttype_data(2, snrVar, :));
            sign_inds_pVals= check_significance(snrVar, :)<.05;
            hold on;
            plot(1:length(all_numCycles), s_sn_ydata, '-ob', 'linew', lw, 'markersize' , mrkSize); % xdata = allModFreq
            plot(1:length(all_numCycles), n_sn_ydata, '-dr', 'linew', lw, 'markersize' , mrkSize);
            plot(find(sign_inds_pVals), 0*find(sign_inds_pVals), '*k', 'linew', lw, 'markersize' , mrkSize+4);
            set(gca, 'xtick', 1:length(allModFreq(inds2use_this_modFreq)), 'xticklabel', all_numCycles_str(inds2use_this_modFreq), 'fontsize', fSize);
            
%             subplot(nSProws, length(SNRs2use), length(SNRs2use)*nSProws+snrVar);
            
            yyaxis right;
            plot(1:length(all_numCycles), std_significance(snrVar, :), '-dk', 'linew', lw, 'markersize' , mrkSize);
            set(gca, 'ycolor', 'k');
            set(gca, 'xtick', 1:length(allModFreq(inds2use_this_modFreq)), 'xticklabel', all_numCycles_str(inds2use_this_modFreq), 'fontsize', fSize);
            grid on;
            ylim([.2*yy_var_lim(modVar) yy_var_lim(modVar)])
            
%             yyaxis left;
%             ax(snrVar+length(SNRs2use))= gca;
%             hold on;
%             plot(1:length(all_numCycles), mu_s_sn_significance(snrVar, :), '-ob', 'linew', lw, 'markersize' , mrkSize);
%             plot(1:length(all_numCycles), mu_s_n_significance(snrVar, :), '-or', 'linew', lw, 'markersize' , mrkSize);
%             ylim([0 .4]);
        end
    end
    subplot(nSProws, length(SNRs2use), length(modFitlerType)*length(SNRs2use)-1);
    xlabel('# cycles');
    
    subplot(nSProws, length(SNRs2use), 1);
    text(.05, .9, sprintf('fMod=%.0f Hz', modFreq), 'units', 'normalized');
    
    subplot(nSProws, length(SNRs2use), length(SNRs2use));
    legend('S&SN', 'N&SN', 'location', 'best');
    linkaxes(ax, 'x');
    linkaxes(ax(1:3), 'y');
    
    subplot(nSProws, length(SNRs2use), 1);
    ylim([-1.5 1.5]);
    xlim([find(inds2use_this_modFreq, 1)-.5 find(inds2use_this_modFreq, 1, 'last')+.5])
    
    if saveFigs
        set(gcf, 'units', 'inches', 'position', [1 1 12 8]);
        figName2Save=sprintf('%sdPrimePlot_fMod%.0fHz_fMax%.0fkHz_FLN', outfigDir, modFreq, FreqCutOff/1e3);
        saveas(gcf, figName2Save, 'tiff');
    end
end

close all;
%%
mrkSize2= 18;
figure(1);
clf;
co= get(gca, 'colororder');
col.nh= co(4,:);
col.hi= co(5,:);
fSize= 20;
lw3= 3;
modFreqs2use= [8 16 32 64 128];
inds2use= ismember(allModFreq, modFreqs2use);
dx= nan(size(final_dPrime_data.corr_sn_n,2), 1);
for snrVar= 1:size(final_dPrime_data.corr_sn_n,2)
    dx(snrVar)= subplot(1, size(final_dPrime_data.corr_sn_n,2), snrVar);
    hold on;
    plot(find(inds2use), final_dPrime_data.corr_sn_s(inds2use, snrVar), '-o', 'color', col.nh, 'linew', lw3, 'markersize', mrkSize2);
    plot(find(inds2use), final_dPrime_data.corr_sn_n(inds2use, snrVar), '-d', 'color', col.hi, 'linew', lw3, 'markersize', mrkSize2);
    plot(find(inds2use), 0*final_dPrime_data.corr_sn_n(inds2use, snrVar), '-k', 'linew', lw3, 'markersize', mrkSize2);
    set(gca, 'xtick', find(inds2use), 'xticklabel', allModFreq_str(inds2use), 'fontsize', fSize);
    grid on;
    ylim([-1.75 1.75]);
    if snrVar==2
        xlabel('Freq_m_o_d (Hz)');
    end
    if snrVar==1
        ylabel('$d''=(\mu _{HI} - \mu _{NH})/ \sigma_{norm}$', 'interpreter', 'latex');
    end
    title(sprintf('SNR=%.0f dB', SNRs2use(snrVar)));
end
linkaxes(dx);
xlim([find(inds2use, 1)-.5 find(inds2use, 1, 'last')+.5]);
legend({'$\emph{corr}(SN,S)$', '$\emph{corr}(SN,N)$'}, 'interpreter', 'latex', 'location', 'southeast');
set(gcf, 'units', 'inches', 'position', [1 1 14 6]);
% 

saveas(gcf, [LatexDir 'AN_intrinsic_dPrime_FLN'], 'epsc');

saveas(gcf, [outfigDir 'AN_intrinsic_dPrime_FLN'], 'png');