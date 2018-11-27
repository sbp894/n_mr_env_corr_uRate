% To run this code, need to 
% 1. convert mfiles to matfiles
% 2. fix bf track unit
% 3. screenDataMat
% 4. chin mr n sEPSM
% 5. generate_bandpass_data
% Have to remove 4 and directly use 5. Codes already there in
% LF_speech_analysis/load_danish_data
clear;
clc;

chinIDs.NH= [321 322 325 338 341 343 346 347 354 355];
chinIDs.HI= [358 360 361 362];

saveFigs= 0;
cache_figHandle=1111;

outfigDir= '/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/outFig/intrinsic_plots/';
loading_Dir='/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/mr_corr_OUTPUT/BP_data/';
modFitlerType= {'BP'};

if ~isdir(outfigDir)
    mkdir(outfigDir);
end

allModFreq=64; % [32 64 128];
allWindows= 64 *1e-3; %[32 64 128]*1e-3;
fSize= 12;
SNRs2use= 0 ;%[-10 -5 0];


% condition_all= get_conditions(allWindows, modFitlerType, allModFreq);


for winVar= 1:length(allWindows)
    TimeResolution=allWindows(winVar);
    
    figure(winVar);
    clf;
    
    for filtVar=1:length(modFitlerType)
        curModFiltType= modFitlerType{filtVar};
        
%         subplot(length(SNRs2use), length(modFitlerType), filtVar);
        cur_win_filttype_data= nan(2, length(SNRs2use), length(allModFreq)); % 2
        
        for modVar=1:length(allModFreq)
            modFreq=allModFreq(modVar);
            
            fName2Load=sprintf('%suRate_uR_f0_%s_%d_win_%.0fms.mat', loading_Dir, curModFiltType, modFreq, TimeResolution*1e3);
            
%             figure(cache_figHandle);
%             clf;
%             ax= nan(length(SNRs2use), 1);
            
            temp= load(fName2Load);
            mr_corr_Data= temp.mr_corr_Data;
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
%                 figure(cache_figHandle);
                ax(snrVar)=subplot(length(SNRs2use), 2, snrVar);
                hold on;
                
                %
                %                 plot3(1-[mr_corr_Data(validINDs.NH).SSNuncorr_sn_n_pos_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_pos_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_n_pos_Final], 'ob');
                %                 plot3(1-[mr_corr_Data(validINDs.NH).SSNuncorr_sn_n_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_n_neg_Final], 'ob');
                %                 plot3(1-[mr_corr_Data(validINDs.HI).SSNuncorr_sn_n_pos_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_sn_pos_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_n_pos_Final], 'dr');
                %                 plot3(1-[mr_corr_Data(validINDs.HI).SSNuncorr_sn_n_neg_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_sn_neg_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_n_neg_Final], 'dr');
                %                 plot3([0 1], [0 1], [0 1], 'k')
                
                subplot(length(SNRs2use), 2, 1);
                hold on;
                plot([mr_corr_Data(validINDs.NH).SSNcorr_s_n_pos_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_pos_Final], 'ob');
                plot([mr_corr_Data(validINDs.NH).SSNcorr_s_n_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_neg_Final], 'ob');
                plot([mr_corr_Data(validINDs.HI).SSNcorr_s_n_pos_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_sn_pos_Final], 'dr');
                plot([mr_corr_Data(validINDs.HI).SSNcorr_s_n_neg_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_sn_neg_Final], 'dr');
                [~,pVal]=ttest( ...
                    [ [mr_corr_Data(validINDs.NH).SSNcorr_s_n_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_n_pos_Final]], ...
                    [ [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_n_pos_Final] ]);
                
                tx= text(0.1, .6, sprintf('$%s |  F_{co}=%.0f Hz | tRes=%.0fms | pVal=%.4f$', curModFiltType, modFreq, TimeResolution*1e3, pVal), 'interpreter', 'latex');
                
                grid on;
                xlabel('corr (S, N)');
                ylabel('corr (S, SN)');
                set(gca, 'fontsize',  fSize);
                xlim([0 .7]);
                ylim([0 .7]);
                title(sprintf('SNR = %d dB', curSNR));
                plot([0 1], [0 1], 'k', 'linew', 2);
                
                subplot(length(SNRs2use), 2, 2);
                hold on;
                plot(1-[mr_corr_Data(validINDs.NH).SSNuncorr_sn_n_pos_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_pos_Final], 'ob');
                plot(1-[mr_corr_Data(validINDs.NH).SSNuncorr_sn_n_neg_Final], [mr_corr_Data(validINDs.NH).SSNcorr_s_sn_neg_Final], 'ob');
                plot(1-[mr_corr_Data(validINDs.HI).SSNuncorr_sn_n_pos_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_sn_pos_Final], 'dr');
                plot(1-[mr_corr_Data(validINDs.HI).SSNuncorr_sn_n_neg_Final], [mr_corr_Data(validINDs.HI).SSNcorr_s_sn_neg_Final], 'dr');
                
                [n_sn_corr_mu, n_sn_corr_dPrime]= calc_dPrime([ (1-[mr_corr_Data(validINDs.NH).SSNuncorr_sn_n_pos_Final]) (1-[mr_corr_Data(validINDs.NH).SSNuncorr_sn_n_neg_Final])], ...
                    [(1-[mr_corr_Data(validINDs.HI).SSNuncorr_sn_n_pos_Final]) (1-[mr_corr_Data(validINDs.HI).SSNuncorr_sn_n_neg_Final])]);
                
                grid on;
                xlabel('corr (SN, N)');
                ylabel('corr (SN, S)');
                set(gca, 'fontsize',  fSize);
                xlim([0 .7]);
                ylim([0 .7]);
                plot([0 1], [0 1], 'k', 'linew', 2);
                
                tx= text(0.1, .6, sprintf('$%s |  F_{co}=%.0f Hz | tRes=%.0fms$', curModFiltType, modFreq, TimeResolution*1e3), 'interpreter', 'latex');
                
                title(sprintf('SNR = %d dB', curSNR));
                grid on;
            end
            
            %%
%             figure(cache_figHandle);
            %             linkaxes(ax, 'xy');
            
            if saveFigs
%                 figure(cache_figHandle);
                set(gcf, 'units', 'inches', 'position', [1 1 14 4]);
                figName2Save=sprintf('%suRate_uR_f0_%s_%d_win_%.0fms_2d', outfigDir, curModFiltType, modFreq, TimeResolution*1e3);
                saveas(gcf, figName2Save, 'tiff');
                saveas(gcf, figName2Save);
            end
%             view([90 0])
        end
    end
end