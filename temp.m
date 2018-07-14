clear;
clc;
figHan.sllRes=2;

SimulationOutputDir=[fileparts(pwd) filesep 'SIM_OUT' filesep];
if ~isdir(SimulationOutputDir)
    mkdir(SimulationOutputDir);
end
saveAllFisg=1;

all_tRes=[.1 .2 .5 1 2 5 10 20 50]*1e-3;
lw=2;

for tResVar=1:length(all_tRes)
    tRes=all_tRes(tResVar);
    allSNRs=[-10 -5 0];
    allCFs=logspace(log10(200), log10(8e3), 25);
    ssnCorrs=nan(length(allSNRs), length(allCFs));
    flnCorrs=nan(length(allSNRs), length(allCFs));
    ssnCorrs_NF=nan(length(allSNRs), length(allCFs));
    flnCorrs_NF=nan(length(allSNRs), length(allCFs));
    
    for snrVar=1:length(allSNRs)
        snr=allSNRs(snrVar);
        for cfVar=1:length(allCFs)
            CF=allCFs(cfVar);
            
            plotAll=0;
            switch snr
                case -10
                    stimDir='../Stimuli/SNR_min10/';
                case -5
                    stimDir='../Stimuli/SNR_min5/';
                case 0
                    stimDir='../Stimuli/SNR_0/';
            end
            
            [sOrg, ~]=audioread('../Stimuli/SNR_min5/SSN_Stim_S_N.wav');
            [fln_Org, ~]=audioread([stimDir '/FLN_Stim' num2str(snr) 'dB_N_P.wav']);
            [s_fln_Org, ~]=audioread([stimDir 'FLN_Stim' num2str(snr) 'dB_SN_P.wav']);
            [ssn_Org, ~]=audioread([stimDir 'SSN_Stim' num2str(snr) 'dB_N_P.wav']);
            [s_ssn_Org, fsOrg]=audioread([stimDir 'SSN_Stim' num2str(snr) 'dB_SN_P.wav']);
            
            stimDur=length(sOrg)/fsOrg;
            
            
            %%
            fs=20e3;
            QbyP=round(fsOrg/fs);
            
            
            s=resample(sOrg, 1, QbyP);
            s_ssn=resample(s_ssn_Org, 1, QbyP);
            ssn=resample(ssn_Org, 1, QbyP);
            s_fln=resample(s_fln_Org, 1, QbyP);
            fln=resample(fln_Org, 1, QbyP);
            
            tStim=(1:length(s))/fs;
            
            N=6;
            F3dB1=CF*.9;
            F3dB2=min(CF*1.05, fs/2);
            
            %%
            d=fdesign.bandpass('N,F3dB1,F3dB2', N,F3dB1,F3dB2, fs);
            dMethod=d.designmethods{1};
            hd=design(d, dMethod);
            % fvtool(hd);
            
            %%
            s_Resp=filter(hd, s);
            s_ssn_Resp=filter(hd, s_ssn);
            ssn_Resp=filter(hd, ssn);
            s_fln_Resp=filter(hd, s_fln);
            fln_Resp=filter(hd, fln);
            
            %%
            sRespEnv=abs(hilbert(s_Resp));
            s_ssn_RespEnv=abs(hilbert(s_ssn_Resp));
            ssn_RespEnv=abs(hilbert(ssn_Resp));
            s_fln_RespEnv=abs(hilbert(s_fln_Resp));
            fln_RespEnv=abs(hilbert(fln_Resp));
            
            %%
            if plotAll
                figure(1);
                clf;
                
                Ashift=max([s_Resp; s_ssn_Resp; ssn_Resp]);
                hold on;
                plot(tStim, 2*Ashift+s_Resp, ':');
                plot(tStim, Ashift+s_ssn_Resp, ':');
                plot(tStim, ssn_Resp, ':');
                plot(tStim, -Ashift+s_fln_Resp, ':');
                plot(tStim, -2*Ashift+fln_Resp, ':');
                
                plot(tStim, 2*Ashift+sRespEnv, 'linew', lw);
                plot(tStim, Ashift+s_ssn_RespEnv, 'linew', lw);
                plot(tStim, ssn_RespEnv, 'linew', lw);
                plot(tStim, -Ashift+s_fln_RespEnv, 'linew', lw);
                plot(tStim, -2*Ashift+fln_RespEnv, 'linew', lw);
                
                xlabel('time (sec)');
            end
            %%
            nSamplesInBin=round(tRes*fs);
            sRespEnv_meanrate=mean(reshape(sRespEnv, nSamplesInBin, length(sRespEnv)/nSamplesInBin),1);
            s_ssn_RespEnv_meanrate=mean(reshape(s_ssn_RespEnv, nSamplesInBin, length(s_ssn_RespEnv)/nSamplesInBin),1);
            ssn_RespEnv_meanrate=mean(reshape(ssn_RespEnv, nSamplesInBin, length(ssn_RespEnv)/nSamplesInBin),1);
            s_fln_RespEnv_meanrate=mean(reshape(s_fln_RespEnv, nSamplesInBin, length(s_fln_RespEnv)/nSamplesInBin),1);
            fln_RespEnv_meanrate=mean(reshape(fln_RespEnv, nSamplesInBin, length(fln_RespEnv)/nSamplesInBin),1);
            
            ssnCorrs(snrVar, cfVar)=corr2(sRespEnv_meanrate, s_ssn_RespEnv_meanrate);
            ssnCorrs_NF(snrVar, cfVar)=corr2(sRespEnv_meanrate, ssn_RespEnv_meanrate);
            flnCorrs(snrVar, cfVar)=corr2(sRespEnv_meanrate, s_fln_RespEnv_meanrate);
            flnCorrs_NF(snrVar, cfVar)=corr2(sRespEnv_meanrate, fln_RespEnv_meanrate);
        end
    end
    
    %%
    figure(figHan.sllRes);
    subplot(3,3,tResVar);
    deltaX=.1;
    plot(nan, nan, 'vb-', 'linew', lw);
    hold on;
    plot(nan, nan, '<c', 'linew', lw);
    plot(nan, nan, '^r-', 'linew', lw);
    plot(nan, nan, 'm>', 'linew', lw);
        
    plot(repmat(allSNRs', 1, size(ssnCorrs,2))-2*deltaX,ssnCorrs, 'bv');
    plot(repmat(allSNRs', 1, size(ssnCorrs_NF,2))-deltaX,ssnCorrs_NF, 'c<');
    plot(repmat(allSNRs', 1, size(flnCorrs,2))+deltaX,flnCorrs_NF, 'm>');
    plot(repmat(allSNRs', 1, size(flnCorrs,2))+2*deltaX,flnCorrs, 'r^');
    plot(allSNRs-.1,nanmean(ssnCorrs,2), 'b-', 'linew', lw);
    plot(allSNRs+.1,nanmean(flnCorrs,2), 'r-', 'linew', lw);
    xlim([-11 1]);
    hold off;
    title(sprintf('tRes=%.1f ms', tRes*1e3));
end
figure(2);
set(gcf, 'units', 'inches', 'position', [1 1 8 8]);
subplot(334); ylabel('S/SN-env corr');
subplot(338); xlabel('SNR')
subplot(331); legend('SSN', 'SSN-NF', 'FLN', 'FLN-NF', 'location', 'southeast'); 
set(gca, 'fontsize', 16);
if saveAllFisg
    saveas(figHan.sllRes, [SimulationOutputDir 'uR_corr_allSimula_allRes'], 'tiff');
end