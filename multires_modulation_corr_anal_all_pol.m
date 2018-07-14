function [corrValsModFreq, corrValsModFreqNF]=multires_modulation_corr_anal(S_rate_envOrg, SN_rate_envOrg, N_rate_envOrg, fsOrg, AllmodFreqs)

plotModFiltSignals=0;
fs=2e3;

corrValsModFreq=nan(length(AllmodFreqs),1);
corrValsModFreqNF=nan(length(AllmodFreqs),1);

S_rate_env=resample(S_rate_envOrg, 1, fsOrg/fs);
SN_rate_env=resample(SN_rate_envOrg, 1, fsOrg/fs);
N_rate_env=resample(N_rate_envOrg, 1, fsOrg/fs);

time=(1:length(S_rate_env))/fs;
stimDur=length(S_rate_env)/fs;

filtered_S_rate_env=cell(length(AllmodFreqs),1);
filtered_SN_rate_env=cell(length(AllmodFreqs),1);
filtered_N_rate_env=cell(length(AllmodFreqs),1);

%% first filter is low pass
N_lp=8;
Fstop_lp=sqrt(AllmodFreqs(1)*AllmodFreqs(2));
lpFilt = designfilt('lowpassiir','FilterOrder',N_lp, ...
    'PassbandFrequency',Fstop_lp,'PassbandRipple',0.2, ...
    'SampleRate',fs);
% fvtool(lpFilt)
filtered_S_rate_env{1}=filtfilt(lpFilt, S_rate_env);
filtered_SN_rate_env{1}=filtfilt(lpFilt, SN_rate_env);
filtered_N_rate_env{1}=filtfilt(lpFilt, N_rate_env);


%% All other filters are band-pass
for bpfVar=2:length(filtered_S_rate_env)
    
    N_bp=20;
    HalfPowerFrequency1=sqrt(AllmodFreqs(bpfVar)*AllmodFreqs(bpfVar-1));
    if bpfVar~=length(filtered_S_rate_env)
        HalfPowerFrequency2=sqrt(AllmodFreqs(bpfVar)*AllmodFreqs(bpfVar+1));
    else
        HalfPowerFrequency2=sqrt(2)*AllmodFreqs(bpfVar);
    end
    
    bpFilt = designfilt('bandpassiir','FilterOrder',N_bp, ...get(p,props)
        'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
        'SampleRate',fs);
    %     fvtool(bpFilt)
    filtered_S_rate_env{bpfVar}=filtfilt(bpFilt, S_rate_env);
    filtered_SN_rate_env{bpfVar}=filtfilt(bpFilt, SN_rate_env);
    filtered_N_rate_env{bpfVar}=filtfilt(bpFilt, N_rate_env);
end


for modFreqVar=1:length(AllmodFreqs)
    curModFreq=AllmodFreqs(modFreqVar);
    curTimeResolution=1/curModFreq;
    curModFiltData_S=filtered_S_rate_env{modFreqVar};
    curModFiltData_SN=filtered_SN_rate_env{modFreqVar};
    curModFiltData_N=filtered_N_rate_env{modFreqVar};
    if rem(stimDur, curTimeResolution) < .5*curTimeResolution
        curCorrVals=nan(floor(stimDur/curTimeResolution), 1);
        curCorrValsNF=nan(floor(stimDur/curTimeResolution), 1);
    else
        curCorrVals=nan(ceil(stimDur/curTimeResolution), 1);
        curCorrValsNF=nan(ceil(stimDur/curTimeResolution), 1);
    end
    for windowVar=1:length(curCorrVals)
        tStart=(windowVar-1)*curTimeResolution;
        indStart=max(1, round(tStart*fs));
        tEnd=windowVar*curTimeResolution;
        indEnd=min(length(curModFiltData_S), round(tEnd*fs));
        validINDs=indStart:indEnd;
        curCorrVals(windowVar)=corr2(curModFiltData_S(validINDs), curModFiltData_SN(validINDs));
        curCorrValsNF(windowVar)=corr2(curModFiltData_S(validINDs), curModFiltData_N(validINDs));
    end
    curCorrVals(curCorrVals<0)=0;
    curCorrValsNF(curCorrValsNF<0)=0;
    corrValsModFreq(modFreqVar)=sqrt(nansum(curCorrVals.^2));
    corrValsModFreqNF(modFreqVar)=sqrt(nansum(curCorrValsNF.^2));
end


if plotModFiltSignals
    lineNum=nan(length(AllmodFreqs), 1);
    nMulShift=3;
    for modFreqVar=1:length(AllmodFreqs)
        curNormFilteredSig=filtered_S_rate_env{modFreqVar};%/max(filteredSignals{modFreqVar});
        hold on;
        plot(time, S_rate_env);
        if AllmodFreqs(modFreqVar)>10 %from Relano-Iborra et. al. 2016
            lineNum(modFreqVar)=plot(time,curNormFilteredSig-nMulShift*modFreqVar, '-.', 'linew', .5);
            plot(time,abs(hilbert(curNormFilteredSig))-nMulShift*modFreqVar, 'color', get(lineNum(modFreqVar), 'color'), 'linew', 2)
        else
            lineNum(modFreqVar)=plot(time,curNormFilteredSig-nMulShift*modFreqVar, '-', 'linew', 2);
        end
    end
    xlabel('time');
    title('meanrate and mod filtered');
    grid on;
    legend(lineNum, num2str(AllmodFreqs'));
end