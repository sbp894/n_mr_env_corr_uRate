function plot_filtered_modulation(vecINOrg, fsOrg)

fs=2e3;
vecIN=resample(vecINOrg, 1, fsOrg/fs);

time=(1:length(vecIN))/fs;
AllmodFreqs=2.^(0:8);

filteredSignals=cell(length(AllmodFreqs),1);

%% first filter is low pass
N_lp=8;
Fstop_lp=sqrt(AllmodFreqs(1)*AllmodFreqs(2));
lpFilt = designfilt('lowpassiir','FilterOrder',N_lp, ...
    'PassbandFrequency',Fstop_lp,'PassbandRipple',0.2, ...
    'SampleRate',fs);
% fvtool(lpFilt)
filteredSignals{1}=filtfilt(lpFilt, vecIN);


%% All other filters are band-pass
for bpfVar=2:length(filteredSignals)
    
    N_bp=20;
    HalfPowerFrequency1=sqrt(AllmodFreqs(bpfVar)*AllmodFreqs(bpfVar-1));
    if bpfVar~=length(filteredSignals)
        HalfPowerFrequency2=sqrt(AllmodFreqs(bpfVar)*AllmodFreqs(bpfVar+1));
    else
        HalfPowerFrequency2=sqrt(2)*AllmodFreqs(bpfVar);
    end
    
    bpFilt = designfilt('bandpassiir','FilterOrder',N_bp, ...get(p,props)
        'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
        'SampleRate',fs);
    hold on;
    %     fvtool(bpFilt)
    filteredSignals{bpfVar}=filtfilt(bpFilt, vecIN);
end

lineNum=nan(length(AllmodFreqs), 1);
nMulShift=2;
for modFreqVar=1:length(AllmodFreqs)
    hold on;
    plot(time, vecIN);
    curNormFilteredSig=filteredSignals{modFreqVar};%/max(filteredSignals{modFreqVar});
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
legend(lineNum, num2str(AllmodFreqs'), 'location', 'eastoutside');