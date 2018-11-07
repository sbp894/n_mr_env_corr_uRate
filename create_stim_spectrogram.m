clear;
clc;
clf;

[x,fs]= audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_08_15-Q361_AN_PTS/Signals/MH/SNRenv/SNR_0/FLN_Stim_S_P.wav');
time= (1:length(x))/fs;

ax(1)=subplot(5,1,5);
plot(time ,x);
xlabel('time (sec)');
ylabel('Stim Amp');

ax(2)=subplot(5,1,1:4);
nfft=1028;
tRes= 40e-3;
tOvLap= 20e-3;
nWind= round(tRes*fs);
novlap= round(tOvLap*fs);
freq2show= logspace(log10(100), log10(10e3), nfft);
spectrogram(x, nWind, novlap, freq2show, fs, 'yaxis');
% set(gca, 'yscale', 'log');
grid on;
xlabel('');
colorbar off;
linkaxes(ax, 'x');

xlim([0 length(x)/fs]);
ylim([0 8]);
set(gcf, 'units', 'inches', 'position', [1 1 10 8]);
saveas(gcf, '/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/outFig/danish_stim_spectrogram', 'tiff');
saveas(gcf, '/media/parida/DATAPART1/Matlab/SNRenv/n_mr_env_corr_uRate/outFig/danish_stim_spectrogram', 'png');