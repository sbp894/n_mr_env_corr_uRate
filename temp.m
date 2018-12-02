clear;
clc;
clf;

[sig, fs]= audioread('/media/parida/DATAPART1/Matlab/SNRenv/SFR_sEPSM/shorter_stim/FLN_Stim_S_P.wav');

% sig =load
% fs = sampling frequency

%
timeResolution= 32e-3;
thresh_dB_minimum2plot=-80;


nRes=round(timeResolution*fs);
NFFT=2^(1 + nextpow2(nRes));



% no need to change
nRes=round(timeResolution*fs);
nOVLap= round(.33*nRes);
[~, F, T, P]= spectrogram(sig, blackman(nRes), nOVLap, NFFT, fs, 'yaxis', 'reassigned');

P(P<(10^(thresh_dB_minimum2plot/10)))=10^(thresh_dB_minimum2plot/10);
imagesc( T, F, P);
colorbar
set(gca,'YDir', 'normal');
ylim([0 200])