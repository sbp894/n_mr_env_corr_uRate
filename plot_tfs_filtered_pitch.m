function plot_tfs_filtered_pitch(vecINOrg_pos, vecINOrg_neg, fsOrg)

fs=2e3;
vecIN_pos=resample(vecINOrg_pos, 1, fsOrg/fs);
vecIN_neg=resample(vecINOrg_neg, 1, fsOrg/fs);
vecIN_tfs=(vecIN_pos-vecIN_neg)/2;
vecIN_env=(vecIN_pos+vecIN_neg)/2;

time=(1:length(vecIN_pos))/fs;
f0Center=128;


%% first filter is low pass
% N_lp=8;
% Fstop_lp=sqrt(AllmodFreqs(1)*AllmodFreqs(2));
% lpFilt = designfilt('lowpassiir','FilterOrder',N_lp, ...
%     'PassbandFrequency',Fstop_lp,'PassbandRipple',0.2, ...
%     'SampleRate',fs);
% % fvtool(lpFilt)
% filteredSignals{1}=filtfilt(lpFilt, vecIN);


%% All other filters are band-pass

N_bp=20;
HalfPowerFrequency1=f0Center/sqrt(2);
HalfPowerFrequency2=f0Center*sqrt(2);

bpFilt = designfilt('bandpassiir','FilterOrder',N_bp, ...get(p,props)
    'HalfPowerFrequency1',HalfPowerFrequency1,'HalfPowerFrequency2',HalfPowerFrequency2, ...
    'SampleRate',fs);
hold on;
%     fvtool(bpFilt)
f0_filtered_tfs=filtfilt(bpFilt, vecIN_tfs);
f0_filtered_env=filtfilt(bpFilt, vecIN_env);

nMulShift=1;
hold on;
plot(time, vecIN_pos);
plot(time, vecIN_neg);
plot(time,f0_filtered_tfs-nMulShift, '-', 'linew', 2);
plot(time,f0_filtered_env-2*nMulShift, '-', 'linew', 2);
xlabel('time');
title('128 Hz filters (octave)');
grid on;

legend('ur+', 'ur-', 'tfs-pitch', 'env-pitch');