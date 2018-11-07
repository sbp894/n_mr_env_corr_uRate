clear;
clf;
clc;

[stim,fs]= audioread('/media/parida/DATAPART1/Matlab/ExpData/MatData/SP-2018_08_15-Q361_AN_PTS/Signals/MH/SNRenv/SNR_0/FLN_Stim_S_P.wav');
stim_time= (1:length(stim))/fs;

