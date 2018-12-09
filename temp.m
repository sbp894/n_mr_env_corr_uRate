%% Coherence Estimate of Two Sequences
% Compute and plot the coherence estimate between two colored noise
% sequences.

%%
% Generate a signal consisting of white Gaussian noise.

clf;

r = randn(16384,1);

%%
% To create the first sequence, bandpass filter the signal. Design a
% 16th-order filter that passes normalized frequencies between 0.2&pi; and
% 0.4&pi; rad/sample. Specify a stopband attenuation of 60 dB. Filter the
% original signal.

dx = designfilt('bandpassiir','FilterOrder',16, ...
    'StopbandFrequency1',0.2,'StopbandFrequency2',0.4, ...
    'StopbandAttenuation',60);
x = filter(dx,r);

%%
% To create the second sequence, design a 16th-order filter that stops
% normalized frequencies between 0.6&pi; and 0.8&pi; rad/sample. Specify a
% passband ripple of 0.1 dB. Filter the original signal.

dy = designfilt('bandstopiir','FilterOrder',16, ...
    'PassbandFrequency1',0.6,'PassbandFrequency2',0.8, ...
    'PassbandRipple',0.1);
y = filter(dy,r);

%%
% Estimate the magnitude-squared coherence of |x| and |y|. Use a 512-sample
% Hamming window. Specify 500 samples of overlap between adjoining segments
% and 2048 DFT points.


plot_mscohere(x, y, 2, 'hamming', 'FracWindow', 512/length(x), 'xscale', 'lin', 'FracOverlap', 500/512, 'nfft', 2048);
[cxy,fc] = mscohere(x,y,hamming(512),500,2048);
hold on;
plot(fc/pi, db(cxy)/2, '--', 'linew', 2)
