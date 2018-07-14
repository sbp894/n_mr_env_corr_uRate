function ExpControlParams=get_ExpControlParams()

ExpControlParams.sentences=3; %1:10; % just 1 constant sentence for data collection
ExpControlParams.BOOTSTRAP_fractiondata=.67; % How many spike lines to go into single estimate of SUMCOR
ExpControlParams.BootstrapLoopMax=24; % How many SUMCORs should be averaged for single estimate of SUMCOR (hence its PSD)
ExpControlParams.nPSDs2Avg=12; % How many SUMCOR-PSDs should be averaged for single estimate of modulation power (BootstrapLoopMax choose nPSDs2Avg)
ExpControlParams.BootstrapLoopReport=60; % How many bootstrap estimates of modulation power (and SNRenv)
ExpControlParams.fixSPL=0; % S=65dB SPL, Should SN be 65 dB SPL or S should be kept constant and different noises be added to create desired SNR. 
ExpControlParams.winCorr0Add1Mul=1; % 1=unbiased
ExpControlParams.fs=100e3; % Sampling Frequency
ExpControlParams.ModFreqs=[1,2,4,8,16,32,64, 128, 256, 512, 1024]; %Which modulation frequencies
ExpControlParams.SCC_onsetIGNORE_sec=50e-3; % Ignore the spikes occuring before this time limit
ExpControlParams.UseSlepian=1; % For averaging (BootstrapLoopMax choose nPSDs2Avg) PSDs to lower variance of the modulation power estimates
ExpControlParams.MAXdelay_fractionStimDur=.75; % Since multiplicativeWindow is used for correction of correlograms, don't go after MAXdelay_fractionStimDur*stimulus delay. 
ExpControlParams.DELAYbinwidth_sec=50e-6; % (Joris, 2003) to constuct correlograms
% ExpControlParams.NoiseOrder={'FLN','SSN'}';
ExpControlParams.ModSumExponent=2; %power to use while estimating SNRenv from modulation powers. DTU group uses 2. 
ExpControlParams.maxSpikes=inf;
ExpControlParams.THRESHOLD_percent_less_than_refractory=5; % threshold percent of ISI that is lower than this before censoring. The code to censor it 'screenDataMAT'.

% ExpControlParams.ohcLoss_dB=0; % not used in data
% ExpControlParams.ihcLoss_dB=0; % not used in data
% ExpControlParams.nRep=25; % not used in data
