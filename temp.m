% clear;
clc;
close all;

x = S_rate_plus;
y = SN_rate_plus;

plot(S_rate_plus)
hold on
plot(SN_rate_plus)


t_start= .32;
t_end= .58;
fs_data= 1/ rate_binWidth;
ind_start= round(fs_data*t_start);
ind_end= round(fs_data*t_end);

s_snap = S_rate_plus(ind_start:ind_end);
sn_snap = SN_rate_plus(ind_start:ind_end);
n_snap = N_rate_plus(ind_start:ind_end);

t_mscohere= .1;
win_mscohere= hamming(round(fs_data*t_mscohere));
win_ovlap= round(fs_data*t_mscohere/2);
nfft= 2^(1+nextpow2(length(win_mscohere)));

s_snap= s_snap-mean(s_snap);
sn_snap= sn_snap-mean(sn_snap);
n_snap= n_snap-mean(n_snap);

[cxy,fc] = mscohere(s_snap, sn_snap, win_mscohere, win_ovlap, nfft, fs_data);

figure
plot(fc,cxy)

