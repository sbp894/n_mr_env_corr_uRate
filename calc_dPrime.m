function [dPrime, delta_mu]=calc_dPrime(x, y)

dPrime=(mean(x)-mean(y))/sqrt(.5*(var(x)+var(y)));
delta_mu= mean(x)-mean(y);