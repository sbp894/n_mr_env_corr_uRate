function dPrime=calc_dPrime(x, y)

dPrime=(mean(x)-mean(y))/sqrt(.5*(var(x)+var(y)));