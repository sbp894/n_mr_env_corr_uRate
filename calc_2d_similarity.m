function corr_2d=calc_2d_similarity(x,y)

corr_2d=sqrt(abs(sum(sum(x.*y)))/norm(x, 'fro')/norm(y, 'fro'));