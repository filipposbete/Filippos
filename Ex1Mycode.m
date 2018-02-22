mean=0
sigma=1
M=mean+sigma*randn(1000, 2)
R=[1, 0.10; 0.10, 1]
L=chol(R)
M=M*L
x=M(:, 1)
y=M(:, 2)

corr(x,y)

myreg = regstats(y,x, 'linear',{'tstat'});
