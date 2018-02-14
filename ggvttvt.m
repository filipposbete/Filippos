mu=0;
sigma=1;
M=mu+sigma*randn(1000, 2);
R=[1, 0.1; 0.1, 1];
L=chol(R);
M=M*L;
x=M(:, 1);
y=M(:, 2);
temp = corr(x, y);

while (temp <0.0995 || temp >0.1005 )
   M=mu+sigma*randn(1000, 2);
    R=[1, 0.1; 0.1, 1];
    L=chol(R);
    M=M*L;
    x=M(:, 1);
    y=M(:, 2);
    temp = corr(x, y);    
end

corr(x, y)

%Regression 
reg1 = regstats(y,x, 'linear');
