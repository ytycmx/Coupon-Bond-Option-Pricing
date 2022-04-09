r0 = 0.00028;
t0 = 0;
T0 = 1;
dt = 1/360;
df0 = 1;
df0_5 = 0.9979;
df1 = 0.9874;
df_MC = [];
df_MC(1) = df0;
for i = 2: T0/dt/2 + 1
    df_MC(i) = df0 - (df0-df0_5)*dt*2*(i-1);
end
for i = T0/dt/2 + 2 : T0/dt +1
    df_MC(i) = df0_5 - (df0_5-df1)*dt*2*(i-(T0/dt/2 +1));
end
kappa = 0.15;
sigma = 0.01/100;
c = 0.02;
K = 99.2;
F = 100;
df2 = 0.9769;
df3 = 0.9576;
df4 = 0.9391;

%df2 = 0.9760;
%df3 = 0.9544;
%df4 = 0.9338;

df_initial = [df2,df3,df4];
V = [];
for q = 1:100
%---------------------------------------obtain the Browian Paths (3*500*360)
rand_seed = 360;
rand_case = 3*500;
BMP = randn(rand_case,rand_seed);
%---------------------------------------MC itself
DP = [];
for i = 1:3
    DP(500*(i-1)+1:i*500,1) = df_initial(i);
end
for i = 1:rand_case
   for j = 2:rand_seed+1
      DP(i,j) =  DP(i,j-1)*df_MC(j-1)/df_MC(j)*exp(-0.5*sigma*sigma*dt+sigma*sqrt(dt)*BMP(i,j-1));
   end
end
%---------------------------------------average the payoff
DP_average = [];
for i = 1:rand_case/500
   DP_average(i) = sum(DP(500*(i-1)+1:i*500, rand_seed+1))/500;
end
V(q) = max([c*F*DP_average(1)*df1+c*F*DP_average(2)*df1+c*F*DP_average(3)*df1+F*DP_average(3)*df1-K*df1, 0]);
end
V_fin = sum(V)/100;
%-----------------------------------------MC for AM option
DP_AM = [];
for j = 2:rand_seed+1
   for i = 1:rand_case/500
      DP_AM(i,j) = sum(DP(500*(i-1)+1:i*500, j))/500;
   end
end
V_AM_INI = [];
for j = 2:rand_seed+1
      V_AM_INI(j) = c*F*DP_AM(1,j)*df1+c*F*DP_AM(2,j)*df1+c*F*DP_AM(3,j)*df1+F*DP_AM(3,j)*df1-K*df1;
end
V_AM = [];
for j = 2:rand_seed+1
    V_AM(j) = max([V_AM_INI(j), 0]);
end
V_AM_FINAL = max(V_AM);
