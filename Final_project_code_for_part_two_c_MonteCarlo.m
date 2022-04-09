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
sigma = 0.031;
c = 0.02;
K = 99.2;
%K = 98;
%K = 99.85;
F = 100;
df2 = 0.9769;
df3 = 0.9576;
df4 = 0.9391;
df5 = 0.9187;

%df2 = 0.9760;
%df3 = 0.9544;
%df4 = 0.9338;
%df5 = 0.9136;

df_initial = [df2,df3,df4,df5];
V = [];
for q = 1:100
%---------------------------------------obtain the Browian Paths (3*1000*360 for 1yr-4yr, 4*1000*360 for 1yr-5yr)
rand_seed = 360;
rand_case = 3*1000;
%rand_case = 4*1000;
BMP = randn(rand_case,rand_seed);
%---------------------------------------MC itself
DP = [];
for i = 1:3
    DP(1000*(i-1)+1:i*1000,1) = df_initial(i);
end
for i = 1:rand_case
   for j = 2:rand_seed+1
      DP(i,j) =  DP(i,j-1)*df_MC(j-1)/df_MC(j)*exp(-0.5*sigma*sigma*dt+sigma*sqrt(dt)*BMP(i,j-1));
   end
end
%---------------------------------------average the payoff
DP_average = [];
for i = 1:rand_case/1000
   DP_average(i) = sum(DP(1000*(i-1)+1:i*1000, rand_seed+1))/1000;
end
V(q) = max([c*F*df1+c*F*DP_average(1)*df1+c*F*DP_average(2)*df1+c*F*DP_average(3)*df1+F*DP_average(3)*df1-K*df1, 0]);
%V(q) = max([c*F*df1+c*F*DP_average(1)*df1+c*F*DP_average(2)*df1+c*F*DP_average(3)*df1+c*F*DP_average(4)*df1+F*DP_average(4)*df1-K*df1, 0]);
end
V_fin = sum(V)/100;








