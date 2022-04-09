r0 = 0.00028;
df0_5 = 0.9979;
df1 = 0.9874;
df1_5 = 0.98215;
df2 = 0.9769;
df2_5 = 0.96725;
df3 = 0.9576;
df3_5 = 0.94835;
df4 = 0.9391;
df4_5 = 0.9289;
df5 = 0.9187;
%df0_5 = 0.9979;
%df1 = 0.9938;
%df1_5 = 0.9857;
%df2 = 0.9760;
%df2_5 = 0.9655;
%df3 = 0.9544;
%df3_5 = 0.9436;
%df4 = 0.9338;
%df4_5 = 0.9241;
%df5 = 0.9136;
df = [df0_5,df1,df1_5,df2,df2_5,df3,df3_5,df4,df4_5,df5];
kappa = 0.15;
sigma = 0.031;
t0 = 0;
T0 = 1;
TN = 4;
dt = 0.5;
dT = 1;
c1 = 0.05;
c2 = 0.03;
c3 = 0.01;
rhigh = 0.1;
rlow = 0.03;
K = 99.2;
F = 100;
%-------------------------------------Monte Carlo for different coupon rate

theta = [0.33255,0.14491,0.15236,0.28243,0.29453,0.29625,0.30958,0.35053,0.36829];
rand_seed1 = 9;
rand_case1 = 10000;
BMP1 = randn(rand_case1,rand_seed1);
R_MC = [];
R_MC(1:rand_case1,1) = r0;
for i = 1:rand_case1
    for j = 1:rand_seed1
        R_MC(i,j+1) = kappa*(theta(j)-R_MC(i,j))*0.5 + sigma*sqrt(0.5)*BMP1(i,j);
    end
end
C = [];
for j = 3:2:9
   for i = 1: rand_case1
       if 1/dT*(1/exp(-R_MC(i,j)*dt-R_MC(i,j+1)*dt)-1)>=rhigh
           C(i,j) = F*c1;
       elseif 1/dT*(1/exp(-R_MC(i,j)*dt-R_MC(i,j+1)*dt)-1)<=rlow
           C(i,j) = F*c3;
       else
           C(i,j) = F*c2;
       end
   end
end


V = zeros(50,50);
%-------------------- EU case
V(1:rand_case1,9) = F;
for i = 1:rand_case1
   for j = 8:-1:3
      if j==8||j==6||j==4
          V(i,j) = exp(-dt*R_MC(i,j))*(V(i,j+1)+C(i,j+1));
      else
          V(i,j) = exp(-dt*R_MC(i,j))*V(i,j+1);
      end
   end
end
for i = 1:rand_case1
  V(i,2) = exp(-dt*R_MC(i,2))*max([V(i,3)+C(i,3)-K,0]);
end
V(1:rand_case1,1) = V(1:rand_case1,2).*exp(-dt*r0);
V0 = sum(V(1:rand_case1,1))/rand_case1;
