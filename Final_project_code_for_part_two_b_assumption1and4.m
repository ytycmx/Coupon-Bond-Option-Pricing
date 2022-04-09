r0 = 0.00028;
df1 = 0.9874;
df2 = 0.9769;
df3 = 0.9576;
df4 = 0.9391;
df5 = 0.9187;
df = [df1,df2,df3,df4,df5];
kappa = 0.15;
sigma = 0.031;
t0 = 0;
T0 = 1;
TN = 4;
dt = 1;
dT = 1;
c1 = 0.05;
c2 = 0.03;
c3 = 0.01;
rhigh = 0.1;
rlow = 0.03;
K = 99.2;
F = 100;
J = floor(1/2/kappa/dt)+1;
deltar = sigma * sqrt(3*dt);
k = [];
k(1) = -J+1;
for j = 2:2*J
   k(j) = j - (J+1);
end
k(2*J+1) = J-1;
eta = [];
p = zeros(50,50);
para1 = 1-kappa*dt;
for i = 1:J
    for j = 1:(i-1)*2+1
       eta(i,j) = (para1*(i-j) - k(i+j+J-(1+2*(j-1)))) * deltar;
    end
end
for i = 1:J
   for j = 1:3:((i-1)*2+1)*3
       p(i,j) = 1/2*((sigma^2*dt + eta(i,floor((j-1)/3)+1)^2)/deltar^2 + eta(i,floor((j-1)/3)+1)/deltar);
   end
   for j = 2:3:((i-1)*2+1)*3
       p(i,j) =  1 - (sigma^2*dt + eta(i,floor((j-2)/3)+1)^2)/deltar^2;
   end
   for j = 3:3:((i-1)*2+1)*3
       p(i,j) = 1/2*((sigma^2*dt + eta(i,floor((j-3)/3)+1)^2)/deltar^2 - eta(i,floor((j-3)/3)+1)/deltar);
   end
end
r = [];
r(1,1) = r0;
for i = 2:5
   for j = 1:2*i-1
      r(i,j) = para1*r(i-1,i-1) + (i-j)*deltar;
   end
end
Q = zeros(50,50);
Q(1,1) = 1;
theta = [];
for i = 2
    para2 = 0;
    for j = 1:3
        Q(i,j) = Q(1,1)*p(i-1,j)*exp(-r(1,1)*dt);
        para2 = para2 + Q(i,j)*exp(-r(i,j)*dt);
    end
    theta(i-1) = (1/(kappa*dt*dt)) * log(para2/df(i));
    for j = 1:3
       r(i,j) = r(i,j) + kappa * theta(i-1) * dt;
    end
end
for i = 3:5
   para2 = 0;
   for j = 1
      Q(i,j) = Q(i-1,j)*p(i-1,j)*exp(-r(i-1,1)*dt);
   end
   for j = 2
       for n = 1:j
          Q(i,j) = Q(i,j) + Q(i-1,n)* p(i-1,j+(n-1)*2)*exp(-r(i-1,n)*dt);
       end
   end
   for j = 3:2*i-1
       Q(i,j) = Q(i-1,j-2)*p(i-1,3*(j-2))*exp(-r(i-1,j-2)*dt) + Q(i-1,j-1)*p(i-1,3*(j-2)+2)*exp(-r(i-1,j-1)*dt) + Q(i-1,j)*p(i-1,3*(j-2)+4)*exp(-r(i-1,j)*dt);
   end
   for j = 1:2*i-1
      para2 = para2 + Q(i,j)*exp(-r(i,j)*dt); 
   end
   theta(i-1) = (1/(kappa*dt*dt)) * log(para2/df(i));
   for j = 1:2*i-1
      r(i,j) = r(i,j) + kappa * theta(i-1) * dt;
   end
end

C = [];
for i = 2:5
   for j = 1:2*i-1
       if 1/dT*(1/exp(-r(i,j)*dt)-1)>=rhigh
           C(i,j) = F*c1;
       elseif 1/dT*(1/exp(-r(i,j)*dt)-1)<=rlow
           C(i,j) = F*c3;
       else
           C(i,j) = F*c2;
       end
   end
end

V = zeros(50,50);
for i = 5
   for j = 1:9
      V(i,j) = F; 
   end
end
for i = 4:-1:2
    for j = 1:2*i-1
       V(i,j) = exp(-dt*r(i,j))*(p(i,3*j-2)*(V(i+1,j)+C(i+1,j))+p(i,3*j-1)*(V(i+1,j+1)+C(i+1,j+1))+p(i,3*j)*(V(i+1,j+2)+C(i+1,j+2))); 
    end
end
for i = 5
   for j = 1:9
      V(i,j) = F + C(i,j); 
   end
end
V(1,1) = exp(-dt*r(1,1))*(p(1,1)*max([V(2,1)-K,0]) + p(1,2)*max([V(2,2)-K,0]) + p(1,3)*max([V(2,3)-K,0]));
