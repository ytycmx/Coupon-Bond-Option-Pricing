x = [0,0.5,2,3,4,5,7,10,20,30];
y = [0.028/100,0.412/100,1.167/100,1.443/100,1.568/100,1.692/100,1.751/100,1.852/100,1.993/100,1.948/100];
xx = [];
for i = 1:61
   xx(i) = 0+0.5*(i-1); 
end
yy = spline(x,y,xx);
plot(x,y,'o',xx,yy);

r_swap = yy;
DF = zeros(1,61);
DF(1) = 1;
DF(2) = 1/(1+r_swap(2)*0.5);
for i = 3:61
    DF(i) = 1/(r_swap(i-1)*0.5 + 1);
   for j =1:i-1
       DF(i) = DF(i) - r_swap(i-1)*0.5/(r_swap(i-1)*0.5 + 1)*(DF(j));
   end
end
plot(xx,DF);
r_spot = zeros(1,61);
r_spot(1) = yy(1);
r_spot(2) = yy(2);
for i = 3:61
   r_spot(i) = -log(DF(i))/xx(i); 
end
plot(xx,r_spot);
