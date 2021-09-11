clear;
clc;

alpha0 = 1; 
betta0 = 0;
gamma0 = -2*log(2);
alpha1 = 1; 
betta1 = 0; 
gamma1 = 0;

a = 0.5; b = 1;

tau = 0.005;
n = (b-a)/tau + 1;
h =  0.005;
x = a:tau:b;

alpha01 = alpha0-betta0/h; 
betta01 = betta0/h; 
gamma01 = gamma0;
alpha02 = -betta1/h; 
betta02 = alpha1+betta1/h; 
gamma02 = gamma1;

X1=-betta01/alpha01;
Z1=gamma01/alpha01;

u = log(x)./x;

a = @(x) 1+h/2*1/x;
b = @(x) -(2+h^2*1/x^2);
c = @(x) 1-h/2*1/x;
t = @(x) h^2*(-2 / x.^3);

[X,Z] = solution(a, b, c, t, x, X1, Z1, n);

yend=(gamma02-alpha02*Z(end))/(betta02+alpha02*X(end));

[Y] = backSolution(yend, X, Z, n);

figure 
plot(x,Y,x, u,'o');
legend("Progonka", "Exact solution");


function [X,Z] = solution(a, b, c, t, x, X1, Z1, n)
  X(1)=X1;
  Z(1)=Z1;
  for k = 1:n - 1
   X(k+1) = -a(x(k))/(b(x(k))+c(x(k))*X(k));
  end
  for k = 1:n - 1
    Z(k+1) = (t(x(k))-c(x(k))*Z(k))/(b(x(k))+c(x(k))*X(k));
  end
end

function [Y] = backSolution(yend, X, Z, n)
  Y(n) = yend;
  for k = n:-1:2
   Y(k-1) = Y(k)*X(k)+Z(k);
  end
end
