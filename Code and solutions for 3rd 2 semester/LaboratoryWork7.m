clear;
clc;

t0 = 1;
t1 = 2;

alpha0 = 1; 
betta0 = 0;
gamma0 = -2*log(2);

alpha1 = 1; 
betta1 = 0; 
gamma1 = 0;

a = 0.5; b = 1;
tau = 0.01;
n = (b-a)/tau + 1;
x = a:tau:b;
U = log(x)./x;

p = @(x) 1/x;
q = @(x) 1 / x^2;
f = @(x) -2 / x.^3;
fSt = @(x) 0;

t1a = (gamma0-betta0*t0)/alpha0;
s1a = t0;

t2a = (gamma0-betta0*t1)/alpha0;
s2a = t1;

[T1, S1] = eulerSolve(p, q, f, x, tau, t1a, s1a, n);
[T2, S2] = eulerSolve(p, q, f, x, tau, t2a, t2a, n);

Z1 = T1; 
Z1der = S1;
Z2 = T2; 
Z2der = S2;

C = (gamma1-alpha1.*Z1(end)-betta1.*Z1der(end))/(alpha1.*(Z2(end)-Z1(end))+betta1.*(Z2der(end)-Z1der(end)));

u = (1-C).*Z1 + C.*Z2 ;

figure 
plot(x,U,x,u,'o');
legend("Exact solution","pristrelka");

function [T, S] = eulerSolve(p, q, f, x, tau, ta, sa, n)
  T(1) = ta;
  S(1) = sa;
  for k = 1:n - 1
    T(k + 1) = T(k) + tau * S(k);
    S(k + 1) = S(k) + tau * (-p(x(k)) * S(k) + q(x(k)) * T(k) + f(x(k)) );
  end
end