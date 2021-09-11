clear;
clc;

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

ta=0;sa=0;
[T_1, S_1] = eulerSolve(p, q, f, x, tau, ta, sa, n);
sa=1;
[T_2, S_2] = eulerSolve(p, q, fSt, x, tau, ta, sa, n);
ta=1;sa=0;
[T_3, S_3] = eulerSolve(p, q, fSt, x, tau, ta, sa, n);

Z = T_1; 
Zdeg = S_1;
Z1 = T_2; 
Z1deg = S_2;
Z2 = T_3; 
Z2deg = S_3;

c10 = alpha0 * Z1(1) + betta0 * Z1deg(1);
c11 = alpha0 * Z2(1) + betta0 * Z2deg(1);
sub1 = gamma0 - alpha0 * Z(1) - betta0 * Zdeg(1);
c21 = alpha1 * Z1(end) + betta1 * Z1deg(end);
c22 = alpha1 * Z2(end) + betta1 * Z1deg(end);
sub2 = gamma1 - alpha1 * Z(end) - betta1 * Zdeg(end);

C = [c10 c11; c21 c22] \ [sub1; sub2];

u = Z + C(1)*Z1 + C(2)*Z2;

figure 
plot(x,U,x, u,'-o');
legend("exact solution", "variable method");

function [T, S] = eulerSolve(p, q, f, x, tau, ta, sa, n)
  T(1) = ta;
  S(1) = sa;
  for k = 1:n - 1
    T(k + 1) = T(k) + tau * S(k);
    S(k + 1) = S(k) + tau * (-p(x(k)) * S(k) + q(x(k)) * T(k) + f(x(k)) );
  end
end