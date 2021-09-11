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
u = log(x)./x;

p = @(x) 1/x;
q = @(x) 1 / x^2;
f = @(x) -2 / x.^3;
U1a = 0;
U2a = -2*log(2);

[U1, U2] = solution(p, q, f, a, b, tau,U1a, U2a, n);
Us = zeros(size(U1));
Us(end) = (gamma1.*U1(end)+betta1*U2(end))/(betta1+alpha1*U1(end));

[U] = reverseSolution(U1, U2, Us, n);


figure 
plot(x,u,x,U,'o');
legend("Exact solution","progonka");


function [U1, U2] = solution(p, q, f, a, b, tau,U1a, U2a, n)
  x = a:tau:b;
  U1(1) = U1a;
  U2(1) = U2a;
  for k = 1:n - 1
    U1(k + 1) = U1(k) + tau * (-U1(k)^2.*q(x(k))+U1(k).*p(x(k)) + 1);
    U2(k + 1) = U2(k) + tau * (-U1(k).*(U2(k).*q(x(k))+f(x(k))));
  end
end

function [U] = reverseSolution(U1, U2, U, n)
 tau = 0.01;
  for k = n:-1:2
   U(k-1) =U(k)-tau.*((U(k) - U2(k))./U1(k));
  end
end
