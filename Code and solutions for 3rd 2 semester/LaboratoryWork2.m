function []= LaboratoryWork2()
tau = 0.1;
t0 = 1;
tk = 2;
b=0.5;
f = @(t,u)(1/t^2-u/t-b*u^2);
T = t0:tau:tk;
U(1) = 2;
[t,u]=ode45(f,[t0,tk],U(1));
for j=1:length(T)-1
    k1=tau*f(T(j), U(j));
    k2=tau*f(T(j)+tau/2, U(j)+k1(1)/2 );
    k3=tau*f(T(j)+tau/2, U(j)+k2(1)/2);
    k4=tau*f(T(j)+tau, U(j)+k3(1));
    dU = (k1+2*k2+2*k3+k4)/6;
    U(j+1)=U(j) + dU(1);
end 
plot(T,U,'r',t,u,'ob'); grid on;
