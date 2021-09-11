function []= LaboratoryWork3()
tau = 0.1;
t0 = 0;
tn = 0.5;
k=14;
a = 1+0.2*k;
uf1 = @(t,u1,u2)(-2*a*t*u1^2+u2^2-t^2-1);
uf2 = @(t,u1,u2)(1/(a*u2^2)-u1-t/u2);
T = t0:tau:tn;
u1s(1) = 1/a;
u2s(1) = 1;
for i=1:length(T)-1
    k1=tau*uf1(T(i), u1s(i),u2s(i));
    m1=tau*uf2(T(i), u1s(i),u2s(i));
    k2=tau*uf1(T(i)+tau/2, u1s(i)+k1(1)/2,u2s(i)+m1(1)/2 );
    m2=tau*uf2(T(i)+tau/2, u1s(i)+k1(1)/2,u2s(i)+m1(1)/2 );
    k3=tau*uf1(T(i)+tau/2, u1s(i)+k2(1)/2,u2s(i)+m2(1)/2);
    m3=tau*uf2(T(i)+tau/2, u1s(i)+k2(1)/2,u2s(i)+m2(1)/2);
    k4=tau*uf1(T(i)+tau, u1s(i)+k3(1),u2s(i)+m3(1));
    m4=tau*uf2(T(i)+tau, u1s(i)+k3(1),u2s(i)+m3(1));
    u1s(i+1)=u1s(i) + (k1+2*k2+2*k3+k4)/6;
    u2s(i+1)=u2s(i) + (m1+2*m2+2*m3+m4)/6;
end 
F=@(t,U)[-2*a*t*U(1)^2+U(2)^2-t^2-1;1/(a*U(2)^2)-U(1)-t/U(2)];
[t1,sol]=ode45(F,[t0,tn],[1/a;1]);
plot(T,u1s,T,u2s, t1,sol(:,1),'g*',t1,sol(:,2),'y*');
legend('R-K1','R-K2', 'Exact solution1','Exact solution1')
grid on;
