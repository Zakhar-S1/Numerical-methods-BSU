function [] = LaboratoryWork4()
f=@(t,u)-14*u(1);
T=10;
T2=9.8;
h=0.01;
h2=0.14;
ts=0:h:T;
t = 0:h:T;
u0=1; 
U1=zeros(1,T/h+1);
U2=zeros(1,T/h+1);
F=zeros(1,T/h+1);
F2=zeros(1,T/h+1);
U1(1)=u0;
U2(1)=u0;
F(1)=f(0,u0);
F2(1)=f(0,u0);

for n=1:3
   k1 = f(t(n), U1(n)).*h;
   k2 = f(t(n) + h / 2, U1(n) + k1 / 2).*h;
   k3 = f(t(n) + h / 2, U1(n) + k2 / 2).*h;
   k4 = f(t(n) + h, U1(n) + k3).*h;
   U1(n + 1) = U1( n) + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4);  
   F( n + 1) = f(t(n + 1), U1(n + 1));
end
U2=U1;
F2=F;
B = h * [-9/24; 37/24; -59/24; 55/24];

for n=4:T/h
    U1(n + 1) = U1(n) + F(n-3 : n) * B;
    F( n + 1) = f(t(n + 1), U1(n + 1));
end

B2 = h2 * [-9/24; 37/24; -59/24; 55/24];
for n=4:T2/h2
    U2(n + 1) = U2(n) + F2(n-3 : n) * B2;
    F2( n + 1) = f(t(n + 1), U2(n + 1));
end
 
accurateF=exp(-14*ts);
 
[t,y] = ode45(f, [0 T], u0,odeset('RelTol',1e-5)); 
startTransformY = repelem(y,4);
endTransformY = startTransformY(1:1001,:);

figure 
semilogy(abs(endTransformY'-accurateF))
legend('ode45 error')
figure 
semilogy(abs(U1-accurateF))
legend('Adams error')
figure 
plot(ts,U1(1,:),'-bo',ts, accurateF,'-r*',t,y(:,1),'-k')
legend('Adams','Accurate','ode45')
figure 
plot(ts,U2(1,:),'-bo',ts, accurateF,'-r*')
legend('Unstable','Accurate')