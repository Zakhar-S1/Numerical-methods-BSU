clear;
clc;

k=4;
a=2+0.3*k;

L=1;
h=0.1;
tau=a/2*h^2;
T=0.02;


[U,X,T]=gridSolution(h,tau,L,T,a);
surf(X,T,U');
title('Grid Method');
xlabel('X');
ylabel('T');

function y=fi(t)
k=4;
a=2+0.3*k;
y=exp(a*t);
end

function y=myu(t)
k=4;
a=2+0.3*k;
y=exp(a*(t-1));
end

function y=nyu(x)
k=4;
a=2+0.3*k;
y=exp(-a*x);
end

function [u,x,t]=gridSolution(h,tau,L,T,a)
    N=L/h;
    K=T/tau;
    
    %{
    %}

    for i=1:N+1
        x(i)=(i-1)*h;
        u(i,1)=fi(x(i));
    end
    
    for j=1:K+1
        t(j)=(j-1)*tau;
        u(1,j)=myu(t(j));
        u(N+1,j)=nyu(t(j));
    end
    s=1/a*tau/h^2;
    for j=1:K
        for i=2:N
        u(i,j+1)=s*u(i+1,j)+s*u(i-1,j)+(1-2*s)*u(i,j);
        end
    end
end

