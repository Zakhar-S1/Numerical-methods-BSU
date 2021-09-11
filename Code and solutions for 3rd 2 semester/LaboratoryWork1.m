function []= LaboratoryWork1()
a=1;
b=2;
n=50;
h=(b-a)/n;
xSt=1:1/n:2;
uSt=xSt.*exp(1).^(1-1./xSt);
F=@(x,y)(y/(x^2)+exp(1)^((x-1)/x));
x(1)=a;
u(1)=1;
for i=1:n-1
    x(i+1)=a+i*h;
    u(i+1)=u(i)+h*F(x(i),u(i));
end 
plot(x,u,'r',xSt,uSt,'.b'); grid on;