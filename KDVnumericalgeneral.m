clc
clear all
M=1.4;
%c=1.4;
A=.48;
B=4.84E-3;
epsilon=0.5*sqrt(A.*M./B);
mu=epsilon.*M.*A;
D=1.5;
d= sqrt((M+epsilon.*D)./mu);
a=-2;
b=2;
c=5;
%d=10;
nx=1000;
nt=1000;
%M=1.2;
dx=(b-a)/nx;
dt=c/nt;
x=a:dx:b;
t=0:dt:c;
u = zeros(nt+1, nx+1);
alpha=d.*dx/2;
H1 = 1./(4*cos(alpha));
H2 = -(cos(alpha)-1)./(2*sin(alpha)).^2;
H3 = ((b.^2)*cos(alpha))./(3*sin(alpha)).^2;
for n = 1:nt+1
    for i= 1:nx+1
        %u(i,n) = exp(x(i:n)+t(i:n));
        %u(i+1,n+1) = 1;
        u(i,n) =(sech(-D+epsilon.*x(i)-mu.*t(n))).^(2);
        %u(i,n) =(sech(epsilon.*x(i)-mu*t(n)+D)).^(2);
        %u(i,n) = exp(1j*d*(x(i)+c*t(n)))+exp(-1j*d*(x(i)+c*t(n)));
    end
end
%u(1,n) = 0;
%u(i+1,n) = 0;
% %u(i:n)=exp(1j.*b(0+i.*dx+M.*(0+n.*dt)))+exp(-1j.*b(0+i.*dx+M.*(0+n.*dt)));
for i = 2:nx
    Ux1= u(i+1,:)-u(i-1,:);
    Uxx2=u(i+1,:)-2.*u(i,:)+u(i-1,:);
end
for i= 1:nx-1
    Ux2= u(i+2,:)-u(i,:);
    Uxx1=u(i+2,:)-2.*u(i+1,:)+u(i,:);
    F1=(epsilon./2).*(u(i+1,:)-Ux1.*H1+Uxx1.*H2).^2+mu.*((b.^2).*Ux1.*H1+Uxx1.*H3);
end
for i= 3:nx+1
    Ux3= u(i,:)-u(i-1,:);
    Uxx3=u(i,:)-2.*u(i-1,:)+u(i-2,:);
end
for i = 2:nx-1
    F2=(epsilon./2).*(u(i,:)+Ux2.*H1+Uxx2.^2)+mu.*(-(b.^2).*Ux2.*H1+Uxx2.*H3);
    F3=(epsilon./2).*(u(i,:)-Ux2.*H1+Uxx2.*H2).^2+mu.*((b.^2).*Ux2.*H1+Uxx2.*H3);
end
for i=2:nx+1
    F4=(epsilon./2).*(u(i-1,:)+Ux3.*H1+Uxx3.*H2).^2+mu.*(-(b.^2).*Ux3.*H1+Uxx3.*H3);
end
% Uxx2=u(i+1:n)-2.*u(i:n)+u(i-1:n);
% Uxx1=u(i+2:n)-2.*u(i+1:n)+u(i:n);
% Uxx3=u(i:n)-2.*u(i-1:n)+u(i-2:n);



F11=0.5.*(F1+F2);
F22=0.5.*(F3+F4);
u=u(i,:)+(dt./dx).*(F11-F22);
plot(x,(u),'o')
hold on
disp([x,u])