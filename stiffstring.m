%Frequency dependent losses not taken into account
%param assign
Fs=10000;
k=1/Fs;;
h=0.03;
gamma=0.05;
lambda=(gamma*k)/h;
kappa=.01;
mu=kappa*k/(h)^2;
hmin=sqrt((gamma^2*k^2)+sqrt((gamma^4*k^4)+16*(kappa^2*k^2)))
x=0:h:1;
t=0:k:2-k;

%raised cosine pulse params
c0=1;
x0=0.5;
xhw=0.1;

%state initialization
u=zeros(round((2/k)),round((1/h))+1);
u(1,:)=(c0/2)*(1+cos(pi*(x-x0)/xhw));
u(1,1:(x0-xhw)/h)=0;
u(1,(x0+xhw)/h:(1/h)+1)=0;
v=zeros(1,length(x));
v(1,round((x0-xhw)/h)+1:round(x0/h)-1)=.01;
v(1,round(x0/h):round((x0+xhw)/h)-1)=.01;
u(2,:)=u(1,:)+(k*v(1,:));

%recursion
a=size(u);
energy=zeros(1,a(1));
for i=2:a(1)
    for j=3:a(2)-2;
        u(i+1,j)=(2*u(i,j))-(u(i-1,j))+(lambda^2)*(u(i,j+1)-(2*u(i,j))+u(i,j-1))-(mu^2)*(u(i,j+2)-(4*u(i,j+1))+(6*u(i,j))-(4*u(i,j-1))+u(i,j-2));
    end
end

%energy calculation
kinetic=zeros(1,a(1));
potential=zeros(1,a(1));

for i = 2:a(1)-1
    kinetic(1,i)=(0.5)*((u(i,:)-u(i-1,:))/k)*((u(i,:)-u(i-1,:))/k)';
    potential(1,i)=0.5*gamma^2*(diff(u(i,:))/h)*(diff(u(i-1,:))/h)'+ 0.5*kappa^2*(1/(h)^4)*diff(u(i,:),2)*diff(u(i-1,:),2)';
end
energy=kinetic+potential;

%plots 

subplot(211)
plot(t,kinetic,t,potential,t,energy);
title('Energy plots')
legend('Numerical KE','Numerical PE','Total Numerical Energy')
xlabel('Time');
ylabel('Energy');

subplot(212)
plot(t(3:length(t)-3),(energy(3:length(energy)-3)-energy(1,3))/energy(1,3));
xlabel('Time')
ylabel('(E-E(3))/E(3)');
title('Energy variation')


