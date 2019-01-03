%params assign
Fs=1000;
k=1/Fs;;
h=0.001;
gamma=1;
lambda=(gamma*k)/h
x=0:h:1;
t=0:k:2-k;;
sigma=0.5;

%raised cosine pulse parameters
c0=1;
x0=0.5;
xhw=0.1;

%initialization
u=zeros((2/k),(1/h)+1);
u(1,:)=(c0/2)*(1+cos(pi*(x-x0)/xhw));
u(1,1:(x0-xhw)/h)=0;
u(1,(x0+xhw)/h:(1/h)+1)=0;
v=zeros(1,length(x));
v(1,((x0-xhw)/h)+1:(x0/h)-1)=.01;
v(1,(x0/h):((x0+xhw)/h)-1)=.01;
u(2,:)=u(1,:)+(k*v(1,:));

%recursion
a=size(u)
energy=zeros(1,a(1));
for i=2:a(1)-1
    for j=2:a(2)-1;
        u(i+1,j)=(2/(1+(sigma*k)))*((1-lambda^2)*u(i,j) + (lambda^2/2)*(u(i,j-1)+u(i,j+1)))-((1-(sigma*k))/(1+(sigma*k)))*u(i-1,j);
    end
end

%energy calculation

kinetic=zeros(1,a(1));
potential=zeros(1,a(1));

for i = 2:a(1)-1
    kinetic(1,i)=(0.5)*((u(i,:)-u(i-1,:))/k)*((u(i,:)-u(i-1,:))/k)';
    potential(1,i)=0.5*gamma^2*(diff(u(i,:))/h)*(diff(u(i-1,:))/h)';
end
energy=kinetic+potential;

%plot
figure()

plot(t,kinetic,t,potential,t,energy)
xlabel('Time')
ylabel('Energies')
title('Energy Plots');
legend('Numerical KE','Numerical PE','Total Numerical Energy')

hold 

figure()

subplot(151)
plot(x,u(1,:))
xlabel('Time')
ylabel('Amplitude')
title('t=0s')

subplot(152)
plot(x,u(300,:))
xlabel('Time')
ylabel('Amplitude')
title('t=0.3s')

subplot(153)
plot(x,u(600,:))
xlabel('Time')
ylabel('Amplitude')
title('t=0.6s')

subplot(154)
plot(x,u(900,:))
xlabel('Time')
ylabel('Amplitude')
title('t=0.9s')

subplot(155)
plot(x,u(1200,:))
xlabel('Time')
ylabel('Amplitude')
title('t=1.2s')



