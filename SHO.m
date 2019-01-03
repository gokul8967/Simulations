Fs=44100;
F0=13000;
W0=2*pi*F0;
k=1/Fs;
%u=zeros(1,400);
u(1)=0;
u(2)=0.0001;
for n=1:500;
    u(n+2)=((2-((k^2)*(W0^2)))*u(n+1)) - u(n);
end
n=0:length(u)-1;
plot(n/Fs,u)
xlabel('Time')
ylabel('Amplitude of oscillator')
