function [D,W,m]=vustrainer
% Vocoder

% Settaggio variabili
fs=10000;       % frequenza di campionamento
fp=300;         % lunghezza del frame per decisione pitch
fv=100;         % lunghezza del frame per decisione v/u/s
durata=2;       % durata della registrazione in secondi
p=14;           % ordine di linear predictive coding
zi=20;          % numero di zeri preposti al segnale
% hh=[0.2057   -0.1271   -0.1006   -0.0837   -0.0735   -0.0675   -0.0641   -0.0624   -0.0619 -0.0618    0.9382   -0.0618   -0.0619   -0.0624   -0.0641   -0.0675   -0.0735   -0.0837 -0.1006   -0.1271    0.2057];

fprintf('Premi un tasto per registrare %g secondi di voce...', durata); 
pause;
s=wavrecord(durata*fs,fs);
a=130*2*pi;
b=200*2*pi;
T=1e-4;
ss=filter([1 -2 1],[1 -2*exp(-a*T)*cos(b*T) exp(-2*a*T)],s);
ss=filter([1 -0.9],1,ss);
ss=[zeros(zi,1); ss];

nfr=floor(durata*fs/fv);
D=[];
for i=0:nfr-1
    n=i*fv+1+zi;
    zc=zerocrossing(ss(n:(n-1+fv)));
    es=logenergy(ss(n:(n-1+fv)));
    if n==1
        prec=ss(1);
    else prec=ss(n-1);
    end
    c=autoc1(ss(n:(n-1+fv)),prec);
    alpha=lpc(ss(n:(n-1+fv)),p-1);
    ep=prederr(es,ss,n,fv,alpha);
    D=[D; [zc es c ep alpha(2)]];   
end
D(:,1)=mediansmooth(D(:,1),10);
D(:,2)=mediansmooth(D(:,2),10);
D(:,3)=mediansmooth(D(:,3),10);
D(:,4)=mediansmooth(D(:,4),10);
D(:,5)=mediansmooth(D(:,5),10);

% Calcolo media e matrice di covarianza

[N,c]=size(D);
d=sum(D);
m=d/N;
c=zeros(5);
for i=1:N
    x=D(i,:);
    c=c+x'*x;
end
W=c/N - m'*m;
figure(1)
plot(D(15:end,1))
figure(2)
plot(D(15:end,2))
figure(3)
plot(D(15:end,3))
figure(4)
plot(D(15:end,4))
figure(5)
plot(D(15:end,5))
figure(6)
plot(s(15:end))
