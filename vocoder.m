% Vocoder

% Silenzio
W1=[65.2552  -57.3959   -0.0913    0.0255   -0.0000;
  -57.3959   55.5802    0.0476   -0.1737    0.0546;
   -0.0913    0.0476    0.0013    0.0039   -0.0014;
    0.0255   -0.1737    0.0039    0.0210   -0.0071;
   -0.0000    0.0546   -0.0014   -0.0071    0.0026];
m1=[49.4225  -51.3214   -0.0396    0.1248   -0.0393];
W11=inv(W1);
% Voiced
W2=[12.3484   -5.9115    0.0642    0.4776   -0.8825;
   -5.9115   21.8853   -0.6446   -0.2708    0.7828;
    0.0642   -0.6446    0.0227    0.0027   -0.0180;
    0.4776   -0.2708    0.0027    0.0319   -0.0425;
   -0.8825    0.7828   -0.0180   -0.0425    0.0969];
m2=[10.2900  -33.2312    0.9078    0.3549   -1.1069];
W21=inv(W2);
% Unvoiced
W3=[143.4694  -18.4661   -2.5371   -3.8762    3.0234;
  -18.4661   43.3894   -0.3503   -0.9416    0.8914;
   -2.5371   -0.3503    0.0601    0.0951   -0.0723;
   -3.8762   -0.9416    0.0951    0.1763   -0.1295;
    3.0234    0.8914   -0.0723   -0.1295    0.1225];
m3=[64.4750  -37.3261   -0.4644   -0.6587    0.6755];
W31=inv(W3);


% Settaggio variabili
fs=10000;       % frequenza di campionamento
fp=300;         % lunghezza del frame per decisione pitch
fv=100;         % lunghezza del frame per decisione v/u/s
durata=2;       % durata della registrazione in secondi
p=14;           % ordine di linear predictive coding
zi=20;          % numero di zeri preposti al segnale

fprintf('Premi un tasto per registrare %g secondi di voce...', durata); 
pause;
s=wavrecord(durata*fs,fs);
a=130*2*pi;
b=200*2*pi;
T=1/fs;
ss=filter([1 -2 1],[1 -2*exp(-a*T)*cos(b*T) exp(-2*a*T)],s);    % filtro high-pass 200 Hz

ss=filter([1 -0.9],1,ss);
ss=[zeros(zi,1); ss];

% Decisione v/u/s

nfr=floor(durata*fs/fv);   % numero dei frame
% inizializzazione vettori
E=[];
D=[];
Dec=[];
lpccoeff=[];
Alpha=[];

% Ciclo principale di decisone v/u/s
for i=0:nfr-1
    n=i*fv+1+zi;
    zc=zerocrossing(ss(n:(n-1+fv)));
    es=logenergy(ss(n:(n-1+fv)));
    if n==1
        prec=ss(1);
    else prec=ss(n-1);
    end
    c=autoc1(ss(n:(n-1+fv)),prec);
    [alpha,e]=lpc(ss(n:(n-1+fv)),p-1);
    lpccoeff=[lpccoeff; alpha];
    ep=prederr(es,ss,n,fv,alpha);
    D=[D; [zc es c ep alpha(2)]];
    Alpha=[Alpha; alpha];
    y1=D(i+1,:)-m1;
    y2=D(i+1,:)-m2;
    y3=D(i+1,:)-m3;
    d1=y1*W11*y1';
    d2=y2*W21*y2';
    d3=y3*W31*y3';
    den=d1*d2+d2*d3+d1*d3;
    P1=(d2*d3)/den;
    P2=(d1*d3)/den;
    P3=(d1*d2)/den;
    [mi,dec]=minim([d1 d2 d3]);
    Dec=[Dec; dec P1 P2 P3];
    E=[E; e];
end
Dec=Dec(20:end,:);
E=E(20:end,:);
Dec(:,1)=mediansmooth(Dec(:,1),6);

% Decisione pitch

for i=1:length(Dec)   
    if Dec(i)==2
        n=i*fv+1+zi+20*fv;
        f=minim([n-1+fp length(ss)]);
        m=allamdf(ss(n:f),fp);
        p(i)=m;
    else p(i)=0;
    end
end

% Sintesi

fprintf('Premi un tasto per cominciare la sintesi'); 
pause
synth=[];
% p=mediansmooth(p,10);
for i=1:length(Dec)
    G=sqrt(E(i));
    switch Dec(i)
        case 1
            synth=[synth zeros(1,fv)];
        case 2
            u=zeros(1,p(i));
            u(1)=1;
            j=1;
            while length(u)<fv   % creazione eccitazione quasi-perodica
                u=[u u];
            end
            u=u(1:fv);
            sy=filter(G,Alpha(i,:),u);
            synth=[synth sy];
        case 3
            u=wgn(fv,1,1);
            u=u/max(u);
            sy=filter(G,Alpha(i,:),u');
            synth=[synth sy];
    end
end
for i=1:length(synth)
    if isreal(synth(i))==0
        synth(i)=real(synth(i));
    end
end
synth=filter(1,[1 -0.9],synth);
wavplay(synth*2,fs)
