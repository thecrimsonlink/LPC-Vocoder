function ep=prederr(es,ss,n,fv,alpha)
% Energia logaritmica dell'errore di previsione
% es=energia logaritmica
% ss=segnale
% fr=numero del frame considerato
% fv=lunghezza del frame
% alpha= coefficienti lpc

c=-phi(ss,n,fv,0);
for i=1:length(alpha)
    c=c+alpha(i)*phi(ss,n,fv,i);
end
ep=es-10*(log10((1e-6)+abs(c)));
    

function f=phi(ss,n,fv,k)
c=0;
for i=0:fv-1
    c=c+ss(n+i)*ss(n+i-k);
end
f=c/fv;
