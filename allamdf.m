function m=allamdf(x,l)
% Trova il minimo dell'AMDF per ogni frame di lunghezza l del segnale x

L=length(x);
N=floor(L/l);  % numero di frame calcolati
for i=0:N-1
    y=x((i*l+1):((i+1)*l));
    a=amdf(y,l);
    m(i+1)=minloc(a);
end

function a=amdf(x,l)
% Average Magnitude Difference Function del segnale x, calcolata in l punti

M=length(x);
for k=1:l
    c=0;
    for m=1+k:M
        c=c+abs(x(m)-x(m-k));
    end
    a(k)=c;
end

function i=minloc(x)
% Trova il primo minimo locale della sequenza x

i=1;
for k=1:length(x)-2
    if x(k)>x(k+1)
        if x(k+1)<x(k+2)
            i=k+1;
            break
        end
    end
end
