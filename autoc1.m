function c=autoc1(x,s)
% Primo coefficiente di autocorrelazione del frame x
% è necessario fornire l'ultimo campione s del frame precedente

l=length(x);
num=x(1)*s;
d2=s^2;
d1=x(1)^2;
for i=2:l
    num=num+x(i)*x(i-1);
    d1=d1+x(i)^2;
    d2=d2+x(i-1)^2;
end
c=num/sqrt(d1*d2);
        
