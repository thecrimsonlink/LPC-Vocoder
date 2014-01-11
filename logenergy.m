function e=logenergy(s)
% Logaritmo dell'energia del frame x

l=length(s);
c=0;
for i=1:l
    c=c+(s(i)^2);
end
e=10*log10((1e-6)+c/l);
