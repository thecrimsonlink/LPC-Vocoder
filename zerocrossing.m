function y=zerocrossing(x)
% Numero di zero crossing del frame x
y=0;
l=length(x);
for i=1:l-1
    if sign(x(i))~=sign(x(i+1))
        y=y+1;
    end
end
