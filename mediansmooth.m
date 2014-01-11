function y=mediansmooth(x,n)
[a,b]=size(x);
if a~=1
    x=x';
end
y(1)=x(1);
for i=2:n
    z=zeros(1,n-i);
    m=[z*x(1) x(1:i)];
    med=median(m);
    y(i)=med;
end
for i=n+1:length(x)
    y(i)=median(x(i-n:i));
end
