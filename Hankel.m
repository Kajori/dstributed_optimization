%H amtrixfunction H=Hankel(y,m,n)if size(y,1)>size(y,2)    y=y.';end    if length(y)<m+n+1    error('large dimension, not enough data');else    H=zeros(m,n);    for j=1:m        H(j,:)=y(j:j+n-1);    endend