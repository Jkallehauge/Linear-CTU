function f=expconv(T, time, a)

    if T==0
        f=a;
    else
       n=numel(time);
       f=zeros(n,1);
       x=(time(2:n)-time(1:n-1))/T;
       da=(a(2:n)-a(1:n-1))./x;
       E=exp(-x);
       E0=1-E;
       E1=x-E0;
       
       add=a(1:n-1).*E0+da.*E1;
       
       for i=1:n-1
           f(i+1)=E(i)*f(i)+add(i);
       end
    end   