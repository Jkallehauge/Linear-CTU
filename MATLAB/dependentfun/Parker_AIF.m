function Cb=Parker_AIF(x) % x is in minutes
A1=0.809;
A2=0.330;
T1=0.17046;
T2=0.365;
Sigma1=0.0563;
Sigma2=0.132;
alpha=1.050;
beta=0.1685;
s=38.078;
tau=0.483;

Cb1=(A1/(Sigma1*sqrt(2*pi)))*exp(-(x-T1).^2/(2*Sigma1^2));
Cb2=(A2/(Sigma2*sqrt(2*pi)))*exp(-(x-T2).^2/(2*Sigma2^2));
Cb3=alpha*exp(-beta*x)./(1+exp(-s*(x-tau)));
 Cb=Cb1+Cb2+Cb3;
% if plot_val
% figure;
% plot(t,Cb1,'-r');
% hold on;
% plot(t,Cb2,'-b');
% hold on;
% plot(t,Cb3,'-g');
% hold on;
% plot(t,Cb,'-k');
% end