function [f] = model_C_TU_expconv(x,par)
Fp=x(1); %Fp
Tp=x(2)/(x(3)+x(1)); %vp/PS
E=x(3)/(x(1)+x(3));%PS/(Fp+PS)

%IRF=(1-E)*exp(-par(:,1)/Tp)+E;
deltaT=diff(par(:,1));
convE=deltaT(1)*Fp*conv(par(:,2),ones(size(par(:,2)))*E);
f_=Fp*(1-E)*Tp*expconv(Tp,par(:,1),par(:,2))+convE(1:numel(par(:,2)))   ;

f =  f_ - par(:,3);