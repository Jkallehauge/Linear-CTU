function [f] = model_TM_v2_expconv(x,par)
T=x(2)/x(1);

f_ = x(1)*T*expconv(T,par(:,1),par(:,2));
f =  f_ - par(:,3);
