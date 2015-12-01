clear;
h1=figure;
addpath([pwd,'\dependentfun']);
dataDir=[pwd,'\Data'];
options = optimset('Display','off','Algorithm','Trust-region-reflective','MaxFunEvals',25,'MaxIter',25,'TolX',1e-8,'TolFun',1e-8);
deltaT_high=0.01/60; %in minutes
startT_high=-20/60;% time before injection
duration=4; %in minutes
t_high=[startT_high:deltaT_high:round(duration/deltaT_high-1)*deltaT_high];
Hct=0.38;
vratio=0.7;
AIFscale=(1-vratio*Hct)/(1-Hct);
Cp_high=AIFscale*Parker_AIF(t_high);
t_high=t_high-startT_high;
CNR=10;
deltaT=2/60; %in minutes
downsample=deltaT/deltaT_high;

%%Plot simulation examples
for runs=1:3
    switch runs
        case 1
            subplot(2,3,1);
            Fp_ref=0.23;
            PS_ref=0.02;
            vp_ref=0.05;
        case 2
            subplot(2,3,2);
            Fp_ref=0.57;
            PS_ref=0.2;
            vp_ref=0.28;
        case 3
            subplot(2,3,3);
            Fp_ref=0.65;
            PS_ref=0.14;
            vp_ref=0.22;
    end
    
    Tp=vp_ref/((Fp_ref+PS_ref));
    E=PS_ref/(Fp_ref+PS_ref);
    C_clean_highE=deltaT_high*Fp_ref*conv(Cp_high,ones(size(Cp_high))*E);
    C_clean_high=Fp_ref*(1-E)*Tp*expconv(Tp,t_high,Cp_high)'+C_clean_highE(1:numel(Cp_high));
    
    t=t_high(1:downsample:numel(t_high));
    Cp=Cp_high(1:downsample:numel(t_high));
    
    clear A B ;
    
    stddev=max(C_clean_high)/CNR;
    noise = stddev*randn(size(C_clean_high));
    C_=C_clean_high+noise;
    C=C_(1:downsample:numel(t_high));
    
    [fitpar_NLLS,resnormNLLS,residual,exitflag,OUTPUT,LAMBDA,J] = lsqnonlin(@model_C_TU_expconv, [Fp_ref vp_ref PS_ref], [0 0 0], [inf inf inf], options, [t; Cp; C]');
    intC=cumtrapz(t/deltaT,C);
    intCp=cumtrapz(t/deltaT,Cp);
    int2Cp=cumtrapz(t/deltaT,intCp);
    
    A(:,1)=-intC;%
    A(:,2)=intCp;
    A(:,3)=int2Cp;
    B=A\C';
    LLS_residual=A*B-C';
    resnormLLS=sqrt(sum(LLS_residual.^2));
    
    fitpar_LLS(1)=B(2)/deltaT;
    fitpar_LLS(2)=(B(2)*B(2)/(B(1)*B(2)-B(3)));
    fitpar_LLS(3)=(B(2)*B(3)/(B(1)*B(2)-B(3)))/deltaT;
    figure(h1);
    
    Tp=fitpar_NLLS(2)/((fitpar_NLLS(1)+fitpar_NLLS(3)));
    E=fitpar_NLLS(3)/(fitpar_NLLS(1)+fitpar_NLLS(3));
    
    C_clean_highE=deltaT*fitpar_NLLS(1)*conv(Cp,ones(size(Cp))*E);
    C_clean_low_fit=fitpar_NLLS(1)*(1-E)*Tp*expconv(Tp,t,Cp)'+C_clean_highE(1:numel(t));
    
    hold on;
    plot(t,C,'ok');
    hold on;
    plot(t_high,C_clean_high,'-k','LineWidth',2);
    
    hold on;
    plot(t,C_clean_low_fit,'-b','LineWidth',2);
    
    hold on;
    plot(t,A*B,'-r','LineWidth',2);
    
    legend('synthetic data','ideal fit','NLLS','LLS');
    title(['L^2_{NLLS}=',num2str(sqrt(sum(residual.^2))),'; L^2_{LLS}=',num2str(sqrt(sum((LLS_residual).^2)))]);
    xlabel('time [min]');
    ylabel('Indicator concentration [mM]');
end
%--------------plot patient data
A=importdata([dataDir,'\Test_patient_raw_data']);
t=A(4:end,1)/60;
deltaT=mean(diff(t));
Cp=A(4:end,2);
start=1;
max_delay=8;
iter=0;
for N=[3 1600 4500]
    iter=iter+1;
    
    
    clear C AMat B delay_map Fp_map PS_map Vp_map Residual_map chi2_map B_1CM;
    %C-TU LLS
    C=A(4:end,N);
    [Fp_LLS, Vp_LLS, PS_LLS, resnormLLS,best_delay,t_C_bestdelay_out,C_best_fit]=LinearC_TU(Cp,C,t*60,start,max_delay,0);
    
    clear C;
    C=A(4+best_delay:end,N);
    %C-TU NLLS
    [fitpar,resnormNLLS,residual,exitflag,OUTPUT,LAMBDA,J] = lsqnonlin(@model_C_TU_expconv, [fitpar_LLS(1) fitpar_LLS(2) fitpar_LLS(3)], [0 0 0], [inf 100 inf], options, [t(1:(end-best_delay)) Cp(1:(end-best_delay)) C]);
    fitpar_NLLS=[fitpar(1) fitpar(2) fitpar(3)];
    Tp=fitpar_NLLS(2)/((fitpar_NLLS(1)+fitpar_NLLS(3))); %vp/(PS+Fp)
    E=fitpar_NLLS(3)/(fitpar_NLLS(1)+fitpar_NLLS(3));%PS/(Fp+PS)
    C_clean_highE=deltaT*fitpar_NLLS(1)*conv(Cp,ones(size(Cp))*E);
    C_clean_low_fit=fitpar_NLLS(1)*(1-E)*Tp*expconv(Tp,t,Cp)'+C_clean_highE(1:numel(t))';
    
    figure(h1);
    subplot(2,3,3+iter);
    
    plot(t,C,'ok'); hold on; plot(t,C_clean_low_fit,'-b','LineWidth',2);hold on; plot(t,C_best_fit,'-r','LineWidth',2);
    legend('measured data','NLLS','LLS');
    title(['L^2_{NLLS}=',num2str(sqrt(sum(residual.^2))),' L^2_{LLS}=',num2str(resnormLLS*numel(C_best_fit))]);
    xlabel('time [min]');
    ylabel('Indicator concentration [mM]');
end


