%% Initialize
clear;
addpath([pwd,'\dependentfun']);
options = optimset('Display','off','Algorithm','Trust-region-reflective','MaxFunEvals',100,'MaxIter',100,'TolX',1e-8,'TolFun',1e-8);
h = waitbar(0,'Please wait running through CNR analysis...');

Resarray=[0.05:0.05:10]/60; %timeresolution in minutes
randt0=1;
if randt0
    numrepeat=100;
else
    numrepeat=1;
end

LLS_timer=zeros(numel(Resarray),numrepeat,3);
NLLS_timer=zeros(numel(Resarray),numrepeat,3);

Accuracy_Fp_NLLS=zeros(numel(Resarray),1);
Accuracy_Fp_LLS=zeros(numel(Resarray),1);
Accuracy_vp_NLLS=zeros(numel(Resarray),1);
Accuracy_vp_LLS=zeros(numel(Resarray),1);
Accuracy_PS_NLLS=zeros(numel(Resarray),1);
Accuracy_PS_LLS=zeros(numel(Resarray),1);

for runs=1:3
    switch runs
        case 1
            Fp_ref=0.23;
            PS_ref=0.02;
            vp_ref=0.05;
        case 2    
            Fp_ref=0.57;
            PS_ref=0.2;
            vp_ref=0.28;
        case 3
            Fp_ref=0.65;
            PS_ref=0.14;
            vp_ref=0.22;         
    end

    
%% generate high resolution AIF based on litterature values from
% Parker GJM, Roberts C, Macdonald A, Buonaccorsi GA, Cheung S, Buckley DL,
% et al. Experimentally-derived functional form for a population-averaged
% high-temporal-resolution arterial input function for dynamic contrast-enhanced MRI.
% MRM. 2006;56(5):993–1000.

    deltaT_high=0.01/60; %in minutes
    startT_high=-20/60;% 
    duration=4; %in minutes
    t_high=[startT_high:deltaT_high:round(duration/deltaT_high-1)*deltaT_high];
    Hct=0.38;
    vratio=0.7;
    AIFscale=(1-vratio*Hct)/(1-Hct);
    Cp_high=AIFscale*Parker_AIF(t_high);

    
%% generate high resolution tissue enhancement    
    t_high=t_high-startT_high;
    
    Tp=vp_ref/((Fp_ref+PS_ref)); %vp/(PS+Fp)
    E=PS_ref/(Fp_ref+PS_ref);%PS/(Fp+PS)
    C_clean_highE=deltaT_high*Fp_ref*conv(Cp_high,ones(size(Cp_high))*E);
    C_clean_high=Fp_ref*(1-E)*Tp*expconv(Tp,t_high,Cp_high)'+C_clean_highE(1:numel(Cp_high));
             
    step_val=0;  
    steps_total=numel(Resarray)*numrepeat;
    
    NLLS_result=zeros(numel(Resarray),numrepeat,3);
    LLS_result=zeros(numel(Resarray),numrepeat,3);
    
    Res_index=0;
    for Res=Resarray
        Res_index=Res_index+1;
        downsample=round(Res/deltaT_high);
        Cp_orig=Cp_high(1:downsample:numel(t_high));

        for repeat=[1:1:numrepeat]
            if randt0
                t0=round(rand(1)*(downsample-1))+1;
            else
                t0=1;
            end    
            t_low=t_high(1:downsample:numel(t_high)-(t0-1));
            t=t_low;
            clear A B;
            step_val=step_val+1;
            waitbar(step_val / steps_total, h ,['calculating Fp=',num2str(Fp_ref),' PS=',num2str(PS_ref), ' Vp=',num2str(vp_ref),' Res=',num2str(Res)]);
            if exist('CNR','var')
                stddev=max(C_clean_high)/CNR;             
                C_=C_clean_high(t0:downsample:numel(C_clean_high));
                Cp_=Cp_high(t0:downsample:numel(C_clean_high));
                noise = stddev*randn(size(C_));
                C=C_+noise;
                noise = stddev*randn(size(Cp_));
                Cp=Cp_+noise;               
            else
                C=C_clean_high(t0:downsample:numel(C_clean_high));
                Cp=Cp_high(t0:downsample:numel(C_clean_high));
            end
            %NLLS                 
            tstart_NLLS=tic;
            [fitpar,resnorm,residual,exitflag,OUTPUT,LAMBDA,J] = lsqnonlin(@model_C_TU_expconv, [Fp_ref vp_ref PS_ref], [0 0 0], [inf inf inf], options, [t; Cp; C]');
            NLLS_timer(Res_index,repeat,runs)=toc(tstart_NLLS);
            fitpar_NLLS=fitpar;
            C_T=C';
            Cp_T=Cp';
            tstart_LLS=tic;
            [param]=mexLinCTU(C_T,Cp_T);

%           intC=cumtrapz(t/Res,C);
%           intCp=cumtrapz(t/Res,Cp);
%           int2Cp=cumtrapz(t/Res,intCp);
%           A(:,1)=-intC;%
%           A(:,2)=intCp;
%           A(:,3)=int2Cp;
%           B=A\C';
            
            fitpar_LLS(1)=param(1)/Res;
            fitpar_LLS(2)=param(2);
            fitpar_LLS(3)=param(3)/Res;
            LLS_timer(Res_index,repeat,runs)=toc(tstart_LLS);
            
            NLLS_result(Res_index,repeat,:)=fitpar_NLLS';
            LLS_result(Res_index,repeat,:)=fitpar_LLS';
        end
        Accuracy_Fp_NLLS(Res_index)=100*mean((squeeze(NLLS_result(Res_index,:,1)))-ones(numrepeat,1)'*Fp_ref)/Fp_ref;
        Accuracy_Fp_LLS(Res_index)=100*mean((squeeze(LLS_result(Res_index,:,1)))-ones(numrepeat,1)'*Fp_ref)/Fp_ref;
        Accuracy_vp_NLLS(Res_index)=100*mean((squeeze(NLLS_result(Res_index,:,2)))-ones(numrepeat,1)'*vp_ref)/vp_ref;
        Accuracy_vp_LLS(Res_index)=100*mean((squeeze(LLS_result(Res_index,:,2))-ones(numrepeat,1)'*vp_ref)/vp_ref);
        Accuracy_PS_NLLS(Res_index)=100*mean((squeeze(NLLS_result(Res_index,:,3)))-ones(numrepeat,1)'*PS_ref)/PS_ref;
        Accuracy_PS_LLS(Res_index)=100*mean((squeeze(LLS_result(Res_index,:,3)))-ones(numrepeat,1)'*PS_ref)/PS_ref;
    end
    
    
    
    if runs==1
        h1=figure;
    else
        figure(h1);
    end
    subplot(3,3,3*runs-2);
    plot(Resarray,Accuracy_Fp_NLLS,'-b','LineWidth',2);
    hold on;
    plot(Resarray,Accuracy_Fp_LLS,'-r','LineWidth',2);
    title(['Fp= ',num2str(Fp_ref),' [min^{-1}]'],'fontweight','bold','fontsize',12);
    legend('NLLS','LLS');
    
    xlabel('Temporal resolution [s]');
    ylabel('Accuracy [%]');
    ylim([-10 10]);
    
    subplot(3,3,3*runs-1);
    plot(Resarray,Accuracy_vp_NLLS,'-b','LineWidth',2);
    hold on;
    plot(Resarray,Accuracy_vp_LLS,'-r','LineWidth',2);
    title(['Vp= ',num2str(vp_ref)],'fontweight','bold','fontsize',12);
    legend('NLLS','LLS');
    xlabel('Temporal resolution [s]');
    ylabel('Accuracy [%]');
    ylim([-10 10]);
    
    subplot(3,3,3*runs);
    plot(Resarray,Accuracy_PS_NLLS,'-b','LineWidth',2);
    hold on;
    plot(Resarray,Accuracy_PS_LLS,'-r','LineWidth',2);
    title(['PS= ',num2str(PS_ref),' [min^{-1}]'],'fontweight','bold','fontsize',12);
    legend('NLLS','LLS');
    xlabel('Temporal resolution [s]');
    ylabel('Accuracy [%]');
    ylim([-10 10]);
    

end

figure;
for runs=1:3
    switch runs
        case 1
            hold on;
            plot(Resarray*60,mean(squeeze(NLLS_timer(:,:,runs))')./mean(squeeze(LLS_timer(:,:,runs))'),'-k','LineWidth',2);    
            Fp_ref=0.23;
            PS_ref=0.02;
            vp_ref=0.05;
            legd{runs}=['Fp= ',num2str(Fp_ref),' [min^{-1}], Vp= ',num2str(vp_ref),' PS= ',num2str(PS_ref),' [min^{-1}]'];
        case 2    
            hold on;
            plot(Resarray*60,mean(squeeze(NLLS_timer(:,:,runs))')./mean(squeeze(LLS_timer(:,:,runs))'),'-r','LineWidth',2);    
            Fp_ref=0.57;
            PS_ref=0.2;
            vp_ref=0.28;
            legd{runs}=['Fp= ',num2str(Fp_ref),' [min^{-1}], Vp= ',num2str(vp_ref),' PS= ',num2str(PS_ref),' [min^{-1}]'];
        case 3
            hold on;
            plot(Resarray*60,mean(squeeze(NLLS_timer(:,:,runs))')./mean(squeeze(LLS_timer(:,:,runs))'),'-b','LineWidth',2);
            Fp_ref=0.65;
            PS_ref=0.14;
            vp_ref=0.22;         
            legd{runs}=['Fp= ',num2str(Fp_ref),' [min^{-1}], Vp= ',num2str(vp_ref),' PS= ',num2str(PS_ref),' [min^{-1}]'];
    end
        
end
legend(legd,'fontweight','bold','fontsize',10);
    title('NLLS/LLS');
    xlabel('Temporal resolution [s]');
    ylabel('Speed up');
    ylim([0 1000]);
