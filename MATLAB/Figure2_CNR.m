% Initialize
clear;
addpath([pwd,'\dependentfun']);
h = waitbar(0,'Please wait running through CNR analysis...');
numrepeat=1000;
CNRarray=[2:0.2:5 5.5:0.5:10 12 15:5:40];

CNR_index=0;
step_val=0;
steps_total=numel(CNRarray)*numrepeat;
NLLS_result=zeros(numel(CNRarray),numrepeat,3);
LLS_result=zeros(numel(CNRarray),numrepeat,3);
exit_result=zeros(numel(CNRarray),numrepeat);

Accuracy_Fp_NLLS=zeros(numel(CNRarray),1);
Accuracy_Fp_LLS=zeros(numel(CNRarray),1);
Accuracy_vp_NLLS=zeros(numel(CNRarray),1);
Accuracy_vp_LLS=zeros(numel(CNRarray),1);
Accuracy_PS_NLLS=zeros(numel(CNRarray),1);
Accuracy_PS_LLS=zeros(numel(CNRarray),1);

precision_Fp_NLLS=zeros(numel(CNRarray),1);
precision_Fp_LLS=zeros(numel(CNRarray),1);
precision_vp_NLLS=zeros(numel(CNRarray),1);
precision_vp_LLS=zeros(numel(CNRarray),1);
precision_PS_NLLS=zeros(numel(CNRarray),1);
precision_PS_LLS=zeros(numel(CNRarray),1);

L2_LLS=zeros(numrepeat,numel(CNRarray));
L2_NLLS=zeros(numrepeat,numel(CNRarray));

DownsampleRes=2/60;% in minutes
randt0=1;

% Literature values from Kallehauge JF, Tanderup K, Duan C, Haack S, Pedersen EM,
% Lindegaard JC, et al. Tracer kinetic model selection for dynamic contrast-enhanced
% magnetic resonance imaging of locally advanced cervical cancer. Acta Oncol (Madr).
% 2014;53(8):1064–72.
Fp_ref=0.57;
PS_ref=0.2;
vp_ref=0.28;


%% generate high resolution AIF based on litterature values from
% Parker GJM, Roberts C, Macdonald A, Buonaccorsi GA, Cheung S, Buckley DL,
% et al. Experimentally-derived functional form for a population-averaged
% high-temporal-resolution arterial input function for dynamic contrast-enhanced MRI.
% MRM. 2006;56(5):993–1000.

DownsampleRes_high=0.01/60; %in minutes
startT_high=-20/60; %in minutes
duration=4; %in minutes
t_high=[startT_high:DownsampleRes_high:round(duration/DownsampleRes_high-1)*DownsampleRes_high];
Hct=0.38;
vratio=0.7;
AIFscale=(1-vratio*Hct)/(1-Hct);
Cp_high=AIFscale*Parker_AIF(t_high);

downsample=DownsampleRes/DownsampleRes_high;

%% generate high resolution tissue enhancement
t_high=t_high-startT_high;

Tp=vp_ref/((Fp_ref+PS_ref));
E=PS_ref/(Fp_ref+PS_ref);
C_clean_highE=DownsampleRes_high*Fp_ref*conv(Cp_high,ones(size(Cp_high))*E);
C_clean_high=(Fp_ref*(1-E)*Tp*expconv(Tp,t_high,Cp_high)+C_clean_highE(1:numel(Cp_high))')';

options = optimset('Display','off','Algorithm','Trust-region-reflective','MaxFunEvals',25,'MaxIter',25,'TolX',1e-8,'TolFun',1e-8);

for CNR=CNRarray
    plot_true=1;
    CNR_index=CNR_index+1;
    for repeat=1:numrepeat;
        noise=zeros(size(C_clean_high));        
        if randt0 % 1 if we vary the onset time within one sampling period
            t0=round(rand(1)*(downsample-1))+1; %Randomly shifting the onset time within one sampling period.
        else
            t0=1;
        end
        t_low=t_high(1:downsample:numel(t_high)-(t0-1));
        t=t_low;
        step_val=step_val+1;
        waitbar(step_val / steps_total, h ,['calculating Fp=',num2str(Fp_ref),' PS=',num2str(PS_ref), ' Vp=',num2str(vp_ref),' CNR=',num2str(CNR)]);
        if exist('CNR','var')
            stddev=max(C_clean_high)/CNR;
            noise = stddev*randn(size(C_clean_high));
            C_clean=C_clean_high(t0:downsample:numel(C_clean_high));
            C_=C_clean_high(t0:downsample:numel(C_clean_high));
            Cp_=Cp_high(t0:downsample:numel(C_clean_high));
            noise = stddev*randn(size(C_));
            C=C_+noise;
            noise = stddev*randn(size(Cp_));
            Cp=Cp_+noise;
        else
            C_clean=C_clean_high(t0:downsample:numel(C_clean_high));
            C=C_clean_high(t0:downsample:numel(C_clean_high));
            Cp=Cp_high(t0:downsample:numel(C_clean_high));
        end
        
        %NLLS
        [fitpar,resnorm,residual,exitflag,OUTPUT,LAMBDA,J] = lsqnonlin(@model_C_TU_expconv, [Fp_ref vp_ref PS_ref], [0 0 0], [inf inf inf], options, [t; Cp; C]');
        fitpar_NLLS=[fitpar(1) fitpar(2) fitpar(3)];
        L2_NLLS(repeat,CNR_index)=sqrt((sum(residual.^2)));
        
        %LLS
        intC=cumtrapz(t/DownsampleRes,C);
        intCp=cumtrapz(t/DownsampleRes,Cp);
        int2Cp=cumtrapz(t/DownsampleRes,intCp);
        
        A(:,1)=-intC; %calculating matrix A (equation 12)
        A(:,2)=intCp;
        A(:,3)=int2Cp;
        B=A\C'; % Matrix inversion
        
        fitpar_LLS(1)=B(2)/DownsampleRes; % estimating parameters via equation 9
        fitpar_LLS(2)=(B(2)*B(2)/(B(1)*B(2)-B(3)));
        fitpar_LLS(3)=(B(2)*B(3)/(B(1)*B(2)-B(3)))/DownsampleRes;
        
        LLS_residual=A*B-C';
        L2_LLS(repeat,CNR_index)=sqrt(sum(LLS_residual.^2));
        
        NLLS_result(CNR_index,repeat,:)=fitpar_NLLS';
        LLS_result(CNR_index,repeat,:)=fitpar_LLS';
        exit_result(CNR_index,repeat)=exitflag;
    end
    Accuracy_Fp_NLLS(CNR_index)=100*mean((squeeze(NLLS_result(CNR_index,:,1)))-ones(numrepeat,1)'*Fp_ref)/Fp_ref;
    Accuracy_Fp_LLS(CNR_index)=100*mean((squeeze(LLS_result(CNR_index,:,1)))-ones(numrepeat,1)'*Fp_ref)/Fp_ref;
    Accuracy_vp_NLLS(CNR_index)=100*mean((squeeze(NLLS_result(CNR_index,:,2)))-ones(numrepeat,1)'*vp_ref)/vp_ref;
    Accuracy_vp_LLS(CNR_index)=100*mean((squeeze(LLS_result(CNR_index,:,2))-ones(numrepeat,1)'*vp_ref)/vp_ref);
    Accuracy_PS_NLLS(CNR_index)=100*mean((squeeze(NLLS_result(CNR_index,:,3)))-ones(numrepeat,1)'*PS_ref)/PS_ref;
    Accuracy_PS_LLS(CNR_index)=100*mean((squeeze(LLS_result(CNR_index,:,3)))-ones(numrepeat,1)'*PS_ref)/PS_ref;
    
    precision_Fp_NLLS(CNR_index)=100*std(squeeze(NLLS_result(CNR_index,:,1)))/Fp_ref;
    precision_Fp_LLS(CNR_index)=100*std(squeeze(LLS_result(CNR_index,:,1)))/Fp_ref;
    precision_vp_NLLS(CNR_index)=100*std(squeeze(NLLS_result(CNR_index,:,2)))/vp_ref;
    precision_vp_LLS(CNR_index)=100*std(squeeze(LLS_result(CNR_index,:,2)))/vp_ref;
    precision_PS_NLLS(CNR_index)=100*std(squeeze(NLLS_result(CNR_index,:,3)))/PS_ref;
    precision_PS_LLS(CNR_index)=100*std(squeeze(LLS_result(CNR_index,:,3)))/PS_ref;
end
figure;
subplot(2,2,1);
[l,p] = boundedline(CNRarray, Accuracy_Fp_LLS, precision_Fp_LLS, '-r', CNRarray, Accuracy_Fp_NLLS,precision_Fp_NLLS, '-b','alpha');

hold on;
l3=plot([min(CNRarray) max(CNRarray)],[0 0],'-k'); l4=plot([10 10],[-50 50],'--k');
title(['Fp= ',num2str(Fp_ref),' [min^{-1}]']);
xlabel('CNR');ylabel('Accuracy [%]');
BL=legend([l(1), p(1),l(2),p(2), l3],'LLS (Accuracy)','LLS (Precision)','NLLS (Accuracy)','NLLS (Precision)','true value');
PatchInLegend = findobj(BL, 'type', 'patch');
set(PatchInLegend, 'facea', 0.5);
ylim([-50 50]);xlim([0 45]);

subplot(2,2,2);
[l,p] = boundedline(CNRarray, Accuracy_vp_LLS, precision_vp_LLS, '-r', CNRarray, Accuracy_vp_NLLS,precision_vp_NLLS, '-b','alpha');
hold on;
l3=plot([min(CNRarray) max(CNRarray)],[0 0],'-k'); l4=plot([10 10],[-50 50],'--k');
BL=legend([l(1), p(1),l(2),p(2), l3],'LLS (Accuracy)','LLS (Precision)','NLLS (Accuracy)','NLLS (Precision)','true value');
PatchInLegend = findobj(BL, 'type', 'patch');
set(PatchInLegend, 'facea', 0.5);
title(['Vp= ',num2str(vp_ref)]);
ylim([-50 50]);xlim([0 45]);
xlabel('CNR');ylabel('Accuracy [%]');

subplot(2,2,3);
[l,p] = boundedline(CNRarray, Accuracy_PS_LLS, precision_PS_LLS, '-r', CNRarray, Accuracy_PS_NLLS, precision_PS_NLLS, '-b','alpha');
hold on;
l3=plot([min(CNRarray) max(CNRarray)],[0 0],'-k');
l4=plot([10 10],[-50 50],'--k');
BL=legend([l(1), p(1),l(2),p(2), l3],'LLS (Accuracy)','LLS (Precision)','NLLS (Accuracy)','NLLS (Precision)','true value');
PatchInLegend = findobj(BL, 'type', 'patch');
set(PatchInLegend, 'facea', 0.5);
title(['PS= ',num2str(PS_ref),' [min^{-1}]']);
ylim([-50 50]);xlim([0 45]);
xlabel('CNR');ylabel('Accuracy [%]');

e_LLS = prctile(L2_LLS,97.5)-prctile(L2_LLS,2.5);
e_NLLS = prctile(L2_NLLS,97.5)-prctile(L2_NLLS,2.5);
subplot(2,2,4);
[l,p] = boundedline(CNRarray, mean(L2_LLS), e_LLS, '-r', CNRarray, mean(L2_NLLS), e_NLLS, '-b','alpha');
title(['Fit quality']);
BL=legend([l(1), p(1),l(2),p(2)],'LLS','95% confidence interval','NLLS','95% confidence interval');
PatchInLegend = findobj(BL, 'type', 'patch');
set(PatchInLegend, 'facea', 0.5);
xlabel('CNR');ylabel('L2-norm');
ylim([0 5]);xlim([0 45]);