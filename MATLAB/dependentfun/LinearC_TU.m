function [Fp, Vp, PS, residual,best_delay,t_C_bestdelay,C_best_fit]=LinearC_TU(Cp,Ct,t,start,max_delay,plot_graph)
% Jesper Kallehauge 04/05-2015
%Denne funktion beregner den extended Tofts model udfra den Murase
%artiklen. Denne funktion tester en rï¿½kke forskellige delays imellem 
%bolus nï¿½r til arterien og nï¿½r vï¿½vet Ct starter med at oplades. Dem som 
%testes er imellem 0 og max_delay
%
%Input skal vï¿½re:
%
%  En plasma funktion                             Cp  (AIF)
%  En opladningskurve                             Ct  (tissue)
%  En tids vektor                                 t   (i sekunder)
%  Tidspunktet hvor bolus kommer til arterien     start
%  Max forskellen imellem tiden nï¿½r bolusen nï¿½r
%  til arterien og opladning til vï¿½vet sker       max_delay 
%  plot er boolean og bestemme om output skal 
%  vises grafisk
%
% eksempel pï¿½ input er 
%   Cp=squeeze(cimage(60,123,3,:));
%   Ct=squeeze(cimage(87,94,3,:));
%   TR=2.9; %sekunder
%   t=TR*(0:(size(Cp,1)-1));    %bestemmes i minutter
%   start=21;
%   max_delay=8;
%   plot_graph=1;
%   
%   [Ktrans, kep, Ve, Vp,residual,delay]=LinearToft(Cp,Ct,t,start,max_delay,plot_graph);

scanlength=numel(Cp);
best_delay=0;
Fp=0;
PS=0;
Vp=0;
residual=0;
if plot_graph
    figure;
    plot(t,Cp);
    legendstr{1}='AIF';
    hold on;
    plot(t,Ct);
    legendstr{2}='Tissue enhancement';
end    
for delay=0:max_delay
    clear A B C
    %shifting the relative start between AIF and Enhancement curve
     A=zeros((scanlength-delay-start+1),3,1);
     B=zeros(3,1); %vector med løsning
     C=zeros((scanlength-delay-start+1),1);
        
     intCp=cumtrapz(t(start:(scanlength-delay)),Cp(start:(scanlength-delay)));  %1st columns in the A matrix
     int2Cp=cumtrapz(t(start:(scanlength-delay)),intCp(1:(scanlength-delay-start+1)));  %1st columns in the A matrix
     Ctissue=-cumtrapz(t((start+delay):scanlength),Ct((start+delay):scanlength,:),1);%integral of the enhandementcurve;
       
     A(:,1)=Ctissue;
     A(:,2)=intCp;
     A(:,3)=int2Cp;
        
     C(1:(scanlength-(delay+start-1)),1)=Ct((start+delay):scanlength,:);%enhancementkurve

    B=A\C;
    %B=GaussJordanElimination(transpose(A)*A,transpose(A)*C)';
    C_temp=A*B;
    C_diff_sqrt(delay+1)=sqrt(sum((C-C_temp).^2))/numel(C);
    if plot_graph
        Fp=B(2)*60;   %bestemmes i minutter
        PS=(B(2)*B(3)/(B(1)*B(2)-B(3))*60);   %bestemmes i minutter
        Vp=(B(2)*B(2)/(B(1)*B(2)-B(3)));
        hold on;plot(t(start+delay:scanlength),C_temp);
        legendstr{delay+3}=['Delay=',num2str(delay),' Fp=',num2str(Fp),' PS=',num2str(PS),' Vp=',num2str(Vp)];
    end
end


[~,min_index]=min(C_diff_sqrt);
best_delay=min_index-1;
 if plot_graph
legend(legendstr);
title(['best delay: ', num2str(best_delay)] );
 end

clear A B C;
A=zeros((scanlength-best_delay-start+1),3,1);
B=zeros(3,1); %vector med løsning
C=zeros((scanlength-best_delay-start+1),1);

intCp=cumtrapz(t(start:(scanlength-best_delay)),Cp(start:(scanlength-best_delay)));  %1st columns in the A matrix
int2Cp=cumtrapz(t(start:(scanlength-best_delay)),intCp(1:(scanlength-best_delay-start+1)));  %1st columns in the A matrix
Ctissue=-cumtrapz(t((start+best_delay):scanlength),Ct((start+best_delay):scanlength,:),1);%integral of the enhandementcurve;

        
     A(:,1)=Ctissue;
     A(:,2)=intCp;
     A(:,3)=int2Cp;
        
C(1:(scanlength-(best_delay+start-1)),1)=Ct((start+best_delay):scanlength,:);%enhancementkurve
B=A\C;

    Fp=B(2)*60;   %bestemmes i minutter

    PS=(B(2)*B(3)/(B(1)*B(2)-B(3))*60);   %bestemmes i minutter

    Vp=(B(2)*B(2)/(B(1)*B(2)-B(3)));
    
residual=C_diff_sqrt(best_delay+1);

A_bestdelay=A;
B_bestdelay=B;

C_bestdelay=C;
C_best_fit=A_bestdelay*B_bestdelay;
t_C_bestdelay=t(start+best_delay:scanlength)/60;
if (plot_graph==1)
    figure;
    plot(t(start:scanlength)/60,Cp,'-r','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',10);
    hold on
    plot(t_C_bestdelay,C_bestdelay(:,1),'-g','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10);    hold on
    plot(t_C_bestdelay,C_best_fit,'-b','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',10);
    h = legend('AIF','enhancement curve', 'best fit',3);
    set(h,'Interpreter','none','Location','NorthEast');
    title(['optimum delay=',num2str(t(start+best_delay)-t(start)),' seconds with residual value =',num2str(residual)]);
    xlabel('Time (Minutes)'); ylabel('Concentration Gd [mmol]');
end;