clear;
dataDir=[pwd,'\Data'];
addpath([pwd,'\dependentfun']);
A=importdata([dataDir,'\Test_patient_raw_data']);
t=A(4:end,1)/60;
deltaT=mean(diff(t));
Cp=A(4:end,2);

load([dataDir,'\Modelfitcomparison.mat']);    
VOI=zeros(176,176,20);
timepoints=numel(A(4:end,1));
xind=squeeze(A(1,3:end));
yind=squeeze(A(2,3:end));
zind=squeeze(A(3,3:end));
DCE=zeros(176,176,20,timepoints);
for Npoints=1:numel(xind)
   DCE(xind(Npoints),yind(Npoints),zind(Npoints),:)=A(4:end,2+Npoints);
end    

index=sub2ind(size(VOI),xind(:),yind(:),zind(:));
VOI_CNR=zeros(size(VOI));
VOI_CNR(index)=CNR;
VOI(index)=1;

VOI_C_TU_Fp_LLS=zeros(size(VOI));
VOI_C_TU_PS_LLS=zeros(size(VOI));
VOI_C_TU_Vp_LLS=zeros(size(VOI));
VOI_C_TU_resnorm_LLS=zeros(size(VOI));

VOI_C_TU_Fp_LLS(index)=Fitpars_LLS(1,:);
VOI_C_TU_PS_LLS(index)=Fitpars_LLS(3,:);
VOI_C_TU_Vp_LLS(index)=Fitpars_LLS(2,:);
temp=Residual_LLS;
VOI_C_TU_resnorm_LLS(index)=sqrt(nansum(temp.^2));

VOI_C_TU_Fp_NLLS=zeros(size(VOI));
VOI_C_TU_PS_NLLS=zeros(size(VOI));
VOI_C_TU_Vp_NLLS=zeros(size(VOI));
VOI_C_TU_resnorm_NLLS=zeros(size(VOI));

VOI_C_TU_Fp_NLLS(index)=Fitpars_NLLS(1,:);
VOI_C_TU_PS_NLLS(index)=Fitpars_NLLS(3,:);
VOI_C_TU_Vp_NLLS(index)=Fitpars_NLLS(2,:);
temp=Residual_NLLS;
VOI_C_TU_resnorm_NLLS(index)=sqrt(nansum(temp.^2));

VOI_C_TU_E=zeros(size(VOI));
VOI_C_TU_E(index)=Eigen_vals_LLS(3,:)./(Eigen_vals_LLS(1,:).*Eigen_vals_LLS(2,:));

[val, ~]=nanmin(VOI);
STATS = regionprops(VOI);

slicenum=round(STATS.Centroid(3));
startY=round(STATS.BoundingBox(1))-5;
startX=round(STATS.BoundingBox(2))-5;
deltaX=round(max(STATS.BoundingBox(4:5)))+10;
deltaY=round(max(STATS.BoundingBox(4:5)))+10;

VOI_slice=zeros(size(VOI));
VOI_slice(:,:,slicenum)=VOI(:,:,slicenum);
index_slice=find(VOI_slice);

figure;
VOI_tmp=zeros(size(VOI));
VOI_tmp(index(~isnan(VOI_C_TU_Fp_LLS(index))))=1;
VOI=VOI_tmp;

ax_col=3;
ax_row=4;

maxlim_row1=0.25*max(VOI_C_TU_Fp_LLS(index_slice));
maxlim_row4=0.15*max(VOI_C_TU_PS_LLS(index_slice));
colormap jet;
h1=subaxis(ax_col,ax_row,5,'Spacing',0.001,'Padding',0,'Margin',0);
him=imagesc(squeeze(VOI_C_TU_Fp_LLS(startX:startX+deltaX,startY:startY+deltaY,slicenum)),[0 maxlim_row1]);
s1Pos = get(h1,'position');
colorbar('location','north');
set(h1,'position',s1Pos);
set(him,'AlphaData',squeeze(VOI(startX:startX+deltaX,startY:startY+deltaY,slicenum)));
axis off;
axis equal;
ylabel('Tofts model');

h2=subaxis(ax_col,ax_row,6,'Spacing',0.001,'Padding',0,'Margin',0);
him=imagesc(squeeze(VOI_C_TU_PS_LLS(startX:startX+deltaX,startY:startY+deltaY,slicenum)),[0 maxlim_row4]);
s2Pos = get(h2,'position');
colorbar('location','north');
set(h2,'position',s2Pos);
set(him,'AlphaData',squeeze(VOI(startX:startX+deltaX,startY:startY+deltaY,slicenum)));
axis off;
axis equal;

h3=subaxis(ax_col,ax_row,7,'Spacing',0.001,'Padding',0,'Margin',0);
him=imagesc(squeeze(VOI_C_TU_Vp_LLS(startX:startX+deltaX,startY:startY+deltaY,slicenum)),[0 1]);
s3Pos = get(h3,'position');
colorbar('location','north');
set(h3,'position',s3Pos);
set(him,'AlphaData',squeeze(VOI(startX:startX+deltaX,startY:startY+deltaY,slicenum)));
axis off;
axis equal;

h4=subaxis(ax_col,ax_row,8,'Spacing',0.001,'Padding',0,'Margin',0);
him=imagesc(squeeze(VOI_C_TU_resnorm_LLS(startX:startX+deltaX,startY:startY+deltaY,slicenum)),[0 1.5]);
s11Pos = get(h4,'position');
colorbar('location','north');
set(h4,'position',s11Pos);
set(him,'AlphaData',squeeze(VOI(startX:startX+deltaX,startY:startY+deltaY,slicenum)));
axis off;
axis equal;

subaxis(ax_col,ax_row,9,'Spacing',0.001,'Padding',0,'Margin',0);
him=imagesc(squeeze(VOI_C_TU_Fp_NLLS(startX:startX+deltaX,startY:startY+deltaY,slicenum)),[0 maxlim_row1]);
set(him,'AlphaData',squeeze(VOI(startX:startX+deltaX,startY:startY+deltaY,slicenum)));
ylabel('Extended Tofts model');
axis off;
axis equal;

subaxis(ax_col,ax_row,10,'Spacing',0.001,'Padding',0,'Margin',0);
him=imagesc(squeeze(VOI_C_TU_PS_NLLS(startX:startX+deltaX,startY:startY+deltaY,slicenum)),[0 maxlim_row4]);
set(him,'AlphaData',squeeze(VOI(startX:startX+deltaX,startY:startY+deltaY,slicenum)));
axis off;
axis equal;

h23=subaxis(ax_col,ax_row,11,'Spacing',0.001,'Padding',0,'Margin',0);
him=imagesc(squeeze(VOI_C_TU_Vp_NLLS(startX:startX+deltaX,startY:startY+deltaY,slicenum)),[0 1]);
s11Pos = get(h23,'position');
set(h23,'position',s11Pos);
set(him,'AlphaData',squeeze(VOI(startX:startX+deltaX,startY:startY+deltaY,slicenum)));
axis off;
axis equal;

subaxis(ax_col,ax_row,12,'Spacing',0.001,'Padding',0,'Margin',0);
him=imagesc(squeeze(VOI_C_TU_resnorm_NLLS(startX:startX+deltaX,startY:startY+deltaY,slicenum)),[0 1.5]);
set(him,'AlphaData',squeeze(VOI(startX:startX+deltaX,startY:startY+deltaY,slicenum)));
axis off;
axis equal;

h = text(2,4, 'LLS'); set(h, 'rotation', 90,'fontsize', 16);
h = text(2,4, 'NLLS'); set(h, 'rotation', 90,'fontsize', 16);
h = text(2,4, 'Fp [min^{-1}]'); set(h, 'rotation', 0,'fontsize', 16);
h = text(2,4, 'PS [min^{-1}]'); set(h, 'rotation', 0,'fontsize', 16);
h = text(2,4, 'Vp'); set(h, 'rotation', 0,'fontsize', 16);
h = text(2,4, '\chi^2/n'); set(h, 'rotation', 0,'fontsize', 16);