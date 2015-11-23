% This program is meant to test the correctness of the fits for the two
% separate models while computing kappa and rho
close all
clear all
addpath('/myOpt/Matlab_AddOns/afni_matlab/matlab');
I = 32;
J = 55;
K = 8;
C = 32;

% FILE NAMES
% ----------
PRJDIR='/data/SFIMJGC/TALK_fMRIClassME/PrcsData/SBJ02/';
BetaFile  =strcat(PRJDIR,'/D03_MEICA/SBJ02_S02Run10.chComp.EXTRA.Beta',num2str(C,'%03d'),'.nii');
cR2File   =strcat(PRJDIR,'/D03_MEICA/SBJ02_S02Run10.chComp.cR2.nii');
FR2File   =strcat(PRJDIR,'/D03_MEICA/SBJ02_S02Run10.chComp.FR2.nii');
lR2File   =strcat(PRJDIR,'/D03_MEICA/SBJ02_S02Run10.chComp.EXTRA.R2Fit',num2str(C,'%03d'),'.nii');

cS0File   =strcat(PRJDIR,'/D03_MEICA/SBJ02_S02Run10.chComp.cS0.nii');
FS0File   =strcat(PRJDIR,'/D03_MEICA/SBJ02_S02Run10.chComp.FS0.nii');
lS0File   =strcat(PRJDIR,'/D03_MEICA/SBJ02_S02Run10.chComp.EXTRA.S0Fit',num2str(C,'%03d'),'.nii');

EchosFile=strcat(PRJDIR,'/D00_OriginalData/SBJ02_S02Run10_Echoes.1D');

I = 32; J = 55; K = 8; C = 32;

% Load the data
[~,Betas,  bInfo,~] = BrikLoad(BetaFile);
[~,Py_cR2s,~    ,~] = BrikLoad(cR2File);
[~,Py_FR2s,~    ,~] = BrikLoad(FR2File);
[~,Py_lR2s,~    ,~] = BrikLoad(lR2File);
[~,Py_cS0s,~    ,~] = BrikLoad(cS0File);
[~,Py_FS0s,~    ,~] = BrikLoad(FS0File);
[~,Py_lS0s,~    ,~] = BrikLoad(lS0File);


TE     = load(EchosFile)  
Beta   = squeeze(Betas(I+1,J+1,K+1,:)); 
Py_cR2 = Py_cR2s(I+1,J+1,K+1,C+1);
Py_FR2 = Py_FR2s(I+1,J+1,K+1,C+1);
Py_lR2 = squeeze(Py_lR2s(I+1,J+1,K+1,:));

Py_cS0 = Py_cS0s(I+1,J+1,K+1,C+1);
Py_FS0 = Py_FS0s(I+1,J+1,K+1,C+1);
Py_lS0 = squeeze(Py_lS0s(I+1,J+1,K+1,:));

ML_cR2 = lscov((TE/mean(TE))',Beta);
ML_cS0 = lscov(ones(length(TE),1),Beta);




disp(['Echo Times     :',num2str(TE)])
disp(['Betas          :',     num2str(Beta')])
disp(['Py cR2(F) | ML :',     num2str(Py_cR2),' (',num2str(Py_FR2),') | ',num2str(ML_cR2)])
disp(['Py R2 Fit      :',     num2str(Py_lR2')])
disp(['Py cS0(F) | ML :',     num2str(Py_cS0),' (',num2str(Py_FS0),') | ',num2str(ML_cS0)])
disp(['Py S0 Fit      :',     num2str(Py_lS0')])




plot(TE/mean(TE),Beta,'k-o')
hold on;
plot(TE/mean(TE),Py_lS0,'r');
plot(TE/mean(TE),Py_lR2,'g');
varTotal = sum(Beta.^2);
varExpS0 = sum(Py_lS0.^2);
varExpR2 = sum(Py_lR2.^2);

% Actual Linear Fit
[a,b]=polyfit((TE/mean(TE))',Beta,1);
LNFit = ((TE/mean(TE))*a(1)) + a(2);
plot(TE/mean(TE),LNFit,'c--');

SSM_LN = sum((LNFit-mean(Beta')).^2);
SSE_LN = sum((Beta'-LNFit).^2);
SST_LN = sum((Beta'-mean(Beta')).^2);
DFM_LN = length(a)-1;
DFE_LN = length(TE)-length(a);
DFT_LN = length(TE)-1;
MSM_LN = SSM_LN/DFM_LN;
MSE_LN = SSE_LN/DFE_LN;
MST_LN = SST_LN/DFT_LN;
F_LN = MSM_LN/MSE_LN; % Agrees with R
pF_LN = 1-cdf('F',F_LN,1,1); % Agrees with R