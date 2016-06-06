function [ gaussStruct  ] = gaussian_Prepare(prepareDataInfo, commonDataInfo)


%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% final version: from external file...
gaussStruct.modelVariability=1; %1=a different model for each cell, 2=same model for flagRegioMat
gaussStruct.typeOfModel=2; %2=REG 3=RF random forest, 4=ann
gaussStruct.pcaFlag=0; %0 means no PCA-noNorm, 1 means PCA-norm
gaussStruct.nPrec=5; %5 for PM, 2 for O3 (nox, voc), 1 for NO2 (nox)
gaussStruct.absDel=1; %absolute(0) or delta(1) values
gaussStruct.arealPoint=0; %0 means areal and point summed up, 1 means areal and point separated
gaussStruct.domain='ineris_7km_deliver_20150831';
gaussStruct.flagReg='ineris_7km_deliver_20150831';
gaussStruct.rad=[50 50 50 50 50]; %roughly 500km
gaussStruct.npat=[1 1 1 1 1];

gaussStruct.aqiFil='PM25_Y';%PM25, O3
gaussStruct.aqi='PM25';%PM25, O3
gaussStruct.nSc=27;
gaussStruct.Ide=[1:19];
gaussStruct.Val=[20:27];%1:26];
gaussStruct.rf=0 %window of cells of training varying F (1=one ring of cells used for training, surrounding the target cell0
gaussStruct.function_type=1; %1 means inverse distance 1/d^2, 2 means exponential
%moved to finalize section...
gaussStruct.nameTest=strcat('20150921_omega2_allDomain_',int2str(gaussStruct.rad(1)),'_rf',int2str(gaussStruct.rf),'tra1_19');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOAD THIS IF YOU ALREADY OPTIMIZED OMEGA
% load('code_optim_per_precursor_cell_rf3_inv/rf3_optim_sce50percRed_alpha.mat');
% load('./code_opt_prec_cel/resultOmegaFinal.mat');
% omega=omegaFinal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINITION OF THE FUNCTION USED FOR EMISSION AGGREGATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X, Y] = meshgrid(-gaussStruct.rad(1):1:gaussStruct.rad(1), -gaussStruct.rad(1):1:gaussStruct.rad(1));
switch gaussStruct.function_type
    case 1
        %load alpha & omega (from netcdf extra parameter)
        gaussStruct.coeff=2; %could vary from 0.5 to 5
        gaussStruct.sigma_x_gaus = NaN; %from 2 to rf/2 both sigma_x_gaus and sigma_x_gaus
        gaussStruct.sigma_y_gaus = NaN; %from 2 to rf/2 both sigma_x_gaus and sigma_x_gaus
        gaussStruct.theta=NaN;      %from 0 to pi
        gaussStruct.centr=gaussStruct.rad(1)+1;
        gaussStruct.F=1./(1+sqrt(X.^2+Y.^2).^gaussStruct.coeff);
    case 2
        %load alpha & omega (from netcdf extra parameter)
        A = 1; x0 = 0; y0 = 0;
        gaussStruct.coeff=NaN;
        gaussStruct.sigma_x_gaus = 6; %from 2 to rf/2 both sigma_x_gaus and sigma_x_gaus
        gaussStruct.sigma_y_gaus = 3; %from 2 to rf/2 both sigma_x_gaus and sigma_x_gaus
        gaussStruct.theta=0,      %from 0 to pi
        a = cos(theta)^2/2/sigma_x_gaus^2 + sin(theta)^2/2/sigma_y_gaus^2;
        b = -sin(2*theta)/4/sigma_x_gaus^2 + sin(2*theta)/4/sigma_y_gaus^2 ;
        c = sin(theta)^2/2/sigma_x_gaus^2 + cos(theta)^2/2/sigma_y_gaus^2;
        gaussStruct.F = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;
end

nameDirOut=strcat('./output/',gaussStruct.domain,'/',gaussStruct.aqiFil,...
    '/absDel_',int2str(gaussStruct.absDel),'/arealPoint_',int2str(gaussStruct.arealPoint),...
    '/',gaussStruct.nameTest,'/rf_',int2str(gaussStruct.rf),'-modTyp',int2str(gaussStruct.typeOfModel),...
    '-modVar',int2str(gaussStruct.modelVariability),'-pca',int2str(gaussStruct.pcaFlag),'/');
mkdir(nameDirOut);
gaussStruct.nameDirOut=nameDirOut;

end