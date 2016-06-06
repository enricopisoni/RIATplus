%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIAL SETTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; warning('off')
tic
%FIX DATA
%sherpa test
prepareDataFile=strcat('./input/prepareGaussRegression_1.txt');
finalizeJobFile=strcat('./input/finalizeRegression.txt');

%riat test
%prepareDataFile=strcat('./input/prepareGaussRegression_1.txt');
%finalizeJobFile=strcat('./input/finalizeRegression.txt');

%prepareDataFile='';
%finalizeJobFile='';

%prepareDataInfo=loadPrepareDataSetting(prepareDataFile);
prepareDataInfo=load_PrepareDataSetting(prepareDataFile);
finalizeJobInfo=load_FinalizeJobSetting(finalizeJobFile);

%[precursors, nx, ny, totalCells]=loadPrecursors();
%indicators=loadIndicators();
load './testdata/Precursor'
load './testdata/Indicator'

%all scenario, all precursors ('til nprec parameter)...
%precursors=Prec(:,:,1,1);
%indicators=Indic(:,:,:);
%first scenario, first precursor...
nPrec=1;
nScen=1;
precursors=Prec(:,:,nPrec,nScen);
%first scenario, first precursor...
indicators=Indic(:,:,nScen);

totCells=nx*ny;

%optimizerValues=0;
%optimizerCondition='';

intermediateResult=interface_prepareData(precursors, indicators, x, y, nx, ny, totCells, prepareDataInfo, optimizerValues, optimizerCondition);

finalInfo=interface_finalizeJob(finalizeJobInfo, prepareDataInfo, intermediateResult, precursors, indicators, x, y, nx, ny, totCells, optimizerValues, optimizerCondition);

    
%display results
%creating scatter
%[corr_reg,mse_reg]=CreateScatter(IndicBC,Indic(:,:,finalInfo.iSc),finalInfo.finalResult,finalInfo.flagRegioMat,finalInfo.iSc,nx,ny,intermediateResult.nameDirOut,intermediateResult.aqi,intermediateResult.absDel);
[corr_reg,mse_reg]=CreateScatter(IndicBC,Indic(:,:),finalInfo.finalResult,finalInfo.flagRegioMat,1,nx,ny,intermediateResult.nameDirOut,intermediateResult.aqi,intermediateResult.absDel);
    
%loop on scenarios
   
%creating maps
%CreateMap(IndicBC,Indic(:,:,finalInfo.iSc),finalInfo.finalResult,finalInfo.flagRegioMat,x,y,finalInfo.iSc,intermediateResult.nameDirOut,intermediateResult.aqi,intermediateResult.absDel,intermediateResult.flagReg,intermediateResult.domain);
CreateMap(IndicBC,Indic(:,:,1),finalInfo.finalResult,finalInfo.flagRegioMat,x,y,1,intermediateResult.nameDirOut,intermediateResult.aqi,intermediateResult.absDel,intermediateResult.flagReg,intermediateResult.domain);
close all

