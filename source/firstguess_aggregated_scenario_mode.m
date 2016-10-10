function [aqi_per_cell]=firstguess_aggregated_scenario_mode(emiTMP, aggregationInfo, commonDataInfo, aqiIndex, seasonIndex)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANNS CONSIDERED
%NN=nnSuperSet(ii).nnSet(jj);
%save('C:\data\work\projects\riat\RiatPlus-v6.1b\datasave\sa_start_fix', 'emiTMP','NN');
precNames=aggregationInfo.extraInfo.pollutantList;%{'NOx';'NMVOC';'NH3';'PM10';'PM25';'SOx'};
[ncel, nx, ny]=calcCellNo('latlon', commonDataInfo);
icells=ncel;
%coordinate=commonDataInfo.domainData.data(:,1:2);
%stepsize=round(max(abs(coordinate(2,1)-coordinate(1,1)),abs(coordinate(2,2)-coordinate(1,2))));
%nxny=round((max(coordinate)-min(coordinate))/stepsize+1);
areal_point=commonDataInfo.optim_flags.areal_point;
flag_optim_dom=commonDataInfo.domainInfo.flag_optim_dom;
pathAR=commonDataInfo.dirs.pathAR;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AREAL CASE
intermediateResult=0;
x=0;
y=0;
indiciMAT=0;
splitResult=0;
commonDataInfo.extraInfo=aggregationInfo.extraInfo;
[s1_NOX]=firstguess_Finalize(aggregationInfo, intermediateResult, commonDataInfo, emiTMP(:,1), 1, ...
    x, y, nx, ny, precNames(1), indiciMAT, splitResult);
[s1_VOC]=firstguess_Finalize(aggregationInfo, intermediateResult, commonDataInfo, emiTMP(:,2), 2, ...
    x, y, nx, ny, precNames(2), indiciMAT, splitResult);
[s1_NH3]=firstguess_Finalize(aggregationInfo, intermediateResult, commonDataInfo, emiTMP(:,3), 3, ...
    x, y, nx, ny, precNames(3), indiciMAT, splitResult);
[s1_PM10]=firstguess_Finalize(aggregationInfo, intermediateResult, commonDataInfo, emiTMP(:,4), 4, ...
    x, y, nx, ny, precNames(4), indiciMAT, splitResult);
[s1_PM25]=firstguess_Finalize(aggregationInfo, intermediateResult, commonDataInfo, emiTMP(:,5), 5, ...
    x, y, nx, ny, precNames(5), indiciMAT, splitResult);
[s1_SO2]=firstguess_Finalize(aggregationInfo, intermediateResult, commonDataInfo, emiTMP(:,6), 6, ...
    x, y, nx, ny, precNames(6), indiciMAT, splitResult);
NH3_all=s1_NH3.resGrid;
%save('C:\data\work\projects\riat\RiatPlus-v6.1b\datasave\sa_NH3_all_fix', 'NH3_all');
NOX_all=s1_NOX.resGrid;
PM10_all=s1_PM10.resGrid;
PM25_all=s1_PM25.resGrid;
SO2_all=s1_SO2.resGrid;
VOC_all=s1_VOC.resGrid;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONSIDER CASE WITH 6 OR 2 PERCURSOR EMISSIONS AS INPUT,
%AREAL CASE. ALWAYS THERE ARE 4 QUADRANTS TO CONSIDER WIND
%DIRETCIONS
if areal_point==0
    %read net and select precursors
    emissioni=[];
    for i_prec=1:size(precNames,1) %areal
        if strcmp(precNames(i_prec,1),'NH3')==1
            emissioni=[emissioni, NH3_all];
        elseif strcmp(precNames(i_prec,1),'NOx')==1
            emissioni=[emissioni, NOX_all];
        elseif strcmp(precNames(i_prec,1),'PM10')==1
            emissioni=[emissioni, PM10_all];
        elseif strcmp(precNames(i_prec,1),'PM25')==1
            emissioni=[emissioni, PM25_all];
        elseif strcmp(precNames(i_prec,1),'SOx')==1
            emissioni=[emissioni, SO2_all];
        elseif strcmp(precNames(i_prec,1),'NMVOC')==1
            emissioni=[emissioni, VOC_all];
        end
    end
    
    
    %             if (size(NN.net,1)==24 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==24))% ANNlinear, 6 input
    %                 emissioni=[NH3_all,NOX_all,PM10_all,PM25_all,SO2_all,VOC_all];
    %             elseif (size(NN.net,1)==8 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==8))%ANNlinear, 2 input
    %                 emissioni=[NOX_all,VOC_all];
    %             end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %IF ALSO POINT SOURCES
elseif areal_point==1
    %point
    [s1_NOXp]=firstguess_Finalize(aggregationInfo, intermediateResult, commonDataInfo, emiTMP(:,7), 1, ...
        x, y, nx, ny, precNames(1), indiciMAT, splitResult);
    [s1_VOCp]=firstguess_Finalize(aggregationInfo, intermediateResult, commonDataInfo, emiTMP(:,8), 2, ...
        x, y, nx, ny, precNames(2), indiciMAT, splitResult);
    [s1_NH3p]=firstguess_Finalize(aggregationInfo, intermediateResult, commonDataInfo, emiTMP(:,9), 3, ...
        x, y, nx, ny, precNames(3), indiciMAT, splitResult);
    [s1_PM10p]=firstguess_Finalize(aggregationInfo, intermediateResult, commonDataInfo, emiTMP(:,10), 4, ...
        x, y, nx, ny, precNames(4), indiciMAT, splitResult);
    [s1_PM25p]=firstguess_Finalize(aggregationInfo, intermediateResult, commonDataInfo, emiTMP(:,11), 5, ...
        x, y, nx, ny, precNames(5), indiciMAT, splitResult);
    [s1_SO2p]=firstguess_Finalize(aggregationInfo, intermediateResult, commonDataInfo, emiTMP(:,12), 6, ...
        x, y, nx, ny, precNames(6), indiciMAT, splitResult);
    NH3_allp=[s1_NH3p];
    %save('C:\data\work\projects\riat\RiatPlus-v6.1b\datasave\sa_NH3_allp_fix', 'NH3_allp');
    NOX_allp=s1_NOXp.resGrid;
    PM10_allp=s1_PM10p.resGrid;
    PM25_allp=s1_PM25p.resGrid;
    SO2_allp=s1_SO2p.resGrid;
    VOC_allp=s1_VOCp.resGrid;
    
    
    
    %create input structure, to be used in the ANNs
    %read net and select precursors
    emissioni=[];
    for i_prec=1:size(precNames,1) %areal
        if strcmp(precNames(i_prec,1),'NH3')==1
            emissioni=[emissioni, NH3_all];
        elseif strcmp(precNames(i_prec,1),'NOx')==1
            emissioni=[emissioni, NOX_all];
        elseif strcmp(precNames(i_prec,1),'PM10')==1
            emissioni=[emissioni, PM10_all];
        elseif strcmp(precNames(i_prec,1),'PM25')==1
            emissioni=[emissioni, PM25_all];
        elseif strcmp(precNames(i_prec,1),'SOx')==1
            emissioni=[emissioni, SO2_all];
        elseif strcmp(precNames(i_prec,1),'NMVOC')==1
            emissioni=[emissioni, VOC_all];
        end
    end
    for i_prec=1:size(precNames,1) %point
        
        if strcmp(precNames(i_prec,1),'NH3')==1
            emissioni=[emissioni, NH3_allp];
        elseif strcmp(precNames(i_prec,1),'NOx')==1
            emissioni=[emissioni, NOX_allp];
        elseif strcmp(precNames(i_prec,1),'PM10')==1
            emissioni=[emissioni, PM10_allp];
        elseif strcmp(precNames(i_prec,1),'PM25')==1
            emissioni=[emissioni, PM25_allp];
        elseif strcmp(precNames(i_prec,1),'SOx')==1
            emissioni=[emissioni, SO2_allp];
        elseif strcmp(precNames(i_prec,1),'NMVOC')==1
            emissioni=[emissioni, VOC_allp];
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANN INPUT DATA
%keep only optim domain
emissioni(find(flag_optim_dom==0),:)=[];
%20130820 - consider only if cell completerly in PAD
%         emissioni(find(flag_optim_dom==0 | flag_optim_dom==2),:)=[];

alpha=aggregationInfo.extraInfo.alpha;
omega=aggregationInfo.extraInfo.omega;
%change dimensions of alpha and omega to be coherent with emissions
alpha=permute(alpha,[2 1 3]);
omega=permute(omega,[2 1 3]);

dimX=size(alpha,1);
dimY=size(alpha,2);
for i=1:size(alpha,3)
    tmp=alpha(:,:,i);
    thisAlpha(:,i)=reshape(tmp,dimX*dimY,1);
end

%remove not used data
thisAlpha(commonDataInfo.flag_region_dom==0,:)=[];
%         new_alpha=reshape(alpha,dimy,dimx);
%emissioni*.aggregationInfo.firstguess.alpha
%in case it is necessary to process quadrant emissions (if too close
%to domain boundary, it is necessary to increment emissions with
%the assumptions that part of the quadrant in which emissions are
%not available, still contain same emission average)
dirs=commonDataInfo.dirs;
domainInfo=commonDataInfo.domainInfo;

%from the SR netcdf, use only NOx(1), NH3(3), PM25(5), SO2(6)
%if jj eq 0 or 1
% if ((aqiIndex == 1) || (aqiIndex == 2)) aqi_per_cell=sum(emissioni(:,[1 3 5 6]).*thisAlpha(:,[1 3 5 6]),2); end
% if (aqiIndex == 6) aqi_per_cell=sum(emissioni(:,[1]).*thisAlpha(:,[1]),2); end
if (aqiIndex == 1) aqi_per_cell=sum(emissioni(:,[1 3 4 6]).*thisAlpha(:,[1 3 4 6]),2); end %nox,nh3,pm10,so2
if (aqiIndex == 2) aqi_per_cell=sum(emissioni(:,[1 3 5 6]).*thisAlpha(:,[1 3 5 6]),2); end %nox,nh3,pm25,so2
if (aqiIndex == 6) aqi_per_cell=sum(emissioni(:,[1]).*thisAlpha(:,[1]),2); end             %nox

%change name to a better one!!! (too similar to caller...)
%save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint8_new','input_rete2');
% 20160421 MM / EP SR / First guess Version
%aqi_per_cell=interface_get_aqipercell(input_rete2, NN, aggregationInfo.mathIntermediateData, commonDataInfo, periodIndex, aqiIndex );
periodIndex=1;
%aqi_val2=firstguess_get_aqipercell(aqi_per_cell, aggregationInfo, commonDataInfo, periodIndex, aqiIndex );
aqi_per_cell=firstguess_get_aqipercell(aqi_per_cell, aggregationInfo, commonDataInfo, periodIndex, aqiIndex );

end