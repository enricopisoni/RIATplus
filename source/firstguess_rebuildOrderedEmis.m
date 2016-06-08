function [resEmis]=firstguess_rebuildOrderedEmis(emis, refInfo, commonDataInfo)

%ncelopt1=length(find(optimizerCondition));
%optimizerCondition=commonDataInfo.optimizerCondition;
%optimizerValues=commonDataInfo.optimizerValues;

%execString=strcat('length(find(', optimizerCondition);
%execString=strcat(execString, '))');
%ncelopt=eval(execString);

%restore first guess emissions
%MOD20160607ET
flag_region_dom=commonDataInfo.domainInfo.flag_optim_dom;
ncelopt=length(find(flag_region_dom==1 | flag_region_dom==2));
%MOD20160607ET

s1_NOX = emis((ncelopt*0)+1:ncelopt*1);

s1_VOC = emis((ncelopt*1)+1:ncelopt*2);

s1_NH3 = emis((ncelopt*2)+1:ncelopt*3);

s1_PM10 = emis((ncelopt*3)+1:ncelopt*4);

s1_PM25 = emis((ncelopt*4)+1:ncelopt*5);

s1_SO2 = emis((ncelopt*5)+1:ncelopt*6);

%create input vector with rigth order of input
resEmis.NH3_all=[s1_NH3];
resEmis.NOX_all=[s1_NOX];
resEmis.PM10_all=[s1_PM10];
resEmis.PM25_all=[s1_PM25];
resEmis.SO2_all=[s1_SO2];
resEmis.VOC_all=[s1_VOC];
%resEmis=emissioni;

end