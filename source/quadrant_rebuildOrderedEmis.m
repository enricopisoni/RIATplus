function [resEmis]=quadrant_rebuildOrderedEmis(emis, refInfo, commonDataInfo)

%ncelopt1=length(find(optimizerCondition));
%optimizerCondition=commonDataInfo.optimizerCondition;
%optimizerValues=commonDataInfo.optimizerValues;

%execString=strcat('length(find(', optimizerCondition);
%execString=strcat(execString, '))');
%ncelopt=eval(execString);

%restore quadrant emissions
%MOD20160607ET
flag_region_dom=commonDataInfo.domainInfo.flag_optim_dom;
ncelopt=length(find(flag_region_dom==1 | flag_region_dom==2));
%MOD20160607ET

s1_NOX = emis((ncelopt*0)+1:ncelopt*1);
s2_NOX = emis((ncelopt*1)+1:ncelopt*2);
s3_NOX = emis((ncelopt*2)+1:ncelopt*3);
s4_NOX = emis((ncelopt*3)+1:ncelopt*4);
s1_VOC = emis((ncelopt*4)+1:ncelopt*5);
s2_VOC = emis((ncelopt*5)+1:ncelopt*6);
s3_VOC = emis((ncelopt*6)+1:ncelopt*7);
s4_VOC = emis((ncelopt*7)+1:ncelopt*8);
s1_NH3 = emis((ncelopt*8)+1:ncelopt*9);
s2_NH3 = emis((ncelopt*9)+1:ncelopt*10);
s3_NH3 = emis((ncelopt*10)+1:ncelopt*11);
s4_NH3 = emis((ncelopt*11)+1:ncelopt*12);
s1_PM10 = emis((ncelopt*12)+1:ncelopt*13);
s2_PM10 = emis((ncelopt*13)+1:ncelopt*14);
s3_PM10 = emis((ncelopt*14)+1:ncelopt*15);
s4_PM10 = emis((ncelopt*15)+1:ncelopt*16);
s1_PM25 = emis((ncelopt*16)+1:ncelopt*17);
s2_PM25 = emis((ncelopt*17)+1:ncelopt*18);
s3_PM25 = emis((ncelopt*18)+1:ncelopt*19);
s4_PM25 = emis((ncelopt*19)+1:ncelopt*20);
s1_SO2 = emis((ncelopt*20)+1:ncelopt*21);
s2_SO2 = emis((ncelopt*21)+1:ncelopt*22);
s3_SO2 = emis((ncelopt*22)+1:ncelopt*23);
s4_SO2 = emis((ncelopt*23)+1:ncelopt*24);

%create input vector with rigth order of input
resEmis.NH3_all=[s1_NH3,s2_NH3,s3_NH3,s4_NH3];
resEmis.NOX_all=[s1_NOX,s2_NOX,s3_NOX,s4_NOX];
resEmis.PM10_all=[s1_PM10,s2_PM10,s3_PM10,s4_PM10];
resEmis.PM25_all=[s1_PM25,s2_PM25,s3_PM25,s4_PM25];
resEmis.SO2_all=[s1_SO2,s2_SO2,s3_SO2,s4_SO2];
resEmis.VOC_all=[s1_VOC,s2_VOC,s3_VOC,s4_VOC];
%resEmis=emissioni;

end