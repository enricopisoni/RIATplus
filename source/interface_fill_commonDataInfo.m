function [ setInfo  ] = interface_fill_commonDataInfo( commonData, aqiIndex, periodIndex)

stringToRun='';
stringToRun=strcat(stringToRun, refInfo.fill_CommonDataInfoFunction);
stringToRun=strcat(stringToRun, '(');
%stringToRun=strcat(stringToRun,'refInfo,');
stringToRun=strcat(stringToRun,'commonData');
stringToRun=strcat(stringToRun,'aqiIndex,');
stringToRun=strcat(stringToRun,'periodIndex');
%stringToRun=strcat(stringToRun,'commonData,');
%stringToRun=strcat(stringToRun,'precursor,');
%stringToRun=strcat(stringToRun,'indicators,');
%stringToRun=strcat(stringToRun,'x,');
%stringToRun=strcat(stringToRun,'y,');
%stringToRun=strcat(stringToRun,'nx,');
%stringToRun=strcat(stringToRun,'ny,');
%stringToRun=strcat(stringToRun,'totalCells,');
%stringToRun=strcat(stringToRun,'optimizerValues,');
%stringToRun=strcat(stringToRun,'optimizerCondition');
stringToRun=strcat(stringToRun,')');
[setInfo]=eval(stringToRun);