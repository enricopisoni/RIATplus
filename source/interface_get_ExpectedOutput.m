function [ finalMatrix ] = interface_get_ExpectedOutput(refInfo, commonData, nx, ny, nPolls)

stringToRun='';
stringToRun=strcat(stringToRun, refInfo.get_expectedOutputFunction);
stringToRun=strcat(stringToRun, '(');
stringToRun=strcat(stringToRun,'refInfo,');
stringToRun=strcat(stringToRun,'commonData,');
stringToRun=strcat(stringToRun,'nx,');
stringToRun=strcat(stringToRun,'ny,');
stringToRun=strcat(stringToRun,'nPolls');
stringToRun=strcat(stringToRun,')');
[finalMatrix]=eval(stringToRun);

end