function [ refInfo  ] = interface_prepareMainData( refInfo, nx, ny, totalCells, pathANN, ...
    optimizerValues, optimizerCondition)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% final version: from external file...

if (refInfo.isCoded)
stringToRun='';
stringToRun=strcat(stringToRun, refInfo.functionName);
stringToRun=strcat(stringToRun, '(');
stringToRun=strcat(stringToRun,'refInfo,');
stringToRun=strcat(stringToRun,'nx,');
stringToRun=strcat(stringToRun,'ny,');
stringToRun=strcat(stringToRun,'totalCells,');
stringToRun=strcat(stringToRun,'pathANN,');
stringToRun=strcat(stringToRun,'optimizerValues,');
stringToRun=strcat(stringToRun,'optimizerCondition');
stringToRun=strcat(stringToRun,')');
[refInfo]=eval(stringToRun);
else
for k=1:refInfo.linesNo
    %remember to assign resGrid somewhere in the lines
    % use evalC to avoid annoying output
    T=evalc(cell2mat(refInfo.codeLines(k)));
    %resGrid=resGrid(optimizerCondition);
end

end

