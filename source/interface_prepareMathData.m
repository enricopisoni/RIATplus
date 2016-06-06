function [ mathInfo  ] = interface_prepareMathData( precursors, indicators, x, y, nx, ny, totalCells, prepareDataInfo, ...
    optimizerValues, optimizerCondition)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% final version: from external file...

if (prepareDataInfo.isCoded)
stringToRun='';
stringToRun=strcat(stringToRun, prepareDataInfo.functionName);
stringToRun=strcat(stringToRun, '(');
stringToRun=strcat(stringToRun,'precursors,');
stringToRun=strcat(stringToRun,'indicators,');
stringToRun=strcat(stringToRun,'x,');
stringToRun=strcat(stringToRun,'y,');
stringToRun=strcat(stringToRun,'nx,');
stringToRun=strcat(stringToRun,'ny,');
stringToRun=strcat(stringToRun,'totalCells,');
stringToRun=strcat(stringToRun,'prepareDataInfo,');
stringToRun=strcat(stringToRun,'optimizerValues,');
stringToRun=strcat(stringToRun,'optimizerCondition');
stringToRun=strcat(stringToRun,')');
[mathInfo]=eval(stringToRun);
else
for k=1:prepareDataInfo.linesNo
    %remember to assign resGrid somewhere in the lines
    % use evalC to avoid annoying output
    T=evalc(cell2mat(prepareDataInfo.codeLines(k)));
    %resGrid=resGrid(optimizerCondition);
end

end

