function [ intermediateInfo  ] = interface_prepare( refInfo, commonData)

internalRefInfo=refInfo.prepareInfo;

if (internalRefInfo.isCoded)
stringToRun='';
stringToRun=strcat(stringToRun, internalRefInfo.functionName);
stringToRun=strcat(stringToRun, '(');
stringToRun=strcat(stringToRun,'refInfo,');
stringToRun=strcat(stringToRun,'commonData');
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
[intermediateInfo]=eval(stringToRun);
else
for k=1:internalRefInfo.linesNo
    %remember to assign resGrid somewhere in the lines
    % use evalC to avoid annoying output
    T=evalc(cell2mat(internalRefInfo.codeLines(k)));
    %resGrid=resGrid(optimizerCondition);
end

end

