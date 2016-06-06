function [interfaceDataInfo]=load_InterfaceDataSetting(fileName)

fid=fopen(fileName);

C1=textscan(fid,'%s');

interfaceDataInfo.functionName=cell2mat(C1{1}(1));      % functionName
interfaceDataInfo.isCoded=cell2mat(C1{1}(2)) == '1';    % function is Embedded (MatLab code)
interfaceDataInfo.linesNo=0;
if (interfaceDataInfo.isCoded ~= 1)
    interfaceDataInfo.linesNo=str2num(cell2mat(C1{1}(3))); % function total lines
    for k=1:interfaceDataInfo.linesNo 
        codeLines(k)=C1{1}(k+3);
    end
    if (interfaceDataInfo.linesNo ~= 0) 
        interfaceDataInfo.codeLines=codeLines;
    end
end

fclose(fid);

end

