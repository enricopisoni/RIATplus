function [loadInfo]=loadPrepareDataSetting(fileName)

 [loadInfo] = load_InterfaceDataSetting(fileName);
% fid=fopen(fileName);
% 
% C1=textscan(fid,'%s');
% 
% prepareData.functionName=cell2mat(C1{1}(1));      % functionName
% prepareData.isCoded=cell2mat(C1{1}(2)) == '1';               % function is Embedded (MatLab code)
% prepareDataInfo.linesNo=str2num(cell2mat(C1{1}(3)));          % function total lines
% for k=1:prepareDataInfo.linesNo 
%     codeLines(k)=C1{1}(k+3);
% end
% if (prepareDataInfo.linesNo ~= 0) 
%     prepareDataInfo.codeLines=codeLines;
% 
% fclose(fid);

end

