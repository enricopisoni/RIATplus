function [loadInfo]=load_FinalizeJobSetting(fileName)

 [loadInfo] = load_InterfaceDataSetting(fileName);
% fid=fopen(fileName);
% 
% C1=textscan(fid,'%s');
% 
% settingInfo.functionName=cell2mat(C1{1}(1));      % functionName
% settingInfo.isCoded=cell2mat(C1{1}(2)) == '1';               % function is Embedded (MatLab code)
% %aggregationInfo.linesNo=cell2mat(C1{1}(3));          % function total lines
% settingInfo.linesNo=str2num(cell2mat(C1{1}(3)));          % function total lines
% for k=1:settingInfo.linesNo 
%     codeLines(k)=C1{1}(k+3);
% end
% if (settingInfo.linesNo ~= 0) 
%     settingInfo.codeLines=codeLines;
% 
% fclose(fid);

end

