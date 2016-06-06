function [type,confInfo]=init_ConfigurationVar(fileName)

fid=fopen(fileName);
C1=textscan(fid,'%s');
% variables
type=cell2mat(C1{1}(1));
confInfo.commonInfoFile=cell2mat(C1{1}(2));
confInfo.geometryInfoFile=cell2mat(C1{1}(3));
confInfo.mathInfoFile=cell2mat(C1{1}(4));
confInfo.domainApplyTypeInfoFile=cell2mat(C1{1}(5));
confInfo.fifthConfInfoFile=cell2mat(C1{1}(6));
confInfo.sixthConfInfoFile=cell2mat(C1{1}(7));

fclose(fid);

end

