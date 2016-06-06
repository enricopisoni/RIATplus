function [loadInfo]=load_AggregationDataSetting(fileName)

 fid=fopen(fileName);
 C1=textscan(fid,'%s');
 % variables
 prepareFile=cell2mat(C1{1}(1));          % prepare file
 finalizeFile=cell2mat(C1{1}(2));          % finalize file
 restoreFile=cell2mat(C1{1}(3));          % restore file
 
 [prepareInfo] = load_InterfaceDataSetting(prepareFile);
 [prepareInfo] = load_InterfaceDataSetting(prepareFile);
 [finalizeInfo] = load_InterfaceDataSetting(finalizeFile);
 [restoreInfo] = load_InterfaceDataSetting(restoreFile);
 
 loadInfo.prepareInfo=prepareInfo;
 loadInfo.finalizeInfo=finalizeInfo;
 loadInfo.restoreInfo=restoreInfo;

end

