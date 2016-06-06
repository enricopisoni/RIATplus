function [commonDataInfo]= quadrant_fill_commonDataInfo(commonDataInfo, periodIndex, aqiIndex)

 % leave here (dynamic setting...)
 pathANN=commonDataInfo.pathANN;
 fName=strtrim(pathANN(periodIndex).ANNs(aqiIndex,:));
 %strtrim(pathANN(k).ANNs(indaqi,:))
 % nn case
 % 20160418: MM & EP net/quadrants input
 [net]=net_read(fName);
 commonDataInfo.radius=net.icells;
 
 [commonDataInfo]=commonDataInfo;

end
