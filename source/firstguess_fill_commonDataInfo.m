function [commonDataInfo]= firstguess_fill_commonDataInfo(commonDataInfo, periodIndex, aqiIndex)

 % leave here (dynamic setting...)
 pathANN=commonDataInfo.pathANN;
 fName=strtrim(pathANN(periodIndex).ANNs(aqiIndex,:));
 %strtrim(pathANN(k).ANNs(indaqi,:))
 % nn case
 % 20160418: MM & EP regression input
 [alpha, omega, radius, flatWeight, pollutantList]=firstguess_read(fName);
 commonDataInfo.extraInfo.pollutantList=pollutantList;
 commonDataInfo.extraInfo.alpha=alpha;
 commonDataInfo.extraInfo.omega=omega;
 [numX, numY, numPoll]=size(commonDataInfo.extraInfo.alpha);
 
 commonDataInfo.extraInfo.radius=radius;
 commonDataInfo.extraInfo.flatWeight=flatWeight;
 
 [commonDataInfo]=commonDataInfo;

end
