function [commonDataInfo]= firstguess_fill_commonDataInfo(commonDataInfo, periodIndex, aqiIndex)

 % leave here (dynamic setting...)
 pathANN=commonDataInfo.pathANN;
 fName=strtrim(pathANN(periodIndex).ANNs(aqiIndex,:));
 %strtrim(pathANN(k).ANNs(indaqi,:))
 % nn case
 % 20160418: MM & EP regression input
 [alpha, omega, radius, flatWeight, pollutantList]=firstguess_read(fName);
 commonDataInfo.firstguess.pollutantList=pollutantList;
 commonDataInfo.firstguess.alpha=alpha;
 commonDataInfo.firstguess.omega=omega;
 [numX, numY, numPoll]=size(commonDataInfo.firstguess.alpha);
 
 commonDataInfo.firstguess.radius=radius;
 commonDataInfo.firstguess.flatWeight=flatWeight;
 
 [commonDataInfo]=commonDataInfo;

end
