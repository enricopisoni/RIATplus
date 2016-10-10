function [naggInfo]=firstguess_updateAggInfo(commonDataInfo, aggInfo, periodIndex, aqiIndex)

fName=strtrim(commonDataInfo.pathANN(periodIndex).ANNs(aqiIndex,:));
naggInfo=aggInfo
naggInfo.extraInfo=0;
if (strcmp(fName,'-999') == 0)
    [alpha, omega, radius, flatWeight, pollutantList]=firstguess_read(fName);
    naggInfo.extraInfo.alpha=alpha;
    naggInfo.extraInfo.omega=omega;
    naggInfo.extraInfo.radius=radius;
    naggInfo.extraInfo.flatWeight=flatWeight;
    naggInfo.extraInfo.pollutantList=pollutantList;
end

end