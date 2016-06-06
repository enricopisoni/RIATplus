function [naggInfo]=firstguess_updateAggInfo(commonDataInfo, aggInfo, periodIndex, aqiIndex)

fName=strtrim(commonDataInfo.pathANN(periodIndex).ANNs(aqiIndex,:));
naggInfo=aggInfo
naggInfo.firstguess=0;
if (strcmp(fName,'-999') == 0)
    [alpha, omega, radius, flatWeight, pollutantList]=firstguess_read(fName);
    naggInfo.firstguess.alpha=alpha;
    naggInfo.firstguess.omega=omega;
    naggInfo.firstguess.radius=radius;
    naggInfo.firstguess.flatWeight=flatWeight;
    naggInfo.firstguess.pollutantList=pollutantList;
end

end