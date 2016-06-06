function [ commonDataInfo, aggregationInfo ] = firstguess_fillcomputesol(commonDataInfo,refInfo, aggregationInfo,periodIndex,aqiIndex)

fName=strtrim(commonDataInfo.pathANN(periodIndex).ANNs(aqiIndex,:));
aggregationInfo.firstguess=0;
if (strcmp(fName,'-999') == 0)
    [alpha, omega, radius, flatWeight, pollutantList]=firstguess_read(fName);
    aggregationInfo.firstguess.alpha=alpha;
    aggregationInfo.firstguess.omega=omega;
    aggregationInfo.firstguess.radius=radius;
    aggregationInfo.firstguess.flatWeight=flatWeight;
    aggregationInfo.firstguess.pollutantList=pollutantList;
end

end
