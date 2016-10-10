function [ commonDataInfo, aggregationInfo ] = firstguess_fillcomputesol(commonDataInfo,refInfo, aggregationInfo,periodIndex,aqiIndex)

fName=strtrim(commonDataInfo.pathANN(periodIndex).ANNs(aqiIndex,:));
aggregationInfo.extraInfo=0;
if (strcmp(fName,'-999') == 0)
    [alpha, omega, radius, flatWeight, pollutantList]=firstguess_read(fName);
    aggregationInfo.extraInfo.alpha=alpha;
    aggregationInfo.extraInfo.omega=omega;
    aggregationInfo.extraInfo.radius=radius;
    aggregationInfo.extraInfo.flatWeight=flatWeight;
    aggregationInfo.extraInfo.pollutantList=pollutantList;
end

end
