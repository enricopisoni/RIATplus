function [extraInfo]=interfaceF_fillAggregationInfoFromFile(type, pathAnn, index1, index2)

extraInfo=[];
if isequal(type, 'FIRSTGUESS')
    %if (isequal(strtrim(aggregationInfo.type),'FIRSTGUESS')==1)
    %fName=strtrim(commonDataInfo.pathANN(1).ANNs(commonDataInfo.aqi_obj+1,:));
    fName=strtrim(pathAnn(index1).ANNs(index2,:));
    if (strcmp(fName,'-999') == 0)
        [alpha, omega, radius, flatWeight, pollutantList]=firstguess_read(fName);
        extraInfo.alpha=alpha;
        extraInfo.omega=omega;
        extraInfo.radius=radius;
        extraInfo.flatWeight=flatWeight;
        extraInfo.pollutantList=pollutantList;
    end
end

end