function [cycleNN, cycleD, cyclebc, extraInfo]=interfaceF_getcompute_aqi_Info(type, aggInfo, pathANN, ii, jj)

extraInfo=[];
if isequal(type, 'FIRSTGUESS')
    aggInfo.mathIntermediateData=0;
    fName=strtrim(pathANN(ii).ANNs(jj,:));
    aggInfo.extraInfo=0;
    if (strcmp(fName,'-999') == 0)
        [alpha, omega, radius, flatWeight, pollutantList]=firstguess_read(fName);
        extraInfo.alpha=alpha;
        extraInfo.omega=omega;
        extraInfo.radius=radius;
        extraInfo.flatWeight=flatWeight;
        extraInfo.pollutantList=pollutantList;
    end
    cycleNN=0;
    cyclebc=0;
    cycleD=quadrant_get_DSuperSet_indexed(aggInfo.geometryIntermediateData, ii, jj);
else
    fName=strtrim(pathANN(ii).ANNs(jj,:));
    [net]=net_read(fName);
    commonDataInfo.radius=net.icells;
    cycleD=quadrant_get_DSuperSet_indexed(aggInfo.geometryIntermediateData, ii, jj);
    cycleNN=neuralNet_get_nnSuperSet_indexed(aggInfo.mathIntermediateData, ii, jj);
    cyclebc=quadrant_get_bcSuperSet_indexed(aggInfo.mathIntermediateData, ii, jj);
end

end