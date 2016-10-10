function [aqiCLEini, aqiCLEred ]=interfaceF_getaqiCLEInfo(type, ii, jj, emiCleRed, aggregationInfo, commonDataInfo)

%if isequal(aggregationInfo.type, 'FIRSTGUESS')
if isequal(type, 'FIRSTGUESS')
    fName=strtrim(commonDataInfo.pathANN(ii).ANNs(jj,:));
    aggregationInfo.extraInfo=0;
    if (strcmp(fName,'-999') == 0)
        [alpha, omega, radius, flatWeight, pollutantList]=firstguess_read(fName);
        aggregationInfo.extraInfo.alpha=alpha;
        aggregationInfo.extraInfo.omega=omega;
        aggregationInfo.extraInfo.radius=radius;
        aggregationInfo.extraInfo.flatWeight=flatWeight;
        aggregationInfo.extraInfo.pollutantList=pollutantList;
    end
    emisP=0;
    NN=0;
    bcSet=0;
    aggregationInfo.mathInfo=0;
    [aqiCLEini]=firstguess_aggregated_scenario_mode(emiCleRed{ii,1}, aggregationInfo, commonDataInfo, ii,jj);
    [aqiCLEred]=firstguess_aggregated_scenario_mode(emiCleRed{ii,2}, aggregationInfo, commonDataInfo, ii,jj);
else
    NN=neuralNet_get_nnSuperSet_indexed(aggregationInfo.mathIntermediateData, ii, jj);
    [aqiCLEini]=MAINaggregated_scenario_mode(emiCleRed{ii,1}, NN,ii,jj, commonDataInfo, aggregationInfo);
    [aqiCLEred]=MAINaggregated_scenario_mode(emiCleRed{ii,2}, NN,ii,jj, commonDataInfo, aggregationInfo);
end

end

