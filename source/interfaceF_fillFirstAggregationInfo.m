function [extraInfo]=interfaceF_fillFirstAggregationInfo(type, commonDataInfo)

%if isequal(aggregationInfo.type, 'FIRSTGUESS')
if isequal(type, 'FIRSTGUESS')
    %         pollN=size(commonDataInfo.pathANN(1).ANNs(:,:), 1);
    %         for i = 1:pollN
    %             fName=strtrim(commonDataInfo.pathANN(1).ANNs(i,:));
    %             aggregationInfo.firstguess=0;
    %             if (strcmp(fName,'-999') == 0)
    %                 [alpha, omega, radius, flatWeight, pollutantList]=firstguess_read(fName);
    %                 aggregationInfo.firstguess.alpha=alpha;
    %                 aggregationInfo.firstguess.omega=omega;
    %                 aggregationInfo.firstguess.radius=radius;
    %                 aggregationInfo.firstguess.flatWeight=flatWeight;
    %                 aggregationInfo.firstguess.pollutantList=pollutantList;
    %                 break;
    %             end
    %         end
    [alpha, omega, radius, flatWeight, pollutantList]=firstguess_getFirstValidInfo(commonDataInfo);
    extraInfo.alpha=alpha;
    extraInfo.omega=omega;
    extraInfo.radius=radius;
    extraInfo.flatWeight=flatWeight;
    extraInfo.pollutantList=pollutantList;
    %A=aggregationInfo.fullA;
    %B=aggregationInfo.fullB;
    %LB=aggregationInfo.fullLB;
    %UB=aggregationInfo.fullUB;
else
    extraInfo=0;
end

end

