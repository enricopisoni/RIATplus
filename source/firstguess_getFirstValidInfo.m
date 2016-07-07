function [alpha, omega, radius, flatWeight, pollutantList]=firstguess_getFirstValidInfo(commonDataInfo)

pollN=size(commonDataInfo.pathANN(1).ANNs(:,:), 1);
for i = 1:pollN
%     fName=strtrim(commonDataInfo.pathANN(1).ANNs(i,:));
    %READ CORRECT FILE FOR OPTIMIZING
    fName=strtrim(commonDataInfo.pathANN(1).ANNs(commonDataInfo.aqi_obj+1,:))
    aggregationInfo.firstguess=0;
    if (strcmp(fName,'-999') == 0)
        [alpha1, omega1, radius1, flatWeight1, pollutantList1]=firstguess_read(fName);
        alpha=alpha1;
        omega=omega1;
        radius=radius1;
        flatWeight=flatWeight1;
        pollutantList=pollutantList1;
        break;
    end
end

end