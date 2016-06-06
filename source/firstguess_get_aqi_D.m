function [ D ]= FG_get_aqi_D(refData, index1)

 %[alpha, omega, radius, flatWeight, pollutantList]=FG_read(fName);
 %refData.first;
 D = refData.DOptSet(index1);
 %aggregationInfo.firstguess.alpha=alpha;
 %aggregationInfo.firstguess.omega=omega;
 %aggregationInfo.firstguess.radius=radius;
 %aggregationInfo.firstguess.flatWeight=flatWeight;
 %aggregationInfo.firstguess.pollutantList=pollutantList;
 %D = refData.DOptSet(index1);

end
