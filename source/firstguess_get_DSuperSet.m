function [ D ]= firstguess_get_DSuperSet(refData, index1)

 %[alpha, omega, radius, flatWeight, pollutantList]=FG_read(fName);
 %refData.first;
 Dtemp = refData.DSuperSet(index1);
 D= Dtemp.DSet;
 %aggregationInfo.firstguess.alpha=alpha;
 %aggregationInfo.firstguess.omega=omega;
 %aggregationInfo.firstguess.radius=radius;
 %aggregationInfo.firstguess.flatWeight=flatWeight;
 %aggregationInfo.firstguess.pollutantList=pollutantList;
 %D = refData.DOptSet(index1);

end
