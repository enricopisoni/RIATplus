function [pollList]= firstguess_getOrderedPoll(commonDataInfo, aggregationInfo, periodIndex, aqiIndex)

 pathANN=commonDataInfo.pathANN;
 fName=strtrim(pathANN(periodIndex).ANNs(aqiIndex,:));
 %strtrim(pathANN(k).ANNs(indaqi,:))
 % nn case
 % 20160418: MM & EP net/quadrants input
 %[net]=net_read(fName);
 % commonDataInfo.radius=net.icells;
 % 20160418: MM & EP regression input
 pollList=commonDataInfo.firstguess.pollutantList;
 %[alpha, omega, radius, flatWeight, pollutantList]=firstguess_read(fName);
 %commonDataInfo.firstguess.pollutantList=pollutantList;
 %commonDataInfo.firstguess.alpha=alpha;
 %commonDataInfo.firstguess.omega=omega;
 %[numX, numY, numPoll]=size(commonDataInfo.firstguess.alpha);
 
 % 5 to 6 polls...
%  newAlpha=zeros(numX, numY, numPoll+1);
%  newAlpha(:,:,1:5)=commonDataInfo.firstguess.alpha;
%  newAlpha(:,:, 6)=commonDataInfo.firstguess.alpha(:,:,5);
%  commonDataInfo.firstguess.alpha=newAlpha;

%  newOmega=zeros(numX, numY, numPoll+1);
%  newOmega(:,:,1:5)=commonDataInfo.firstguess.omega;
%  newOmega(:,:, 6)=commonDataInfo.firstguess.omega(:,:,5);
%  commonDataInfo.firstguess.omega=newOmega;
 
 %commonDataInfo.firstguess.radius=radius;
 %commonDataInfo.firstguess.flatWeight=flatWeight;
 % end MM
 
 [commonDataInfo]=commonDataInfo;

end
