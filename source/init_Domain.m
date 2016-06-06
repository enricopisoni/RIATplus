%function [domainData, extraInfo]=init_Domain(commonDataInfo, pathFOD)
function [commonDataInfo]=init_Domain(commonDataInfo)

pathFOD=commonDataInfo.dirs.pathFOD;
% !!! include pad file info
%ntm_tm=1 (means technical measure....0 means NTM)
%NOTE: AR and RE are in %

%rows are quadruples * 2 (low and high). Quadruples are related to
%GAINS activities, without NOC

%select optimization domain: %1 means cell in optimization domain, 0 means
%outside (2 is border cell, partly inside and partly outside)
domainData=importdata(pathFOD);
if size(domainData.data,2)==6
    domainInfo.flag_region_dom=domainData.data(:,3); %regional domain
    domainInfo.flag_optim_dom=domainData.data(:,4);%optimization domain
    domainInfo.flag_aqi_dom=domainData.data(:,5); %aqi computation domain
    domainInfo.pop=domainData.data(:,6); %population
else
    domainInfo.flag_region_dom=domainData.data(:,3);
    domainInfo.flag_optim_dom=domainData.data(:,3);%optimization domain
    domainInfo.flag_aqi_dom=domainData.data(:,4); %aqi computation domain
    domainInfo.pop=domainData.data(:,5); %population
end

% load domain coordinate
coordinate=domainData.data(:,1:2);
domainData.coordinate=coordinate;
domainInfo.coordinate.xutm=coordinate(:,1);
domainInfo.coordinate.yutm=coordinate(:,2);

domainInfo.special.xutm=coordinate(:,1);
domainInfo.special.yutm=coordinate(:,2);

domainInfo.special.nx=length(unique(domainInfo.special.xutm));
domainInfo.special.ny=length(unique(domainInfo.special.yutm));

domainInfo.special.ncell=domainInfo.special.nx*domainInfo.special.ny;

%domainInfo.special.xutm=coordinate(:,1);
%domainInfo.special.yutm=coordinate(:,2);

%compute nx and ny from coordinate
% stepsize=max(coordinate(2,1)-coordinate(1,1),coordinate(2,2)-coordinate(1,2));
%EPcorrection
stepsize=round(max(abs(coordinate(2,1)-coordinate(1,1)),abs(coordinate(2,2)-coordinate(1,2))));

nxny=round((max(coordinate)-min(coordinate))/stepsize+1);
domainInfo.nx=nxny(1,1);
domainInfo.ny=nxny(1,2);
domainInfo.ncel=domainInfo.nx*domainInfo.ny;

commonDataInfo.domainData=domainData;
commonDataInfo.domainInfo=domainInfo;
%commonDataInfo.optimizerValues=domainInfo.flag_optim_dom;

[commonDataInfo] = commonDataInfo;

end