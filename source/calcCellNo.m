function [ncel, nx, ny]=calcCellNo(type, commonDataInfo)

nx=commonDataInfo.domainInfo.special.nx;
ny=commonDataInfo.domainInfo.special.ny;
ncel=commonDataInfo.domainInfo.special.ncell; 
% if strcmp(type, 'utm')
%     % utm way to calc number of cells
%     coordinate=commonDataInfo.domainData.data(:,1:2);
%     %commonDataInfo.coordinate=coordinate;
%     %coordinate=commondDataInfo.coordinate;
%     stepsize=round(max(abs(coordinate(2,1)-coordinate(1,1)),abs(coordinate(2,2)-coordinate(1,2))));
%     nxny=round((max(coordinate)-min(coordinate))/stepsize+1);
%     nx=nxny(1,1);
%     ny=nxny(1,2);
%     ncel=nx*ny;
%  end
%  if strcmp(type, 'latlon') 
%     % latlon / ncdf 
%     ncdfFileName=commonDataInfo.pathANN(1).ANNs(1,:);
%     lat=ncread(ncdfFileName,'lat');
%     lon=ncread(ncdfFileName,'lon');
%     ny=length(lat);
%     nx=length(lon);
%     ncel=nx*ny;
%     % read info from ncdf file
%  end

end

