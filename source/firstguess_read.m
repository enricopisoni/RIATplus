function [alpha, omega, radius, flatWeight, pollutantList]=firstguess_read(fileName)

% global attroib, check...
%  radius=30; %ncread(fileName,'Radius of influence');
% True polls list (from nc file)
%pollutantList={'NOx';'NMVOC';'NH3';'PPM';'SOx'};

flatWeight=[];%ncread(fileName,'flatWeight');

alpha=ncread(fileName,'alpha');
omega=ncread(fileName,'omega');

% workaround to put PM25 AND PM10 instead of PPM
nDim=size(alpha);
alphaNew(:,:,1:3)=alpha(:,:,1:3);
alphaNew(:,:,4)=alpha(:,:,4);
alphaNew(:,:,5)=alpha(:,:,4);
alphaNew(:,:,6)=alpha(:,:,5);

omegaNew(:,:,1:3)=omega(:,:,1:3);
omegaNew(:,:,4)=omega(:,:,4);
omegaNew(:,:,5)=omega(:,:,4);
omegaNew(:,:,6)=omega(:,:,5);

alpha=alphaNew;
omega=omegaNew;
pollutantList={'NOx';'NMVOC';'NH3';'PM10';'PM25';'SOx'};

ncid=netcdf.open(fileName,'NC_NOWRITE');
varid=netcdf.getConstant('GLOBAL');
[num_dims, num_vars, num_global_atts, unlimdimID] = netcdf.inq(ncid);
for attid = 0:num_global_atts-1
    attName = netcdf.inqAttName(ncid,varid,attid);
    if strcmp(attName,'Radius_of_influence')==1 | strcmp(attName,'Radius of influence')==1
        radius=netcdf.getAtt(ncid,varid,attName);
    end
end
%try
%    radius=netcdf.getAtt(ncid,varid,'Radius_of_influence');
%catch exception
%    radius=netcdf.getAtt(ncid,varid,'Radius of influence');
%end
%radius=netcdf.getAtt(ncid,varid,'Radius_of_influence');
tuePollutantList=netcdf.getAtt(ncid,varid,'Order_Pollutant');

%remain = pollutantList;

%polls = regexp(pollutantList,'([^ ,:]*)','tokens');

%ncread(fileName,'Order_Pollutant');%ncread(fileName,'radius');

netcdf.close(ncid);

%nc close?

end