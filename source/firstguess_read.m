function [alpha, omega, radius, flatWeight, pollutantList]=firstguess_read(fileName)

 % global attroib, check...
%  radius=30; %ncread(fileName,'Radius of influence');
 %pollutantList=ncread(fileName,'Order_Pollutant');%ncread(fileName,'radius');
 %pollutantList={'NOx';'NMVOC';'NH3';'PM25';'SOx'};
 pollutantList={'NOx';'NMVOC';'NH3';'PM10';'PM25';'SOx'};
  
 flatWeight=[];%ncread(fileName,'flatWeight');

 alpha=ncread(fileName,'alpha');
 omega=ncread(fileName,'omega');
 
 ncid=netcdf.open(fileName,'NC_NOWRITE');
 varid=netcdf.getConstant('GLOBAL');
 radius=netcdf.getAtt(ncid,varid,'Radius of influence');
 netcdf.close(ncid);
  
 %nc close?

end