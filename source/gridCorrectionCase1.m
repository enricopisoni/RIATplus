function [domainData, domainInfo, emi]=gridCorrectionCase1(origDomainData, origDomainInfo)

 domainData=origDomainData;
 domainInfo=origDomainInfo;
 
 domainData.ny=domainData.ny-3; %problem in RIAT about number of cells
 domainData.ncel=domainData.nx*domainData.ny;
 [dataerase]=find(domainData.coordinate(:,1)>874.28 | domainData.coordinate(:,2)>5210.7);
 domainInfo.flag_optim_dom=domainInfo.flag_optim_dom(dataerase,:)=[];
 domainData.coordinate(dataerase,:)=[];
 domainInfo.pop(dataerase,:)=[];
 domainInfo.flag_aqi_dom(dataerase,:)=[];
 %redefine xutm and yutm
 domainInfo.special.xutm=domainData.coordinate(:,1);
 domainInfo.special.yutm=domainData.coordinate(:,2);
 %process files to get the following format:
 % xutm1 yutm1
 % xutm1 yutm2
 % ...
 %while now format is
 % xutm1 yutm1
 % xutm2 yutm1
 % ...
 %new regional flag
 flag_optim_dom1=reshape(domainInfo.flag_optim_dom,domainInfo.nx,domainInfo.ny)';
 flag_optim_dom2=reshape(flag_optim_dom1,domainInfo.ncel,1);
 domainInfo.flag_optim_dom=flag_optim_dom2;
        
 %new pop file
 pop1=reshape(domainInfo.pop,nx,ny)';
 pop2=reshape(pop1,nx*ny,1);
 domainInfo.pop=pop2;
        
 %new aqi compuation domain file
 flag_aqi_dom1=reshape(domainInfo.flag_aqi_dom,domainInfo.nx,domainInfo.ny)';
 flag_aqi_dom2=reshape(flag_aqi_dom1,domainInfo.ncel,1);
 domainInfo.flag_aqi_dom=flag_aqi_dom2;
        
 %new xutm
 xutm1=reshape(domainInfo.xutm,domainInfo.nx,domainInfo.ny)';
 xutm2=reshape(xutm1,domainInfo.ncel,1);
        
 %new xutm
 yutm1=reshape(domainInfo.yutm,nx,ny)';
 yutm2=reshape(yutm1,nx*ny,1);
        
 %new coordinates
 domainData.coordinate=[xutm2 yutm2];
 domainInfo.xutm=xutm2;
 domainInfo.yutm=yutm2;
 %clipping data
 emi(dataerase,:)=[];
 %taking into account both aggregated and detailed emissions - this part
 %has been added to manage "boundary cells" (cells partly inside and
 %partly outside the optimization domain)
 emiAGG=emi(:,1);
 emiDET=emi(:,2);
 emiAGG1=reshape(emiAGG,nx,ny)';
 emiAGG2=reshape(emiAGG1,nx*ny,1);
 emiDET1=reshape(emiDET,nx,ny)';
 emiDET2=reshape(emiDET1,nx*ny,1);
 emi=[];
 emi=[emiAGG2 emiDET2];

end