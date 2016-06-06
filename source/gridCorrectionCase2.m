function [domainData, domainInfo, emi]=gridCorrectionCase2(origDomainData, origDomainInfo)

 domainData=origDomainData;
 domainInfo=origDomainInfo;
 
 flag_optim_dom1=reshape(domainInfo.flag_optim_dom,domainInfo.ny,domainInfo.nx);
 flag_optim_dom2=reshape(flipud(flag_optim_dom1),domainInfo.ncel,1);
 domainInfo.flag_optim_dom=flag_optim_dom2;
    
 %new pop file
 pop1=reshape(domainInfo.pop,domainInfo.ny,domainInfo.nx);
 pop2=reshape(flipud(pop1),domainInfo.ncel,1);
 domainInfo.pop=pop2;
 
 %new aqi compuation domain file
 flag_aqi_dom1=reshape(domainInfo.flag_aqi_dom,domainInfo.ny,domainInfo.nx);
 flag_aqi_dom2=reshape(flipud(flag_aqi_dom1),domainInfo.ncel,1);
 domainInfo.flag_aqi_dom=flag_aqi_dom2;
 
 %new xutm
 xutm1=reshape(xutm,domainInfo.ny,domainInfo.nx);
 xutm2=reshape(flipud(xutm1),domainInfo.ncel,1);
 
 %new xutm
 yutm1=reshape(yutm,domainInfo.ny,domainInfo.nx);
 yutm2=reshape(flipud(yutm1),domainInfo.ncel,1);
 
 %new coordinates
 coordinate=[xutm2 yutm2];
 domainInfo.xutm=xutm2;
 domainInfo.yutm=yutm2;
 domainInfo.special.xutm=xutm2;
 domainInfo.special.yutm=yutm2;
 
 emiAGG=emi(:,1);
 emiDET=emi(:,2);
 emiAGG1=reshape(emiAGG,domainInfo.ny,domainInfo.nx);
 emiAGG2=reshape(flipud(emiAGG1),domainInfo.ncel,1);
 emiDET1=reshape(emiDET,domainInfo.ny,domainInfo.nx);
 emiDET2=reshape(flipud(emiDET1),domainInfo.ncel,1);
 emi=[];
 emi=[emiAGG2 emiDET2];

end