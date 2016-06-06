function [ finalInfo ] = quadrant_Finalize( refInfo, intermediateResult, commonData, precursor, indicators, ...
    x, y, nx, ny, precName, indiciMAT, splitResult)
%    x, y, nx, ny, totalCells, inq_deb, indiciMAT, optimizerValues, optimizerCondition)

dimx=nx;
dimy=ny;
valori_emi=precursor;
r=commonData.radius;

prova_emi=valori_emi(:,1);
emi=reshape(prova_emi,dimy,dimx);

%enlarged emission, with increased zeros all around (used to be able to
%perform simple products among emi and factor matrix
emi_enlarged=zeros(size(emi,1)+2*(r-2),size(emi,2)+2*(r-2));
ZeroMap=emi_enlarged;
emi_enlarged(r-2+1:r-2+1+size(emi,1)-1,r-2+1:r-2+1+size(emi,2)-1)=emi;
[sa sb]=size(emi_enlarged);
indiciMAT_enlarged=zeros(size(emi,1)+2*(r-2),size(emi,2)+2*(r-2));
indiciMAT_enlarged(r-2+1:r-2+1+size(indiciMAT,1)-1,r-2+1:r-2+1+size(indiciMAT,2)-1)=indiciMAT;
%
%left quadrant
% %factors to perform products of cell emissions
 fattore=tril(ones(r+1,r+1))-diag(repmat(0.5,r+1,1));
 sotto=flipud(fattore);
 fattore=[fattore;sotto(2:end,:)];
 fattore(r+1,r+1)=0.25;

%initialize variable, and set dimensions
spicchio3=zeros(dimy,dimx);

%down quadrant
 %factors to perform products of cell emissions
 fattoreD=rot90(fattore,3);

%initialize variable, and set dimensions
spicchio1=zeros(dimy,dimx);

%right quadrant
% %factors to perform products of cell emissions
 fattoreR=rot90(fattore,2);

%initialize variable, and set dimensions
spicchio4=zeros(dimy,dimx);

%up quadrant
% %factors to perform products of cell emissions
 fattoreU=rot90(fattore,1);

%initialize variable, and set dimensions
spicchio2=zeros(dimy,dimx);

% % % initialize fattore
%  fattoreS1=zeros(((sa-r)*((sa-r)-r+2-1))+((sb-r)-r+2),dimy*dimx);
%  fattoreS2=fattoreS1;
%  fattoreS3=fattoreS1;
%  fattoreS4=fattoreS1;

if sum(valori_emi)~=0
    
    %loop to aggregate emissions
    %create index matrix
    %INDx=[];
    INDx=zeros(size(r+1:sb-r,2)*((sa-r)-(r)),2);
    cic=1;
    for iv=r+1:(sa-r)
        
        i1=(size(r+1:sb-r,2)*(cic-1))+1;
        i2=((size(r+1:sb-r,2)*cic));
        INDx(i1:i2,:)=[ones(size(r+1:sb-r,2),1)*iv,(r+1:sb-r)'];
        cic=cic+1;
    end
    
    for i=1:size(INDx,1)
        if indiciMAT_enlarged(INDx(i,1),INDx(i,2))>0
            %left quadrant
            %spicchio3(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1)-r:INDx(i,1)+r,INDx(i,2)-r:INDx(i,2)))); %rectangle
            spicchio3(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1)-r:INDx(i,1)+r,INDx(i,2)-r:INDx(i,2)).*fattore)); %quadrant

            %up quadrant (down as matrix, up geographically speaking)
            %spicchio1(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1)-r:INDx(i,1),INDx(i,2)-r:INDx(i,2)+r))); %rectangle
            spicchio1(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1)-r:INDx(i,1),INDx(i,2)-r:INDx(i,2)+r).*fattoreD));

            %right quadrant
            %spicchio4(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1)-r:INDx(i,1)+r,INDx(i,2):INDx(i,2)+r)));%rectangle
            spicchio4(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1)-r:INDx(i,1)+r,INDx(i,2):INDx(i,2)+r).*fattoreR));

            %down quadrant (up as matrix, down geographically speaking)
            %spicchio2(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1):INDx(i,1)+r,INDx(i,2)-r:INDx(i,2)+r))); %rectangle
            spicchio2(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1):INDx(i,1)+r,INDx(i,2)-r:INDx(i,2)+r).*fattoreU));

            %areal sum
            %spicchio(INDx(i,1)-r+2,INDx(i,2)-r+2)=sum(sum(emi_enlarged(INDx(i,1)-r:INDx(i,1)+r,INDx(i,2)-r:INDx(i,2)+r)));
        end
    end
    
    
end

%ORIGINAL
spiw1=reshape(spicchio1,dimx*dimy,1);
spiw2=reshape(spicchio2,dimx*dimy,1);
spiw3=reshape(spicchio3,dimx*dimy,1);
spiw4=reshape(spicchio4,dimx*dimy,1);

if (splitResult == 1)
    par1_Wrong=[spiw1,spiw2,spiw3,spiw4];
else
    par1_Wrong=[spiw1;spiw2;spiw3;spiw4];
end
%save('C:\data\work\projects\riat\par1_Wrong', 'par1_Wrong');
finalInfo.resGrid=par1_Wrong;

end

