function [ finalInfo ] = quadrant_Finalize_( refInfo, intermediateResult, commonData, precursor, indicators, ...
    x, y, nx, ny, totalCells, optimizerValues, optimizerCondition)

prova_emi=precursor(:,1);
emi=reshape(prova_emi,ny,nx);

%enlarged emission, with increased zeros all around (used to be able to
%perform simple products among emi and factor matrix
emi_enlarged=zeros(size(emi,1)+2*(intermediateResult.radius-2),size(emi,2)+2*(intermediateResult.radius-2));
emi_enlarged(intermediateResult.radius-2+1:intermediateResult.radius-2+1+size(emi,1)-1,intermediateResult.radius-2+1:intermediateResult.radius-2+1+size(emi,2)-1)=emi;
[sa sb]=size(emi_enlarged);

%
%left quadrant
%factors to perform products of cell emissions
fattore=tril(ones(intermediateResult.radius+1,intermediateResult.radius+1))-diag(repmat(0.5,intermediateResult.radius+1,1));
sotto=flipud(fattore);
fattore=[fattore;sotto(2:end,:)];
fattore(intermediateResult.radius+1,intermediateResult.radius+1)=0.25;

%initialize variable, and set dimensions
spicchio3=zeros(ny,nx);

%down quadrant
%factors to perform products of cell emissions
fattoreD=rot90(fattore,3);

%initialize variable, and set dimensions
spicchio1=zeros(ny,nx);

%right quadrant
%factors to perform products of cell emissions
fattoreR=rot90(fattore,2);

%initialize variable, and set dimensions
spicchio4=zeros(ny,nx);

%up quadrant
%factors to perform products of cell emissions
fattoreU=rot90(fattore,1);

%initialize variable, and set dimensions
spicchio2=zeros(ny,nx);

%loop to aggregate emissions
for i=intermediateResult.radius+1:sa-intermediateResult.radius
    for j=intermediateResult.radius+1:sb-intermediateResult.radius
        %left quadrant
        spicchio3(i-intermediateResult.radius+2,j-intermediateResult.radius+2)=sum(sum(emi_enlarged(i-intermediateResult.radius:i+intermediateResult.radius,j-intermediateResult.radius:j).*fattore));
        %up quadrant (down as matrix, up geographically speaking)
        spicchio1(i-intermediateResult.radius+2,j-intermediateResult.radius+2)=sum(sum(emi_enlarged(i-intermediateResult.radius:i,j-intermediateResult.radius:j+intermediateResult.radius).*fattoreD));
        %right quadrant
        spicchio4(i-intermediateResult.radius+2,j-intermediateResult.radius+2)=sum(sum(emi_enlarged(i-intermediateResult.radius:i+intermediateResult.radius,j:j+intermediateResult.radius).*fattoreR));
        %down quadrant (up as matrix, down geographically speaking)
        spicchio2(i-intermediateResult.radius+2,j-intermediateResult.radius+2)=sum(sum(emi_enlarged(i:i+intermediateResult.radius,j-intermediateResult.radius:j+intermediateResult.radius).*fattoreU));
    end
end

spi1=reshape(spicchio1,nx*ny,1);
spi2=reshape(spicchio2,nx*ny,1);
spi3=reshape(spicchio3,nx*ny,1);
spi4=reshape(spicchio4,nx*ny,1);

execString=strcat('spi1(', optimizerCondition);
execString=strcat(execString, ')');
spi1=eval(execString);

execString=strcat('spi2(', optimizerCondition);
execString=strcat(execString, ')');
spi2=eval(execString);

execString=strcat('spi3(', optimizerCondition);
execString=strcat(execString, ')');
spi3=eval(execString);

execString=strcat('spi4(', optimizerCondition);
execString=strcat(execString, ')');
spi4=eval(execString);

%resGrid=[spi1;spi2;spi3;spi4];
execString=strcat(optimizerCondition);

finalInfo.resGrid=[spi1;spi2;spi3;spi4];

end

