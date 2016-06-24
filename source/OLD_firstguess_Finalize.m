function [ finalInfo ] = firstguess_Finalize( refInfo, intermediateResult, commonData, precursor, indicators, ...
    x, y, nx, ny, precName, indiciMAT, splitResult)
%    x, y, nx, ny, totalCells, inq_deb, indiciMAT, optimizerValues, optimizerCondition)

dimx=nx;
dimy=ny;
valori_emi=precursor;
rad=200;%commonData.firstguess.radius;
pollIndex=indicators;
alpha=commonData.firstguess.alpha;
omega=commonData.firstguess.omega;
%change dimensions of alpha and omega to be coherent with emissions
alpha=permute(alpha,[2 1 3]);
omega=permute(omega,[2 1 3]);

prova_emi=valori_emi(:,1);
emi=reshape(prova_emi,dimy,dimx);

[X, Y] = meshgrid(-rad(1):1:rad(1), -rad(1):1:rad(1));
%testA=unique(alpha(:,:,pollIndex(1)));
%testO=unique(omega(:,:,pollIndex(1)));
vecPrecompF=unique(omega(:,:,pollIndex(1)));
vecPrecompF(isnan(vecPrecompF))=[];
%if (isempty(vecPrecompF)==0) check=0
%else vecPrecompF(1)=[];

if isempty(vecPrecompF)==0
    for i=1:length(vecPrecompF)
        coeff=vecPrecompF(i); %could vary from 0.5 to 5_gaus and sigma_x_gaus
        sigma_y_gaus = NaN; %from 2 to rf/2 both sigma_x_gaus and sigma_x_gaus
        theta=NaN;      %from 0 to pi
        centr=rad(1)+1;
        Ftmp(:,:,i)=1./((1+sqrt(X.^2+Y.^2)).^coeff);
        vw=30;
        tmp1=Ftmp(:,:,i);
        minVal=tmp1(rad(1)+1-vw,rad(1)+1);
        tmp1(tmp1>=minVal)=NaN;
        aveVal=nanmean(nanmean(tmp1));
        
        tmp1=Ftmp(:,:,i);
        tmp1(tmp1<minVal)=aveVal;
        Ftmp(:,:,i)=tmp1;
        
        F_TBP(i)=aveVal;
        coeff_TBP(i)=coeff;        
    end
end
%%

%enlarged emission, with increased zeros all around (used to be able to
%perform simple products among emi and factor matrix
emi_enlarged=zeros(size(emi,1)+2*(rad-2),size(emi,2)+2*(rad-2));
ZeroMap=emi_enlarged;
emi_res=emi;
emi_enlarged(rad-2+1:rad-2+1+size(emi,1)-1,rad-2+1:rad-2+1+size(emi,2)-1)=emi;
[sa sb]=size(emi_enlarged);

if sum(valori_emi)~=0
    
    for indX=1+rad:sa-rad
        for indY=1+rad:sb-rad
            omegaXY=omega(indX-rad,indY-rad,pollIndex);
%             emi_res(indX-rad, indY-rad)=NaN;
            emi_res(indX-rad, indY-rad)=0;
            if (not(isnan(omegaXY)) & omegaXY>0)
                ind=find(vecPrecompF==omegaXY);
                
                mask=Ftmp(:,:,ind);
                landMask=isfinite(mask);

                mask=mask(landMask);
                gridPortion=emi_enlarged(indX-rad:indX+rad, indY-rad:indY+rad);
                gridPortion=gridPortion(landMask);
                %disp(pollIndex);
%                 disp(indX);
%                 disp(indY);
                %disp(indX-rad);
                %disp(indY-rad);
                %disp('grid');
                %disp(size(gridPortion));
                %disp(length(gridPortion));
                %disp('mask');
                %disp(size(mask));
                %disp(length(mask));
                %mask=reshape(mask,[61,61]);
                
                emi_res(indX-rad, indY-rad)=sum(sum(gridPortion.*mask));
            end
        end
    end
    
end

%save('C:\data\work\projects\riat\gauss_agg', 'gauss_agg');
finalInfo.resGrid=reshape(emi_res,nx*ny,1);

end

