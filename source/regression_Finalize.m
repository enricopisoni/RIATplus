function [ finalInfo  ] = regression_Finalize( finalizeJobInfo, prepareDataInfo, intermediateResult, precursors, indicators, ...
    x, y, nx, ny, totalCells, optimizerValues, optimizerCondition)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% final version: from external file...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIAL SETTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%moved from prepare section...

%restore regression info
%nameRegFile=strcat(intermediateResult.nameDirOut,'regression.mat');

aa=size(precursors);
precs=precursors;
if (length(aa) == 2)
  if (aa(1) == 1) || (aa(2) == 1) precs=reshape(precursors,ny(1),nx(1));
  end
end

Prec2=zeros(ny+intermediateResult.rad(1)*2,nx+intermediateResult.rad(1)*2);
%transpose on rows/column...
Prec2(intermediateResult.rad+1:end-intermediateResult.rad,intermediateResult.rad+1:end-intermediateResult.rad)=precs(:,:);
PrecOld=precs;
precursors=[];
precursors=Prec2;

%load Indicator
% [Indic,IndicBC]=ReadIndicIneris7(nSc,nPrec,domain,aqiFil,aqi,absDel); %0=abs, 1=delta
% load(strcat('./input/',domain,'/flag_7km_ms.mat'));
%load(strcat('./input/','ineris_7km_deliver_20150831','/flagRegioMat.mat'));
% load poland_flag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IDENTIFICATION BEGIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s%%%%%%%%%%%%%%%%%%%%%%%%%%
    %alpha=NaN(ny,nx,intermediateResult.nPrec);
    alpha=NaN(ny,nx);
    
    %create empty variable only if omega does not exist...otherwise it
    %means that an optimized omega variable has been already loaded
    if exist('omega','var')==0
        %omega=NaN(ny,nx,intermediateResult.nPrec);
        omega=NaN(ny,nx);
    end
    
    if exist('sigX','var')==0
        %sigX=NaN(ny,nx,intermediateResult.nPrec);
        sigX=NaN(ny,nx);
    end
    
    if exist('sigY','var')==0
        %sigY=NaN(ny,nx,intermediateResult.nPrec);
        sigY=NaN(ny,nx);
    end
    
    if exist('thet','var')==0
        %thet=NaN(ny,nx,intermediateResult.nPrec);
        thet=NaN(ny,nx);
    end
    
    %XMin=NaN(ny,nx,intermediateResult.nPrec);
    %XMax=NaN(ny,nx,intermediateResult.nPrec);
    XMin=NaN(ny,nx);
    XMax=NaN(ny,nx);
    yMin=NaN(ny,nx);
    yMax=NaN(ny,nx);
    %bInt=NaN(ny,nx,intermediateResult.nPrec,2);
    bInt=NaN(ny,nx,2);
    
    for ic=1:nx
        disp(strcat('Creating regression on x: _',int2str(ic),' of _',int2str(nx)));
        for ir=1:ny
            %if flagRegioMat(ir,ic)==1
                
                
                %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                 %%%%%%% OMEGA CHANCING PER CELL, OPTIMIZED IN THIS CODE%%%%
                %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                 %                 x0=[0 0 0 0 0 1.5];
                %                 %                 lb=[-1 -1 -1 -1 -1 1];
                %                 %                 ub=[1 1 1 1 1 2];
                %                 opts = statset('TolFun',10^-6,'Display','off');
                %                 PrecDummyQuad=squeeze(Prec(ir:ir+rad+rad,ic:ic+rad+rad,Ide,:));
                %                 [IndicEq]=EquaIndic(ic,ir,rf,nx,ny,size(Ide,2),Indic(:,:,Ide)); %indicator
                %                 %                 [mdl,rrorfitted]=lsqcurvefit(@(beta,PrecDummyQuad) InvDistN_opt_NoNorm(beta,PrecDummyQuad,rad(1),size(Ide,2),nPrec),x0,PrecDummyQuad,double(IndicEq),lb,ub,opts);                %                 toc
                %                 mdl=nlinfit(PrecDummyQuad,double(IndicEq),@(beta,PrecDummyQuad) InvDistN_opt_NoNorm(beta,PrecDummyQuad,rad(1),size(Ide,2),nPrec),x0,opts);
                %                 omega(ir,ic,:)=mdl(6);
                %                 alpha(ir,ic,:)=mdl(1:5);
                %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                %
                
                
                
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%% OMEGA CHANCING PER CELL AND/OR PRECURSOR,
%                 %%%%%% OPTIMIZED OUTSIDE THIS CODE %%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 PrecDummyQuad=squeeze(Prec(ir:ir+rad+rad,ic:ic+rad+rad,Ide,:));
%                 
%                 switch function_type
%                     case 1
%                         coeff=squeeze(omega(ir,ic,:)); %could vary from 0.5 to 5
%                         centr=rad(1)+1;
%                         for poll=1:length(coeff)
%                             F(:,:,1,poll)=1./(1+sqrt(X.^2+Y.^2).^coeff(poll));
%                             %                                             F(centr,centr,1,poll)=1;
%                         end
%                         
%                     case 2
%                         A = 1; x0 = 0; y0 = 0;
%                         coeff=NaN;
%                         sigma_x_gaus = squeeze(sigX(ir,ic,:));
%                         sigma_y_gaus = squeeze(sigY(ir,ic,:));
%                         theta=squeeze(thet(ir,ic,:));
%                         for poll=1:nPrec
%                             a = cos(theta(poll))^2/2/sigma_x_gaus(poll)^2 + sin(theta(poll))^2/2/sigma_y_gaus(poll)^2;
%                             b = -sin(2*theta(poll))/4/sigma_x_gaus(poll)^2 + sin(2*theta(poll))/4/sigma_y_gaus(poll)^2 ;
%                             c = sin(theta(poll))^2/2/sigma_x_gaus(poll)^2 + cos(theta(poll))^2/2/sigma_y_gaus(poll)^2;
%                             F(:,:,1,poll) = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;
%                         end
%                 end
%                 Fproc=(repmat(F,1,1,length(Ide),1));
%                 %                 Fproc=(repmat(F,1,1,length(Ide),nPrec));
%                 PrecPatch=squeeze(sum(sum(PrecDummyQuad.*Fproc,2),1));
%                 %                                                 [PrecPatch]=EquaPrec(ic,ir,rf,nx,ny,nSc,nPrec,Prec(:,:,Ide,:),rad,Fproc); %FOR INCREASED RF
%                 [IndicEq]=EquaIndic(ic,ir,rf,nx,ny,size(Ide,2),Indic(:,:,Ide)); %indicator
%                 [modelCoef, xmin, xmax, ymin, ymax,bint]=Regression_MinMax(PrecPatch,IndicEq,pcaFlag); %YES PCA YES NORM
%                 alpha(ir,ic,:)=modelCoef;
%                 sigX(ir,ic,:)=sigma_x_gaus;
%                 sigY(ir,ic,:)=sigma_x_gaus;
%                 thet(ir,ic,:)=theta;
%                 XMin(ir,ic,:)=xmin';
%                 XMax(ir,ic,:)=xmax';
%                 yMin(ir,ic)=ymin;
%                 yMax(ir,ic)=ymax;
%                 bInt(ir,ic,:,:)=bint;
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 
                
                
                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%% CASE FOR OMEGA FIXED%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %PrecDummyQuad=squeeze(precursors(ir:ir+intermediateResult.rad+intermediateResult.rad,ic:ic+intermediateResult.rad+intermediateResult.rad,intermediateResult.Ide,:));
                                PrecDummyQuad=squeeze(precursors(ir:ir+intermediateResult.rad+intermediateResult.rad,ic:ic+intermediateResult.rad+intermediateResult.rad));
                                %Fproc=(repmat(intermediateResult.F,1,1,length(intermediateResult.Ide),intermediateResult.nPrec));
                                Fproc=(repmat(intermediateResult.F,1));
                                PrecPatch=squeeze(sum(sum(PrecDummyQuad.*Fproc,2),1));
                                %                                 [PrecPatch]=EquaPrec(ic,ir,rf,nx,ny,nSc,nPrec,Prec(:,:,Ide,:),rad,Fproc); %FOR INCREASED RF
                                %[IndicEq]=EquaIndic(ic,ir,intermediateResult.rf,nx,ny,size(intermediateResult.Ide,2),indicators(:,:,intermediateResult.Ide)); %indicator
                                [IndicEq]=EquaIndic(ic,ir,intermediateResult.rf,nx,ny,indicators(:,:)); %indicator
                                [modelCoef, xmin, xmax, ymin, ymax,bint]=Regression_MinMax(PrecPatch,IndicEq,intermediateResult.pcaFlag); %YES PCA YES NORM
                                alpha(ir,ic,:)=modelCoef;
                                omega(ir,ic,:)=intermediateResult.coeff;
                                sigX(ir,ic,:)=intermediateResult.sigma_x_gaus;
                                sigY(ir,ic,:)=intermediateResult.sigma_y_gaus;
                                thet(ir,ic,:)=intermediateResult.theta;
                                XMin(ir,ic,:)=xmin';
                                XMax(ir,ic,:)=xmax';
                                yMin(ir,ic)=ymin;
                                yMax(ir,ic)=ymax;
                                bInt(ir,ic,:,:)=bint;
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    %save(nameRegFile,'alpha','omega','sigX','sigY','thet','XMin','XMax','yMin','yMax','bInt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE TO NECTDF BEGIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for k=1:5
%     alphaF(:,:,k)=fliplr(alpha(:,:,k)');
%     
%     omegaF(:,:,k)=fliplr(omega(:,:,k)');
%     
%     sigXF(:,:,k)=fliplr(sigX(:,:,k)');
%     sigYF(:,:,k)=fliplr(sigY(:,:,k)');
%     thetF(:,:,k)=fliplr(thet(:,:,k)');
%     
%     x_min(:,:,k)=fliplr(XMin(:,:,k)');
%     x_max(:,:,k)=fliplr(XMax(:,:,k)');
% end
% y_min=fliplr(yMin');
% y_max=fliplr(yMax');
% 
% latitude=fliplr(y');
% longitude=fliplr(x');
% 
% mode = netcdf.getConstant('NETCDF4');
% ncid = netcdf.create(strcat(intermediateResult.nameDirOut,'SR_PM25_annual_20150909_rf',int2str(intermediateResult.rad(1)),'.nc'),mode);
% 
% latDimId = netcdf.defDim(ncid,'latitude',448);
% lonDimId = netcdf.defDim(ncid,'longitude',384);
% pollDimId = netcdf.defDim(ncid,'pollutant',5);
% 
% varid_lat = netcdf.defVar(ncid,'lat','double',[lonDimId latDimId]);
% varid_lon = netcdf.defVar(ncid,'lon','double',[lonDimId latDimId]);
% varid_alpha = netcdf.defVar(ncid,'alpha','double',[lonDimId latDimId pollDimId]);
% varid_omega = netcdf.defVar(ncid,'omega','double',[lonDimId latDimId pollDimId]);
% 
% varid_sigX = netcdf.defVar(ncid,'sigX','double',[lonDimId latDimId pollDimId]);
% varid_sigY = netcdf.defVar(ncid,'sigY','double',[lonDimId latDimId pollDimId]);
% varid_thet = netcdf.defVar(ncid,'thet','double',[lonDimId latDimId pollDimId]);
% 
% varid_x_min = netcdf.defVar(ncid,'x_min','double',[lonDimId latDimId pollDimId]);
% varid_x_max = netcdf.defVar(ncid,'x_max','double',[lonDimId latDimId pollDimId]);
% varid_y_min = netcdf.defVar(ncid,'y_min','double',[lonDimId latDimId]);
% varid_y_max = netcdf.defVar(ncid,'y_max','double',[lonDimId latDimId]);
% varid = netcdf.getConstant('GLOBAL');
% 
% netcdf.putVar(ncid,varid_lat,latitude);
% netcdf.putVar(ncid,varid_lon,longitude);
% netcdf.putVar(ncid,varid_alpha,alphaF);
% netcdf.putVar(ncid,varid_omega,omegaF);
% 
% netcdf.putVar(ncid,varid_sigX,sigXF);
% netcdf.putVar(ncid,varid_sigY,sigYF);
% netcdf.putVar(ncid,varid_thet,thetF);
% 
% netcdf.putVar(ncid,varid_x_min,x_min);
% netcdf.putVar(ncid,varid_x_max,x_max);
% netcdf.putVar(ncid,varid_y_min,y_min);
% netcdf.putVar(ncid,varid_y_max,y_max);
% netcdf.putAtt(ncid,varid,'Order_Pollutant','NOx, NMVOC, NH3, PM25, SOx');
% netcdf.putAtt(ncid,varid,'Equation','Exponential emission aggregation')
% netcdf.putAtt(ncid,varid,'Equation information','Note that emissions and concentrations have to be normalized between min and max values')
% netcdf.putAtt(ncid,varid,'Radius of influence',intermediateResult.rad(1))
% 
% netcdf.putAtt(ncid,varid,'Normalization is used? 0=noNorm, 1=yesNorm',intermediateResult.pcaFlag)
% 
% netcdf.close(ncid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VALIDATION BEGIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load(strcat('./input/',intermediateResult.domain,'/flagRegioMat.mat'));

    finalResult=nan(ny,nx);
    for ic=1:nx
        for ir=1:ny %use same routine for training and validation, with 0 for no sliding F
            %             [ir ic]
                
%                 %%%%%%%%%%%%%%%%%% OPTIMZED OMEGA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 [X, Y] = meshgrid(-rad:1:rad, -rad:1:rad);
%                 
%                 switch function_type
%                     case 1
%                         F=[];
%                         coeff=squeeze(omega(ir,ic,:));  %could vary from 0.1 to 3
% %                         centr=rad(1)+1;
%                         for poll=1:length(coeff)
%                             F(:,:,poll)=1./(1+sqrt(X.^2+Y.^2).^coeff(poll));
%                             %                                             F(centr,centr,poll)=1;
%                         end
%                     case 2
%                         A = 1; x0 = 0; y0 = 0;
%                         sigma_x_gaus = squeeze(sigX(ir,ic,:));
%                         sigma_y_gaus = squeeze(sigY(ir,ic,:));
%                         theta=squeeze(thet(ir,ic,:));
%                         for poll=1:nPrec
%                             a = cos(theta(poll))^2/2/sigma_x_gaus(poll)^2 + sin(theta(poll))^2/2/sigma_y_gaus(poll)^2;
%                             b = -sin(2*theta(poll))/4/sigma_x_gaus(poll)^2 + sin(2*theta(poll))/4/sigma_y_gaus(poll)^2 ;
%                             c = sin(theta(poll))^2/2/sigma_x_gaus(poll)^2 + cos(theta(poll))^2/2/sigma_y_gaus(poll)^2;
%                             F(:,:,poll) = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;
%                         end
%                 end
%                 
%                 PrecDummyQuad=squeeze(Prec(ir:ir+rad+rad,ic:ic+rad+rad,iSc,:));
%                 PrecPatch=squeeze(sum(sum(PrecDummyQuad.*F,2),1));
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%% NOT OPTIMIZED OMEGA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %PrecDummyQuad=squeeze(precursors(ir:ir+intermediateResult.rad+intermediateResult.rad,ic:ic+intermediateResult.rad+intermediateResult.rad,iSc,:));
                                PrecDummyQuad=squeeze(precursors(ir:ir+intermediateResult.rad+intermediateResult.rad,ic:ic+intermediateResult.rad+intermediateResult.rad));
                                %PrecPatch=squeeze(sum(sum(PrecDummyQuad.*(repmat(intermediateResult.F,1,1,intermediateResult.nPrec)),2),1)); %IN CASE NOT OPTIMIZED OMEGA
                                %PrecPatch=squeeze(sum(sum(PrecDummyQuad.*(repmat(intermediateResult.F,1,1)),2),1)); %IN CASE NOT OPTIMIZED OMEGA
                                PrecPatch=squeeze(sum(sum(PrecDummyQuad.*(repmat(intermediateResult.F,1,1)),2),1)); %IN CASE NOT OPTIMIZED OMEGA
                                %                 PrecPatch=squeeze(sum(sum(PrecDummyQuad.*F,2),1));   %IN CASE OPTIMIZED OMEGA
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if intermediateResult.pcaFlag==0
                    %finalResult(ir,ic) = PrecPatch'*squeeze(alpha(ir,ic,:));
                    finalResult(ir,ic) = PrecPatch'*squeeze(alpha(ir,ic));
                elseif intermediateResult.pcaFlag==1
                    %                     XNorm=(PrecPatch - squeeze(XMean(ir,ic,:))) ./ squeeze(XStd(ir,ic,:));
                    %                     output(ir,ic) = (XNorm'*squeeze(alpha(ir,ic,:))).*yStd(ir,ic)+yMean(ir,ic);
                    %XNorm=(PrecPatch - squeeze(XMin(ir,ic,:))) ./ ( squeeze(XMax(ir,ic,:)) - squeeze(XMin(ir,ic,:)));
                    XNorm=(PrecPatch - squeeze(XMin(ir,ic))) ./ ( squeeze(XMax(ir,ic)) - squeeze(XMin(ir,ic)));
                    XNorm(isnan(XNorm))=0; %remove nan
                    XNorm(isfinite(XNorm)==0)=0; %remove nan
                    %finalResult(ir,ic) = (XNorm'*squeeze(alpha(ir,ic,:))).* (yMax(ir,ic)-yMin(ir,ic)) + yMin(ir,ic);
                    finalResult(ir,ic) = (XNorm'*squeeze(alpha(ir,ic))).* (yMax(ir,ic)-yMin(ir,ic)) + yMin(ir,ic);
                    
                end
                
        end
    end
    
%finalInfo.iSc=iSc;
%finalInfo.flagRegioMat=flagRegioMat;
finalInfo.resGrid=finalResult;

end

