function [corr_reg,mse_reg] = CreateScatter(bctarget,target,output,flagRegioMat,iSc,nx,ny,nomeDir,aqi,absdel)
%independent scenario validation scatter

%create scatter for abs and delta values
nameFile={'absValScatter-sce-n','delValScatter-sce-n'};
if absdel==0
    scatterTarget={'target','bctarget-target'};
    scatterOutput={'output','bctarget-output'};
elseif absdel==1
    scatterTarget={'bctarget-target','target'};
    scatterOutput={'bctarget-output','output'};    
end
if strfind(aqi,'PM10')>0
    axisBound={[0 25 0 25],[-5 25 -5 25]}; %for opera
%     axisBound={[0 50 0 50],[-5 40 -5 40]}; %only for riat-lomb
elseif strfind(aqi,'pm10_year_avg')>0
    axisBound={[0 30 0 30],[-5 15 -5 15]};
elseif strfind(aqi,'O3')>0
        axisBound={[0 80 0 80],[-10 10 -10 10]};
elseif strfind(aqi,'o3_year_avg')>0
    axisBound={[0 80 0 80],[-30 30 -30 30]};
elseif strfind(aqi,'PM25')>0
    axisBound={[0 50 0 50],[-5 35 -5 35]};
%     axisBound={[0 20 0 20],[0 10 0 10]}; %only for riat-lomb
elseif strfind(aqi,'NO2')>0
    axisBound={[0 40 0 40],[-5 40 -5 40]};
elseif strfind(aqi,'no2')>0 
    axisBound={[0 40 0 40],[-5 40 -5 40]};
elseif strfind(aqi,'SOMO35')>0
    axisBound={[3000 11000 3000 11000],[-1200 1200 -1200 1200]};
elseif strfind(aqi,'AOT40')>0
    axisBound={[30000 90000 30000 90000],[-20000 20000 -20000 20000]};
elseif strfind(aqi,'MAX8H')>0
    axisBound={[90 130 90 130],[-10 10 -10 10]};
end

for i=1:2
    %scatter
    h=figure;
    xgraph=eval(scatterTarget{i});
    ygraph=eval(scatterOutput{i});
    xgraph(flagRegioMat==0)=[];
    ygraph(flagRegioMat==0)=[];
    xgraph=reshape(xgraph,size(xgraph,1)*size(xgraph,2),1);
    ygraph=reshape(ygraph,size(ygraph,1)*size(ygraph,2),1);
    plot(xgraph,ygraph,'r*');
    axis(axisBound{i});
    grid on;
    hold on;
%         plot([axisBound])
    plot([min(min([xgraph ygraph])) max(max([xgraph ygraph]))],[min(min([xgraph ygraph])) max(max([xgraph ygraph]))],'b--');
%     plot([0 max(max([xgraph ygraph]))/2],[0 max(max([xgraph
%     ygraph]))],'b--');
%     plot([0 max(max([xgraph ygraph]))],[0 max(max([xgraph ygraph]))/2],'b--');
    xlabel('CTM model','FontSize', 20);
    ylabel('SR model','FontSize', 20);
    
    %statistics
    indnan=isnan(ygraph);
    xgraph(indnan==1)=[];
    ygraph(indnan==1)=[];
    
    corr=corrcoef(xgraph,ygraph);
    corr_reg=corr(1,2);
    mse_reg=mean((ygraph-xgraph).^2);
    RMSE = sqrt(mean((ygraph-xgraph).^2));  
    percBias=mean((ygraph-xgraph)./xgraph*100);
    percStd=(std(ygraph)-std(xgraph))/std(xgraph)*100;
    text(2,23,strcat('RMSE=',num2str(RMSE)),'FontSize',16);
    text(2,22,strcat('Corr=',num2str(corr_reg)),'FontSize',16);
    text(2,21,strcat('percBias=',num2str(percBias)),'FontSize',16);
    text(2,20,strcat('percStd=',num2str(percStd)),'FontSize',16);
    
%     title(strcat('corr=',num2str(corr_reg),' rmse=',num2str(RMSE),' percBias=',...
%         num2str(percBias),' percStd=',num2str(percStd)),'FontSize', 14);
    set(gca, 'FontSize', 16);
    
    %save
    nameScatter=strcat(nomeDir,nameFile{i},int2str(iSc));
    print(h,'-dpng',nameScatter);
    close all
end



