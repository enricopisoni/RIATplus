function [] = CreateMap(bctarget,target,output,flagRegioMat,x,y,iSc,nomeDir,aqi,absdel,flagReg,domain)
%creating maps

target(flagRegioMat==0)=NaN;
output(flagRegioMat==0)=NaN;

%create graphs for target, output, bias, percentage bias
nameFile={'targetMap-sce-n','outputMap-sce-n','biasMap-sce-n','biasPercMap-sce-n','output-potencial-sce-n','target-potencial-sce-n'};
if absdel==0
    mapInfo={'target','output','output-target','(output-target)./target*100','output_DeltaC','target_DeltaC'};
elseif absdel==1
    mapInfo={'bctarget-target','bctarget-output',...
        '(bctarget-output)-(bctarget-target)',...
        '((bctarget-output)-(bctarget-target))./(bctarget-target)*100','output*2','target*2'};
end
if strfind(aqi,'PM10')>0
    %     levels={[0:1:16],[0:1:16],[-5:0.5:5],[-10:1:10]}; %for opera
    levels={[0:1:25],[0:1:25],[-2:0.2:2],[-5:0.5:5]}; %only for riat-lomb
    Range={[0 25],[0 25],[-2 2],[-5 5]};
elseif strfind(aqi,'pm10_year_avg')>0
    levels={[0:1:24],[0:1:24],[-5:0.5:5],[-10:1:10]};
elseif strfind(aqi,'o3_year_avg')>0
    levels={[0:2:80],[0:2:80],[-10:1:10],[-20:1:20]};
elseif strfind(aqi,'O3')>0
    levels={[0:2:80],[0:2:80],[-10:1:10],[-20:1:20]};
    Range={[30 50],[30 50],[-2 2],[-10 10]};
elseif strfind(aqi,'PM25')>0
    levels={[0:1:20],[0:1:20],[-2:0.2:2],[-15:1:15]};
    Range={[0 25],[0 25],[-2 2],[-10 10],[0 10],[0 10]};
elseif strfind(aqi,'NO2')>0
    levels={[0:1:20],[0:1:20],[-5:0.5:5],[-30:3:30]};
elseif strfind(aqi,'no2')>0
    levels={[0:1:20],[0:1:20],[-5:0.5:5],[-30:3:30]};
elseif strfind(aqi,'SOMO35')>0
    levels={[3000:500:11000],[3000:500:11000],[-400:50:400],[-8:1:8]};
elseif strfind(aqi,'AOT40')>0
    levels={[30000:5000:90000],[30000:5000:90000],[-10000:1000:10000],[-8:1:8]};
elseif strfind(aqi,'MAX8H')>0
    levels={[90:2:130],[90:2:130],[-5:1:5],[-4:0.5:4]};
end


for i=1:6
    h=figure;
    %     contourf(x,y,eval(mapInfo{i}));
        
    if strcmp(domain,'tno_data')==1
        flagRegRed=importdata(strcat('./input/tno_data/flag_',flagReg,'.mat'));
        contour(x,y,flagRegRed,1,'k');
        hold on
        geoshow(y,x,eval(mapInfo{i}),'DisplayType','texturemap');
        contourcmap(levels{i},'jet','colorbar','on');
    elseif strcmp(domain,'50km_chimere')==1
        %create figures
        h=figure;
        %europe map
        worldmap('Belgium');%worldmap('Belgium');%worldmap('Europe');
        %country boundaries
        geoshow(strcat('./input/',domain,'/Cntry02/cntry02.shp'),'FaceColor','white');        %map to be plotted
        geoshow(y,x,eval(mapInfo{i}),'DisplayType','texturemap');
        %title and colorbar
        colorbar
        %fix bounds for visualization
        set(gca,'CLim',Range{i});
        %transparency, and NaN in white
        %         alpha(.8);
        %         ColorData = get(h,'Cdata');
        %         set(h,'AlphaData',...
        %             double(~isnan(ColorData)),'FaceAlpha',...
        %             'texturemap','FaceColor','texture',...
        %             'AlphaDataMapping','none');
    elseif (strcmp(domain,'ineris_28km')==1 || strcmp(domain,'ineris_7km_9bins')==1 ...
            || strcmp(domain,'ineris_7km_20150731_fixedSce6')==1 || strcmp(domain,'ineris_7km_deliver_20150831')==1)
        %create figures
        h=figure;
        %europe map
%         worldmap('Europe');
        %country boundaries
        geoshow(y,x,eval(mapInfo{i}),'DisplayType','texturemap');
        geoshow(strcat('./input/',domain,'/Cntry02/cntry02.shp'),'FaceColor','none');        %map to be plotted
%         geoshow(strcat('./input/',domain,'/NUTS_2013_SHP/data/NUTS_BN_01M_2013.shp'),'FaceColor','none');        %map to be plotted
        tmp=gca;
        tmp.XLim=[min(min(x)) max(max(x))];
        tmp.YLim=[min(min(y)) max(max(y))];

        %title and colorbar
        colorbar
        %fix bounds for visualization
        set(gca,'CLim',Range{i});

    else
%         worldmap('italy');
        worldmap([43.5 45.5],[8.5 13.5]);
        geoshow('./input/tno_data/_shape_nuts/NUTS_2010_03M_SH/Data/NUTS_BN_03M_2010.shp');
        for ii=1:size(x,1)
            for jj=1:size(x,2)
                [Lat(ii,jj) Lon(ii,jj)]=utm2deg(x(ii,jj)+2500,y(ii,jj)+2500,'32 T');
            end
        end
        geoshow(Lat,Lon,eval(mapInfo{i}),'DisplayType','texturemap');
        colorbar
        %fix bounds for visualization
        set(gca,'CLim',Range{i});
        alpha(.8);
%         contourf(x,y,eval(mapInfo{i}));
%         contourcmap(levels{i},'jet','colorbar','on');
    end
%     shading flat
    nameScatter=strcat(nomeDir,nameFile{i},int2str(iSc));
    print(h,'-dpng',nameScatter);
    
end
