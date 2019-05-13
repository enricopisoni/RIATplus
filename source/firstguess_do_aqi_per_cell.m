function [aqi_per_cell]= firstguess_do_aqi_per_cell( emissioni, NN, aggregationInfo, commonDataInfo, periodIndex, aqiIndex)

        %ToDO fill: aggregationInfo.alpha;
        %ToDO fill: aggregationInfo.flatweight;
        %use aggregationInfo.firstguess
        
%         %MOD20160607ET
%         flag_optim_dom=commonDataInfo.domainInfo.flag_optim_dom;
%         flag_region_dom=commonDataInfo.domainInfo.flag_region_dom;
%         flag_opt_filt=flag_optim_dom(flag_region_dom==1 | flag_region_dom==2,1);
%         input_rete2=emissioni;
%         if strcmp(NN.Class,'Delta')==0 %if not delta filter on optim
%                 %load network
%                 input_rete2=input_rete2(flag_opt_filt==1 | flag_opt_filt==2,:);
%         end
%         input_rete2=input_rete2';
%         %MOD20160607ET     
        
        flag_optim_dom=commonDataInfo.domainInfo.flag_optim_dom;
        rowToKeep=find(flag_optim_dom==1 | flag_optim_dom==2); %rows related to PAD

        %commonDataInfo.domainInfo.flag_aqi_dom;
        %commonDataInfo.domainInfo.flag_optim_dom;
        alpha=aggregationInfo.extraInfo.alpha;
        omega=aggregationInfo.extraInfo.omega;
        %change dimensions of alpha and omega to be coherent with emissions
        alpha=permute(alpha,[2 1 3]);
        omega=permute(omega,[2 1 3]);
        
        dimX=size(alpha,1);
        dimY=size(alpha,2);
        for i=1:size(alpha,3)
            tmp=alpha(:,:,i);
            thisAlpha(:,i)=reshape(tmp,dimX*dimY,1);
        end
            
        %remove not used data
        thisAlpha(flag_optim_dom==0,:)=[];
%         new_alpha=reshape(alpha,dimy,dimx);
        %emissioni*.aggregationInfo.firstguess.alpha
        %in case it is necessary to process quadrant emissions (if too close
        %to domain boundary, it is necessary to increment emissions with
        %the assumptions that part of the quadrant in which emissions are
        %not available, still contain same emission average)
        dirs=commonDataInfo.dirs;
        domainInfo=commonDataInfo.domainInfo;
        
        %from the SR netcdf, use only NOx(1), NH3(3), PM25(5), SO2(6)
        %if jj eq 0 or 1
        % FF 20190402 Terraria Mod: aggiunto VOC(2) nel calcolo delle
        % concentrazioni
%         if (aqiIndex == 1) aqi_per_cell=sum(emissioni(:,[1 3 4 6]).*thisAlpha(:,[1 3 4 6]),2); end %nox,nh3,pm10,so2
%         if (aqiIndex == 2) aqi_per_cell=sum(emissioni(:,[1 3 5 6]).*thisAlpha(:,[1 3 5 6]),2); end %nox,nh3,pm25,so2
%         if (aqiIndex == 6) aqi_per_cell=sum(emissioni(:,[1]).*thisAlpha(:,[1]),2); end             %nox
        if (aqiIndex == 1) aqi_per_cell=sum(emissioni(:,[1 2 3 4 6]).*thisAlpha(:,[1 2 3 4 6]),2); end %nox,nh3,pm10,so2
        if (aqiIndex == 2) aqi_per_cell=sum(emissioni(:,[1 2 3 5 6]).*thisAlpha(:,[1 2 3 5 6]),2); end %nox,nh3,pm25,so2
        if (aqiIndex == 6) aqi_per_cell=sum(emissioni(:,[1]).*thisAlpha(:,[1]),2); end             %nox
        
        %change name to a better one!!! (too similar to caller...)
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint8_new','input_rete2');
        % 20160421 MM / EP SR / First guess Version
        %aqi_per_cell=interface_get_aqipercell(input_rete2, NN, aggregationInfo.mathIntermediateData, commonDataInfo, periodIndex, aqiIndex );
        aqi_per_cell=firstguess_get_aqipercell(aqi_per_cell, aggregationInfo, commonDataInfo, periodIndex, aqiIndex );
        % 20160421 quadrant / NN version 
        %aqi_per_cell=interface_get_aqipercell(input_rete2, NN, aggregationInfo.mathIntermediateData, commonDataInfo, periodIndex, aqiIndex );
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint9_new','aqi_per_cell');       
        %define if lin o net

end
