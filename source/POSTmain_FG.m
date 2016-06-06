%20160517
%increased dlmwrite precision, to 5 decimals, file maps_cost
%always write 11 macrosectors, in line 671

% new (temporary) post to manage First Guess/SR mode
function []=POSTmain_FG(f2,DSuperSet,aggregationInfo,pathANN,flag_region_dom,flag_optim_dom,flag_aqi_dom,...
    coordinate,emi,global_data,base_emi_low,base_emi_high,CLE,...
    base_emi_low_noc,base_emi_high_noc,pathOCM,sSet,areal_point,...
    flag_reg_net,pathOUF,pop,nocID,pathAR,pm10aveToExceed,...
    cell_threshold_set,fidStatus,fidExit, commonDataInfo)


%UNIBS(ET)20131001 - added flag_region_dom in input
%20130417 - ADDED THE "MODEL BIAS" CAPABILITY, TO CORRECT THE AQIS USING
%THE CTM MODEL BIAS
%postprocessing message

strStatus='PROGRESSION: Starting the post-processing...';
disp(strStatus);
fprintf(fidStatus, '%s\n',strStatus);
% fprintf(fidStatus, '%s\n','Starting the post-processing...');

%POST PROCESSING OF RESULTS
format short

% %number of considered GHG=4;
nghg=4;

%UNIBS(ET)20131001 -define number of cells of regional domain
%to compute AQI for maps
%ncelopt=length(find(flag_region_dom==1 | flag_region_dom==2));
flag_optim_dom=commonDataInfo.domainInfo.flag_optim_dom;
ncelopt=length(find(flag_optim_dom==1 | flag_optim_dom==2));

%number of pareto points
npoints=size(sSet,2);

%load file containing all the paths
fid=fopen(f2);
c1=textscan(fid,'%s');
detailresults=cell2mat(c1{1}(1));       %0=no figure; 1=figures
detailresults=str2num(detailresults);   %trasform it to number
pathghg=cell2mat(c1{1}(2));             %ghg db file
pathextpop=cell2mat(c1{1}(3));          %population data to compute external costs
pathextoth=cell2mat(c1{1}(4));          %coefficient to compute external costs
pathprb=cell2mat(c1{1}(5));             %regional boundaries file
pathpci=cell2mat(c1{1}(6));             %cities file
pathMB=cell2mat(c1{1}(7));              %model bias
fclose(fid);

%define some variables to save results
nomeaqi={'PM10','PM25','AOT40','SOMO35','MAX8H','NO2','PM10dailyExceed'};
nomeagg={'_MEAN','_THRE','_WEIGMEAN'};
nomesea={'TP1','TP2','TP3'};
ylab={'[microg/m3]','[microg/m3]','[microg/m3*h]','[microg/m3*d]','[microg/m3]','[microg/m3]','[# days]';...
    '[#overThres]','[#overThres]','[#overThres]','[#overThres]','[#overThres]','[#overThres]','[#overThres]';...
    '[microg/m3]','[microg/m3]','[microg/m3*h]','[microg/m3*d]','[microg/m3]','[microg/m3]','[# days]'};

%ENR20130417 - create data structure for model bias
model_bias_data=importdata(pathMB);
model_bias=struct('AQIs',{});

%20150304 ET - MODEL BIAS COLUMN CORRECTION
cols=size(model_bias_data.data,2);
if cols<23
    model_bias_data.data(:,cols+1:23)=0;
end
for indexSea = 1:size(nomesea,2), %available season
    for j = 1:size(pathANN(1).ANNs,1), %aqis
        %UNIBS(ET)20131001 - changed  flag_optim_dom with flag_region_dom
        rowToKeep=find(flag_region_dom==1 | flag_region_dom==2); %rows related to PAD
        model_bias(indexSea).AQIs{j}=model_bias_data.data(rowToKeep,2+j+7*(indexSea-1));
    end
end

%check if the regional boundaries file / cities file have been provided
%if one of the two files is missing, put detailresults=0
if strcmp(pathprb,'-1')
    detailresults=0;
end
if strcmp(pathprb,'-1')
    detailresults=0;
end

pathopp=pathOUF;

%load data for external costs
%evaluate the case in which external costs are not computed
if (strcmp(pathextpop,'-1') | strcmp(pathextoth,'-1'))
    extpop=-999;
    extoth=-999;
else
    extpopTmp=importdata(pathextpop);
    extpop=extpopTmp.data(:,3:end);
    extoth=importdata(pathextoth);
end
% 	- "external_cost_pop_data"contains
% 		% pop per age, in % (9 columns: 0-5 6-10 11-14 15-19 20-24 25-29 30-59 60-64 >65)
% 		% mortality rate per age (10-18 columns: 0-5 6-10 11-14 15-19 20-24 25-29 30-59 60-64 >65)
% 		% asmathics rate in % (19 column)
% 		% fraction of population > 30 years (20 column)
% 	- "external_cost_oth_data": contains impact and cost coefficients (2 columns)
% 		ASMATICI
% 		Adulti
% 		 	1  uso di broncodilatatore
% 			2  tosse
% 			3  affaticamento nella respirazione (sintomi)
% 		Bambini
% 			4  uso di broncodilatatore
% 			5  tosse
% 			6  affaticamento nella respirazione (sintomi)
% 		OLTRE I 65 ANNI
% 			7  infarto
% 		BAMBINI
% 			8  tosse cronica
% 		ADULTI
% 			9  giorni di attività ridotta
% 			10 bronchite cronica
% 		POPOLAZ. TOT
% 			11 mortalità cronica
% 			12 ricoveri per problemi respiratori
% 			13 ricoveri per problemi cardiocircolatori
% 		SOPRA I 30 ANNI
% 			14 anni di vita persi

% define domain coordinate
xutm=coordinate(:,1);
yutm=coordinate(:,2);

%create directories for saving results
warning off
mkdir(pathopp)
outdir=pathopp;
mkdir(strcat(pathopp,'\emi_cost'));
mkdir(strcat(pathopp,'\emi_cost\m'));
mkdir(strcat(pathopp,'\emi_cost\msat'));
mkdir(strcat(pathopp,'\maps_aqi'));
mkdir(strcat(pathopp,'\maps_cost'));
mkdir(strcat(pathopp,'\pareto'));
mkdir(strcat(pathopp,'\maps_emi'));

%GHG BUDGET COMPUTATION
for indexGHG=1:npoints
    xsolv=sSet(indexGHG).X;
    %writing CO2, CH4, N20, Fgas [kton/year]
    
    %if file for GHG is not provided, do not compute GHG
    if strcmp(pathghg,'-1')
        emi_rem_tot(indexGHG,1:nghg)=-999;
    else
        
        emi_rem_tot(indexGHG,:)=POSTcomputeGHG(xsolv,indexGHG,outdir,...
            pathOCM,CLE,global_data,pathghg,emi);
    end
end
budgetCO2=[emi_rem_tot];

%COMPUTE AQIS (SEASONAL AND YEARLY) AND EMISSION (YEARLY) MAPS
%define levels to create aqi maps
% OPERA ALSACE
%livelliAQI={0:1:18,0:1:16,0:2000:40000,3000:200:6200,80:2:110,0:2:36,0:1:50};
% RIAT
livelliAQI={0:5:45,0:5:45,10000:5000:90000,1000:250:6000,0:5:45,0:5:45,0:1:50};
% OPERA EMR
% livelliAQI={0:1:25,0:1:20,30000:2500:70000,5000:250:9000,90:2.5:120,0:2.5:35,0:2.5:35};

% managing the fact that PAD and ACD are different
indpad=(flag_optim_dom==1 | flag_optim_dom==2);
indacd=flag_aqi_dom(indpad);


%create aqi maps
for indexPareto=1:npoints
    for indexSea = 1:size(nomesea,2), %available season
        for j = 1:size(pathANN(1).ANNs,1), %aqis
            %compute the requested AQI, only if related ANN is available
            if (isequal(strtrim(pathANN(indexSea).ANNs(j,:)),'-999')==0)
                %display stage of the work
                disp(strcat('saving aqi maps, aqi',nomeaqi{j},...
                    '-season',nomesea{indexSea},'-point ',int2str(indexPareto)));
                %compute emission and aqi maps
                %20140403 (ET) - Added bcSuperSet as input
                % inside use
                % emis=FG_rebuildOrderedEmis(E_full, aggregationInfo.geometryIntermediateData, commonDataInfo);
                %20160422 (MM) - Quadrant/NN Version
                %[emi_low,emi_high,aqi_val]=POSTcompute_aqi(sSet(indexPareto).X,nnSuperSet(indexSea).nnSet(j),...
                %    DSuperSet(indexSea).DSet(j),bcSuperSet(indexSea).bcSet(j));
                SRInfo=aggregationInfo.firstguess;
                % calling:
                % function [emi_low,emi_high,aqi_val2]=POST_FG_compute_aqi(x,SR,DD,indexSea, aqiIndex)
                [emi_low,emi_high,aqi_val]=POST_FG_compute_aqi(sSet(indexPareto).X,SRInfo,...
                    DSuperSet, indexSea, j, emi, commonDataInfo);
                %DSuperSet(indexSea).DSet(j), indexSea, j);
                
                %ENR20130417 - adding the model bias
                aqi_val=aqi_val'+model_bias(indexSea).AQIs{j};
                
                %compute external costs only for PM10, and for annual
                %values
                extCost(indexPareto)=-999*1000000; % 10^6 is used because of Meuro conversion;
                extCost_morb(indexPareto)=-999*1000000;
                extCost_yoll(indexPareto)=-999*1000000;
                if j==1 && indexSea==1
                    
                    %if there are no external cost data, do not perform
                    %computations
                    if extpop==-999
                        extCost(indexPareto)=-999*1000000; % 10^6 is used because of Meuro conversion;
                        extCost_morb(indexPareto)=-999*1000000;
                        extCost_yoll(indexPareto)=-999*1000000;
                    else
                        %UNIBS(ET)20131001 - added  flag_region_dom in
                        %input
                        [extCost_tot_morb,extCost_tot_yoll, resp_prob,card_prob,yoll_prob]=...
                            POSTexternalcosts(aqi_val,extpop,extoth,pop,flag_optim_dom,flag_region_dom);
                        
                        extCost_morb(indexPareto)=extCost_tot_morb;
                        extCost_yoll(indexPareto)=extCost_tot_yoll;
                        
                        %save some impacts infos
                        %                         POSTcomputeMaps(0:5:40,resp_prob,...
                        %                             strcat('maps_aqi/impact_respiratory',nomesea{indexSea}),strcat('impact_respiratory',nomesea{indexSea},{'cases/year'}),...
                        %                             indexPareto,outdir,pathprb,pathpci,flag_optim_dom,detailresults);
                        
                        %                         POSTcomputeMaps(0:5:60,card_prob,...
                        %                             strcat('maps_aqi/impact_cardiovasculary',nomesea{indexSea}),strcat('impact_cardiovasculary',nomesea{indexSea},{'cases/year'}),...
                        %                             indexPareto,outdir,pathprb,pathpci,flag_optim_dom,detailresults);
                        %UNIBS(ET)20131002 - flag_region_dom in input used
                        %to print the map for the entire regional domain
                        POSTcomputeMaps(0:1:12,yoll_prob,...
                            strcat('maps_aqi/impact_yoll',nomesea{indexSea}),strcat('impact_yoll',nomesea{indexSea},{'cases/year'}),...
                            indexPareto,outdir,pathprb,pathpci,flag_region_dom,detailresults);
                    end
                end
                
                if(j == 7) %case for PM10 number of exceedances over threshold
                    % EU level relation, fron Philippe Thunis, from Chimere/Urban stations
                    aqi_val = ...
                        pm10aveToExceed(1).*aqi_val-pm10aveToExceed(2);
                    aqi_val (aqi_val < 0) = 0;
                end
                
                %create map of aqi
                %UNIBS(ET)20131002 - flag_region_dom in input used
                %to print the map for the entire regional domain
                POSTcomputeMaps(livelliAQI{j},aqi_val,...
                    strcat('maps_aqi/',nomeaqi{j},nomesea{indexSea}),strcat(nomeaqi{j},nomesea{indexSea},ylab(1,j)),...
                    indexPareto,outdir,pathprb,pathpci,flag_region_dom,detailresults);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %UPDATE SSET WITH MODEL BIAS
                aqi_per_cell=aqi_val;
                aqi_per_cell(indacd==0)=[];
                %PATCH 20160601 EP
                indisnan=find(isnan(aqi_per_cell));
                if isempty(indisnan)==0
                    aqi_per_cell(isnan(aqi_per_cell))=[];
                end
                
                %structure of sSet:
                %sSet(npoints).AQIs(yea-win-sum, AQInum, mean-thre-pwa)
                %average
                sSet(indexPareto).AQIs(indexSea,j,1)=mean(aqi_per_cell);
                % number of cells over thresold
                T = ones(size(aqi_per_cell)) * cell_threshold_set(j);
                sSet(indexPareto).AQIs(indexSea,j,2)=sum(aqi_per_cell > T);
                % population average
                indpop=pop;
                %UNIBS(ET)20131002 - Erase cells outside regional domain
                indpop(flag_region_dom==0)=[];
                indpop(indacd==0)=[];
                %PATCH 20160601 EP                
                if isempty(indisnan)==0
                    indpop(indisnan)=[];
                end
                sSet(indexPareto).AQIs(indexSea,j,3)=sum(aqi_per_cell.*indpop)/sum(indpop);
                %END OF MODEL BIAS UPDATE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end
        end
        
        %depict emission maps - 6 precursors
        %         livelliemi={0:100:2500,0:100:3000,0:50:800,0:10:100,0:10:100,0:100:2000}; %RIAT
        livelliemi={0:5:150,0:5:150,0:5:100,0:2.5:60,0:2.5:60,0:2.5:100}; %OPERA-aspa
        nomeemi={'NOX','VOC','NH3','PM10','PM25','SO2'};
        nomeUM={'[ton/year]'};
        %         if detailresults==1
        if indexSea==1 %only for yearly values, create emissions
            for poll=1:6
                %display stage of work
                disp(strcat('saving emission maps, emi',nomeemi{poll},'-point ',int2str(indexPareto)));
                %create map
                %UNIBS(ET)20131002 - flag_region_dom in input used
                %to print the map for the entire regional domain
                POSTcomputeMaps(livelliemi{poll},emi_low(:,poll)+emi_high(:,poll),...
                    strcat('maps_emi/',nomeemi{poll},nomesea{indexSea}),strcat(nomeemi{poll},nomesea{indexSea},nomeUM),...
                    indexPareto,outdir,pathprb,pathpci,flag_region_dom,detailresults);
            end
        end
    end
end

%External costs in Meuro
extCost_morb=extCost_morb/1000000;
extCost_yoll=extCost_yoll/1000000;

%CREATE PARETO CURVE AS TXT FILE, ADDING CO2 BUDGET AND EXTERNAL COSTS
[cost_final_over_cle,aqi_final]=POSTparetocurve(sSet,budgetCO2,extCost_morb,extCost_yoll);

%CREATE PARETO CURVE AS FIGURES (IF REQUIRED)
if detailresults==1 %in case detailed results are requestd
    for k=1:size(nomeagg,2) %two types of taggregation (mean or threshold)
        for indexSea = 1:size(nomesea,2), %available season (year, winter, summer)
            for j = 1:size(pathANN(1).ANNs,1), %7 aqis defined
                %compute the requested AQI, only if related ANN is available
                if (isequal(strtrim(pathANN(indexSea).ANNs(j,:)),'-999')==0)
                    %                     if k==4 & j>1 %for the aggregation method related to dailyPM10 threshold exceedances, only consider PM10 ANNs
                    %                         continue
                    %                     end
                    x=0:cost_final_over_cle(end);
                    %index to read in the "aqi_final" variable
                    %7 are the AQIs, 21 are the AQIs*temporal periods
                    %dimensionss
                    indread=j+7*(indexSea-1)+21*(k-1);
                    %create figures
                    %choose to depict only different costs (case of lot of
                    %points on pareto curve)
                    [a , difcost, c]=unique(round(cost_final_over_cle),'rows');
                    y = interp1(cost_final_over_cle(difcost),aqi_final(difcost,indread),x,'pchip');
                    plot(x,y);
                    set(gca,'fontsize',18)
                    hold on;
                    plot(cost_final_over_cle(difcost),aqi_final(difcost,indread),'bo');
                    set(gca,'fontsize',18)
                    xlabel('cost [meuro]');
                    ylabel(strcat({'Air quality index '},ylab(k,j)));
                    hold off;
                    print(gcf,'-dpng',strcat('./',outdir,'/pareto/',strcat(nomeaqi{j},nomeagg{k},nomesea{indexSea}),'.png'));
                end
            end
        end
    end
end

% %COMPUTE EMISSION AND COST TABLES, FOR EACH PARETO CURVE POINT
for indexPareto=1:npoints
    xsolv=sSet(indexPareto).X;
    disp(strcat('saving emission and cost details, point ',int2str(indexPareto)));
    POSTcomputeTables(xsolv,indexPareto,outdir,pathOCM,CLE,global_data)
end

strStatus='PROGRESSION: Post-procesing part finished.';
disp(strStatus);
fprintf(fidStatus, '%s\n',strStatus);
% fprintf(fidStatus, '%s\n','Post-procesing part finished.');
fclose(fidStatus)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FUNCTION TO COMPUTE PARETO CURVE
    function [cost_final_over_cle,aqi_final]=POSTparetocurve(sSet,budgetCO2,extCost_morb,extCost_yoll)
        
        %number of pareto points
        npoints=size(sSet,2);
        
        %variables initialization
        cost_final_over_cle=[];
        aqi_final=[];
        
        %solution set is written in sSet variable
        %         sSet =  struct('X', {}, 'X_free', {}, 'AQIs', {}, ...
        %             'COST', {}, 'BUDGET', {}, ...
        %             'COSTPERMACROSECTOR', {}, ...
        %             'BUDGETPERMACROSECTOR', {});
        %         for i=1:npoints
        %             sSet{i}=solutionSet(i).SOL_LIST;
        %         end
        
        %for each point create vector of elements, for cost and aqi
        %in AQI_final, write a line for each pareto point, and a column for each of
        %the following AQIs
        
        %yearly-mean: pm10, pm25, aot, somo, max8h, n02, pm10DailyExceed
        %winter-mean: pm10, pm25, aot, somo, max8h, n02, pm10DailyExceed
        %summer-mean: pm10, pm25, aot, somo, max8h, n02, pm10DailyExceed
        %yearly-thre: pm10, pm25, aot, somo, max8h, n02, pm10DailyExceed
        %winter-thre: pm10, pm25, aot, somo, max8h, n02, pm10DailyExceed
        %summer-thre: pm10, pm25, aot, somo, max8h, n02, pm10DailyExceed
        %yearly-weigmean: pm10, pm25, aot, somo, max8h, n02, pm10DailyExceed
        %winter-weigmean: pm10, pm25, aot, somo, max8h, n02, pm10DailyExceed
        %summer-weigmean: pm10, pm25, aot, somo, max8h, n02, pm10DailyExceed
        
        % aqi(ii,jj,kk)
        % ii=season
        % jj=aqi
        % kk=aggregation type
        
        for i=1:npoints
            cost_final_over_cle=[cost_final_over_cle; ...
                sSet(i).COST-sSet(1).COST];
            d1=size(sSet(i).AQIs(:,:,1),1);
            d2=size(sSet(i).AQIs(:,:,1),2);
            aqi_mean=reshape(sSet(i).AQIs(:,:,1)',1,d1*d2);
            aqi_thre=reshape(sSet(i).AQIs(:,:,2)',1,d1*d2);
            aqi_weigmean=reshape(sSet(i).AQIs(:,:,3)',1,d1*d2);
            
            aqi_final=[aqi_final; aqi_mean aqi_thre aqi_weigmean];
        end
        
        res=[cost_final_over_cle aqi_final budgetCO2 extCost_morb' extCost_yoll'];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %MANAGING OPTIMIZATION ERRORS
        %Test over solution costs
        TestData=res(:,1);
        Check=ones(size(res,1),1);
        for i_test=1:size(res,1)-1
            Check(i_test+1)=TestData(i_test)<TestData(i_test+1);
        end
        
        if  TestData(1)==TestData(end)
            display('There is no possible solution for the selected configuration.')
            Error=1;
        elseif sum(Check)~=size(res,1)
            
            %   for j_test=2:size(res,1)
            %       if Check(j_test)==0
            %           res(j_test,:)=[];
            %       end
            %   end
            
            display('For at least one of the points of the pareto curve, the solution is not robust.')
            Error=2;
        else
            display('The solution is robust')
            Error=0;
        end
        fprintf(fidExit, int2str(Error));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %write results on CSV file, aggregated values to depict pareto curve
        disp('saving pareto curve')
        file=strcat(outdir,'/pareto/resultscostvsaqis.csv');
        fid=fopen(file,'wt');
        delim=',';
        %fprintf(fid,'%s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s\n',...
        fprintf(fid,'%s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s\n',...
            'Cost Over CLE [Meuro]', delim, 'pm10_mean_TP1 [microg/m3]', delim, ...
            'pm25_mean_TP1 [microg/m3]', delim, 'aot40_mean_TP1 [microg/m3*h]', delim, ...
            'somo35_mean_TP1 [microg/m3*d]', delim, 'max8h_mean_TP1 [microg/m3*h]', delim, ...
            'no2_mean_TP1 [microg/m3*d]', delim, 'pm10_dailyExceed_TP1 [#]', delim,...
            'pm10_mean_TP2 [microg/m3]', delim, ...
            'pm25_mean_TP2 [microg/m3]', delim, 'aot40_mean_TP2 [microg/m3*h]', delim, ...
            'somo35_mean_TP2 [microg/m3*d]', delim, 'max8h_mean_TP2 [microg/m3*h]', delim, ...
            'no2_mean_TP2 [microg/m3*d]', delim, 'pm10_dailyExceed_TP2 [#]', delim,...
            'pm10_mean_TP3 [microg/m3]', delim, ...
            'pm25_mean_TP3 [microg/m3]', delim, 'aot40_mean_TP3 [microg/m3*h]', delim, ...
            'somo35_mean_TP3 [microg/m3*d]', delim, 'max8h_mean_TP3 [microg/m3*h]', delim, ...
            'no2_mean_TP3 [microg/m3*d]', delim, 'pm10_dailyExceed_TP3 [#]', delim,...
            'pm10_thre_TP1 [microg/m3]', delim, ...
            'pm25_thre_TP1 [microg/m3]', delim, 'aot40_thre_TP1 [microg/m3*h]', delim, ...
            'somo35_thre_TP1 [microg/m3*d]', delim, 'max8h_thre_TP1 [microg/m3*h]', delim, ...
            'no2_thre_TP1 [microg/m3*d]',delim,'pm10_dailyExceed_thre_TP1 [#]', delim,...
            'pm10_thre_TP2 [microg/m3]', delim, ...
            'pm25_thre_TP2 [microg/m3]', delim, 'aot40_thre_TP2 [microg/m3*h]', delim, ...
            'somo35_thre_TP2 [microg/m3*d]', delim, 'max8h_thre_TP2 [microg/m3*h]', delim, ...
            'no2_thre_TP2 [microg/m3*d]',delim,'pm10_dailyExceed_thre_TP2 [#]', delim,...
            'pm10_thre_TP3 [microg/m3]', delim, ...
            'pm25_thre_TP3 [microg/m3]', delim, 'aot40_thre_TP3 [microg/m3*h]', delim, ...
            'somo35_thre_TP3 [microg/m3*d]', delim, 'max8h_thre_TP3 [microg/m3*h]', delim, ...
            'no2_thre_TP3 [microg/m3*d]',delim,'pm10_dailyExceed_thre_TP3 [#]', delim,...
            'pm10_weigmean_TP1 [microg/m3]', delim, ...
            'pm25_weigmean_TP1 [microg/m3]', delim, 'aot40_weigmean_TP1 [microg/m3*h]', delim, ...
            'somo35_weigmean_TP1 [microg/m3*d]', delim, 'max8h_weigmean_TP1 [microg/m3*h]', delim, ...
            'no2_weigmean_TP1 [microg/m3*d]', delim, 'pm10_dailyExceed_weigmean_TP1 [#]', delim,...
            'pm10_weigmean_TP2 [microg/m3]', delim, ...
            'pm25_weigmean_TP2 [microg/m3]', delim, 'aot40_weigmean_TP2 [microg/m3*h]', delim, ...
            'somo35_weigmean_TP2 [microg/m3*d]', delim, 'max8h_weigmean_TP2 [microg/m3*h]', delim, ...
            'no2_weigmean_TP2 [microg/m3*d]', delim, 'pm10_dailyExceed_weigmean_TP2 [#]', delim,...
            'pm10_weigmean_TP3 [microg/m3]', delim, ...
            'pm25_weigmean_TP3 [microg/m3]', delim, 'aot40_weigmean_TP3 [microg/m3*h]', delim, ...
            'somo35_weigmean_TP3 [microg/m3*d]', delim, 'max8h_weigmean_TP3 [microg/m3*h]', delim, ...
            'no2_weigmean_TP3 [microg/m3*d]',delim,'pm10_dailyExceed_weigmean_TP3 [#]', delim,...
            'CO2 [kton/year]',delim,'CH4 [kton/year]',...
            delim,'N20 [kton/year]',delim,'Fgases [kton/year]',delim,...
            'External Cost Morbidity [Meuro/year]',delim,'External Cost Yoll [Meuro/year]');
        fclose(fid);
        dlmwrite(file,res,'-append','roffset', 0, 'precision', 5);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [emi_rem_tot]=POSTcomputeGHG(x1,indexGHG,outdir,pathOCM,CLE,global_data,pathghg,emi)
        
        resTMP=0;
        
        %load files with technologies DB
        out_clemfr=load(pathOCM);
        out_clefr_ghg=load(pathghg);
        
        %load infos for GHG
        reGHG=out_clefr_ghg(:,6:9);
        unefGHG=out_clefr_ghg(:,10:13);
        
        %inizialization of variables
        emicosttab_msat_tot=0;
        emicosttab_msa_tot=0;
        emicosttab_m_tot=0;
        emi_alpha_noc_tot=0;
        emi_rem_tot=0;
        %cv as a vector
        if (size(x1,1)>size(x1,2))
            x=x1;
        else
            x=x1';
        end
        
        %define low emissions matrix (to be used for emission computation)
        x_alpha=[x x x x];
        
        %define CLE matrix, for results computation
        x_alphaCLE=[CLE CLE CLE CLE];
        result=[];
        icel=0;
        
        %put very small numbers to zero (approximation issue)
        x(find(x<1e-6))=0;
        
        %add noc sec-act-tec data
        tmp=out_clemfr(find(out_clemfr(:,4)==nocID),:);
        global_data=[global_data out_clemfr(find(out_clemfr(:,4)~=nocID),24) ...
            out_clemfr(find(out_clemfr(:,4)~=nocID),25); tmp];
        
        xtmp=zeros(size(tmp,1),4);
        x_alpha=[x_alpha; xtmp];
        x_alphaCLE=[x_alphaCLE; xtmp];
        xtmp=zeros(size(tmp,1),1);
        x=[x; xtmp];
        CLE=[CLE; xtmp];
        
        %select sector-activity pairs to be considered for aggregated results
        tmp_sec_act=global_data(:,[2 3]);
        [sec_act bsa c]=unique(tmp_sec_act,'rows');
        
        %select macrosectors to be considered for aggregated results
        tmp_ms=global_data(:,[1]);
        [ms bms c]=unique(tmp_ms,'rows');
        
        %unit costs
        %         uc=global_data(:,18);
        
        %loop over cells to create emissions file
        for i=1:length(flag_optim_dom)
            i;
            re=[];
            emi_alpha=[];
            
            %for cells inside the domain
            if (flag_optim_dom(i)==1 | flag_optim_dom(i)==2)
                icel=icel+1;
                
                %coordinate
                xutm=coordinate(i,1)*1000;
                yutm=coordinate(i,2)*1000;
                
                %load emissions
                emi_alpha=[];
                emi_alpha=emi{i,2};
                
                actlev=emi_alpha(:,7);
                emi_alpha_fin=[];
                emi_alpha_fin=repmat(actlev,1,4).*unefGHG;
                re=reGHG/100; %removal efficiencies
                
                %emission reductions for single cells, detailed description
                emi_red_msat_cel_opt{i}=repmat(actlev,1,4).*unefGHG.*re.*(x_alpha);
                emi_red_msat_cel_cle{i}=repmat(actlev,1,4).*unefGHG.*re.*(x_alphaCLE);
                
                %results for single cells, aggregated for sec-act
                %do not consider difference between low and high
                for j=1:size(sec_act,1)
                    ind=find(global_data(:,2)==sec_act(j,1) & global_data(:,3)==sec_act(j,2));
                    ind1=find(global_data(:,2)==sec_act(j,1) & global_data(:,3)==sec_act(j,2) & global_data(:,5)==1);
                    ind2=find(global_data(:,2)==sec_act(j,1) & global_data(:,3)==sec_act(j,2) & global_data(:,5)==2);
                    
                    emi_red_msa_cel_opt{i}(j,1:4)=sum(emi_red_msat_cel_opt{i}(ind,:),1);
                    emi_red_msa_cel_cle{i}(j,1:4)=sum(emi_red_msat_cel_cle{i}(ind,:),1);
                    
                    %compute initial emissions, change if you take into account low or high level emissions
                    if (~isempty(ind1) & ~isempty(ind2))
                        emi_ini_msa_cel{i}(j,1:4)=emi_alpha_fin(ind1(1),:)+emi_alpha_fin(ind2(1),:);
                    elseif (~isempty(ind1) & isempty(ind2))
                        emi_ini_msa_cel{i}(j,1:4)=emi_alpha_fin(ind1(1),:);
                    elseif (isempty(ind1) & ~isempty(ind2))
                        emi_ini_msa_cel{i}(j,1:4)=emi_alpha_fin(ind2(1),:);
                    elseif (isempty(ind1) & isempty(ind2))
                        emi_ini_msa_cel{i}(j,1:4)=0;
                    end
                    
                    %if considering the remaining emissions (ini-red)
                    %emi_rem_msa_cel{i}(j,:)=emi_ini_msa_cel{i}(j,:)-emi_red_msa_cel_opt{i}(j,:);
                    %if considering reduced emissions over cle (opt-cle)
                    emi_rem_msa_cel{i}(j,:)=emi_red_msa_cel_opt{i}(j,:)-emi_red_msa_cel_cle{i}(j,:);
                    %emi_rem_msa_cel{i}(j,:)=emi_red_msa_cel_opt{i}(j,:);
                end
                
                resTMP=resTMP+emi_alpha_fin;
                emi_rem_tot=emi_rem_tot+sum(emi_rem_msa_cel{i});
            end
        end
        
        emi_rem_msa_cel;
        
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FUNCTION TO COMPUTE EMISSIONS AND COSTS TABLES
    function []=POSTcomputeTables(x,index,outdir,pathOCM,CLE,global_data)
        %,flag_optim_dom,coordinate,...
        %   emi)
        %compute emission and cost tables, aggregated of ms, ms-s-a, ms-s-a-t
        
        %v1.1.2 EP - now the variable out_clemfr is correctly defined.
        %otherwise it was not working, if feeded with a file not named
        %"out_clemfr.txt"
        %measures DB
        out_clemfr=load(pathOCM);
        
        %format of data
        %         format bank
        
        %inizialization of variables
        emicosttab_msat_tot=0;
        emicosttab_msa_tot=0;
        emicosttab_m_tot=0;
        emi_alpha_noc_tot=0;
        
        %define low emissions matrix (to be used for emission computation)
        x_alpha=[x x x x x x];
        
        %define CLE matrix, for results computation
        x_alphaCLE=[CLE CLE CLE CLE CLE CLE ];
        result=[];
        icel=0;
        
        %put very small numbers to zero (approximation issue)
        x(x<1e-6)=0;
        
        %add noc sec-act-tec data
        tmp=out_clemfr(out_clemfr(:,4)==nocID,:);
        global_data=[global_data out_clemfr(find(out_clemfr(:,4)~=nocID),24) ...
            out_clemfr(find(out_clemfr(:,4)~=nocID),25); tmp];
        
        xtmp=zeros(size(tmp,1),6);
        x_alpha=[x_alpha; xtmp];
        x_alphaCLE=[x_alphaCLE; xtmp];
        xtmp=zeros(size(tmp,1),1);
        x=[x; xtmp];
        CLE=[CLE; xtmp];
        
        %select sector-activity pairs to be considered for aggregated results
        tmp_sec_act=global_data(:,[2 3]);
        [sec_act , ~, c]=unique(tmp_sec_act,'rows');
        
        %select macrosectors to be considered for aggregated results
        tmp_ms=global_data(:,[1]);
        [ms , ~, c]=unique(tmp_ms,'rows');
        
        %unit costs
        uc=global_data(:,18);
        
        %loop over cells to create emissions file
        for i=1:length(flag_optim_dom)
            i;
            re=[];
            emi_alpha=[];
            
            %20120213 - EP
            %variables to save costs per cell
%             costOverCle_ms_cel_save(i,1:size(ms,1))=-999;
            %always write 11 macrosectors, in line 671
            costOverCle_ms_cel_save(i,1:11)=-999;
            
            %for cells inside the domain
            if (flag_optim_dom(i)==1 | flag_optim_dom(i)==2)
                icel=icel+1;
                
                %coordinate
                xutm=coordinate(i,1)*1000;
                yutm=coordinate(i,2)*1000;
                
                %load emissions
                emi_alpha=[];
                emi_alpha=emi{i,2};
                
                %remove last lines, related to emissions where no possible
                %reductions are appliable
                %emi_alpha_noc=emi_alpha(length(emi_alpha)-15:length(emi_alpha),:);
                %emi_alpha(length(emi_alpha)-15:length(emi_alpha),:)=[];
                emi_alpha_fin=emi_alpha(:,1:6);
                actlev=emi_alpha(:,7);
                %noc are sec-act where no control is available. Considers both
                %emi_alpha_noc and sa_noc
                %sa_noc=out_clemfr(length(emi_alpha)-15:length(emi_alpha),:);
                
                %both emissions and demoval efficiencies are ordered as: nox, voc,
                %nh3, pm10, pm25, so2
                re=global_data(:,6:11)/100; %removal efficiencies
                
                %emission reductions for single cells, detailed description
                emi_red_msat_cel_opt{i}=emi_alpha_fin.*re.*(x_alpha);
                emi_red_msat_cel_cle{i}=emi_alpha_fin.*re.*(x_alphaCLE);
                
                %compute cost for single cells, detailed description
                costTotal_msat_cel{i}=uc.*actlev.*(x);
                costCle_msat_cel{i}=uc.*actlev.*(CLE);
                costOverCle_msat_cel{i}=costTotal_msat_cel{i}-costCle_msat_cel{i};
                
                %data to be saved for detailed description
                data_msat_cel=global_data(:,1:5);
                
                %final results for single cells, detailed description
                %ms s a t low-high emired costcle costovercle
                emicosttab_msat_cel{i}=[data_msat_cel x emi_red_msat_cel_cle{i} emi_red_msat_cel_opt{i} ...
                    costTotal_msat_cel{i} costCle_msat_cel{i} costOverCle_msat_cel{i}];
                
                %results for single cells, aggregated for sec-act
                %do not consider difference between low and high
                for j=1:size(sec_act,1)
                    ind=find(global_data(:,2)==sec_act(j,1) & global_data(:,3)==sec_act(j,2));
                    ind1=find(global_data(:,2)==sec_act(j,1) & global_data(:,3)==sec_act(j,2) & global_data(:,5)==1);
                    ind2=find(global_data(:,2)==sec_act(j,1) & global_data(:,3)==sec_act(j,2) & global_data(:,5)==2);
                    
                    data_msa_cel(j,1:3)=[global_data(ind(1),1) sec_act(j,:)];
                    emi_red_msa_cel_opt{i}(j,1:6)=sum(emi_red_msat_cel_opt{i}(ind,:),1);
                    emi_red_msa_cel_cle{i}(j,1:6)=sum(emi_red_msat_cel_cle{i}(ind,:),1);
                    
                    costTotal_msa_cel{i}(j)=sum(costTotal_msat_cel{i}(ind));
                    costCle_msa_cel{i}(j)=sum(costCle_msat_cel{i}(ind));
                    costOverCle_msa_cel{i}(j)=costTotal_msa_cel{i}(j)-costCle_msa_cel{i}(j);
                    
                    %compute initial emissions, change if you take into account low or high level emissions
                    if (~isempty(ind1) && ~isempty(ind2))
                        emi_ini_msa_cel{i}(j,1:6)=emi_alpha_fin(ind1(1),:)+emi_alpha_fin(ind2(1),:);
                    elseif (~isempty(ind1) && isempty(ind2))
                        emi_ini_msa_cel{i}(j,1:6)=emi_alpha_fin(ind1(1),:);
                    elseif (isempty(ind1) && ~isempty(ind2))
                        emi_ini_msa_cel{i}(j,1:6)=emi_alpha_fin(ind2(1),:);
                    elseif (isempty(ind1) && isempty(ind2))
                        emi_ini_msa_cel{i}(j,1:6)=0;
                    end
                    
                    emi_rem_msa_cel{i}(j,1:6)=emi_ini_msa_cel{i}(j,1:6)-emi_red_msa_cel_opt{i}(j,1:6);
                    
                end
                
                %save results for single cells, aggregated for sec-act
                emicosttab_msa_cel{i}=[data_msa_cel emi_ini_msa_cel{i} emi_red_msa_cel_cle{i}...
                    emi_red_msa_cel_opt{i} emi_rem_msa_cel{i} costTotal_msa_cel{i}' costCle_msa_cel{i}' costOverCle_msa_cel{i}'];
                
                %results for single cells, aggregated for macrosector
                for j=1:size(ms,1)
                    ind=find(data_msa_cel(:,1)==ms(j));
                    
                    data_m_cel(j)=ms(j);
                    emi_red_m_cel_opt{i}(j,1:6)=sum(emi_red_msa_cel_opt{i}(ind,:),1);
                    emi_red_m_cel_cle{i}(j,1:6)=sum(emi_red_msa_cel_cle{i}(ind,:),1);
                    
                    emi_ini_m_cel{i}(j,1:6)=sum(emi_ini_msa_cel{i}(ind,:),1);
                    emi_rem_m_cel{i}(j,1:6)=emi_ini_m_cel{i}(j,1:6)-emi_red_m_cel_opt{i}(j,1:6);
                    
                    costTotal_m_cel{i}(j)=sum(costTotal_msa_cel{i}(ind));
                    costCle_m_cel{i}(j)=sum(costCle_msa_cel{i}(ind));
                    costOverCle_m_cel{i}(j)=costTotal_m_cel{i}(j)-costCle_m_cel{i}(j);
                end
                
                %save results for single cells, aggregated for macrosector
                emicosttab_m_cel{i}=[data_m_cel' emi_ini_m_cel{i} emi_red_m_cel_cle{i} ...
                    emi_red_m_cel_opt{i} emi_rem_m_cel{i} ...
                    costTotal_m_cel{i}' costCle_m_cel{i}' costOverCle_m_cel{i}'];
                
                %20120213 - EP
                %save cost over cle, per cell
                costOverCle_ms_cel_save(i,1:size(ms,1))=costOverCle_m_cel{i};
                
                %compute total over the optimization domain
                emicosttab_msat_tot=emicosttab_msat_tot+emicosttab_msat_cel{i}(:,7:21);
                emicosttab_msa_tot=emicosttab_msa_tot+emicosttab_msa_cel{i}(:,4:30);
                emicosttab_m_tot=emicosttab_m_tot+emicosttab_m_cel{i}(:,2:28);
                
                %sum part of emissions, for border cell, that cannot be
                %reduced
                %                 emi_ini_m_cel{i}(j,1:6)=sum(emi_ini_msa_cel{i}(ind,:),1)+sum(emi{i,1}); %sum also aggregated emissions for border cells
            end
        end
        
        emicosttab_msat_tot(abs(emicosttab_msat_tot)<1e-6)=0;
        
        %20120213 - EP
        %save costs per cell per macrosector
        % xutum yutm costms1 costms2 ... costms10 (cost over cle)
        nomefile=strcat(outdir,'/maps_cost/m','_point',int2str(index),'_costPerCellPerMs.csv');
        fid=fopen(nomefile,'wt');
        delim=',';
        fprintf(fid,'%s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s \n',...
            'xutm', delim, 'yutm', delim, 'CostOverCle_ms1 [Meuro/year]', delim, 'ms2',delim,...
            'ms3',delim,'ms4',delim,'ms5',...
            delim,'ms6',delim,'ms7',delim,'ms8',delim,...
            'ms9',delim,'ms10',delim,'ms11');
        fclose(fid);
        results=[coordinate costOverCle_ms_cel_save];
        dlmwrite(nomefile,results,'-append','roffset', 0, 'precision', 5);
        
        %data per msat, global data over domain
        re=global_data(:,6:11);
        arREF=global_data(:,25);
        arCLEPOT=global_data(:,19:20);
        %20130610 - save arOPT instead of cv; multiply arOPT x 100
        emicosttab_msat_tot=[data_msat_cel re arREF arCLEPOT x.*100 emicosttab_msat_tot];
        nomefile=strcat(outdir,'/emi_cost/msat/msat','_point',int2str(index),'.csv');
        fid=fopen(nomefile,'wt');
        delim=',';
        %20130610 - save arOPT instead of cv; multiply arOPT x 100
        fprintf(fid,'%s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s\n',...
            'ms', delim, 'sec', delim, 'act', delim, 'tec',delim,...
            'lowHigh',delim,...
            'reNox',delim,'reVoc',...
            delim,'reNH3',delim,'rePM10',delim,'rePM25',delim,...
            'reSo2',delim,'arREF',delim,'arCLE',delim,'arPOT',delim,...
            'arOPT',delim,'redCLENOx[ton]',delim,'redCLEVoc[ton]',delim,'redCLENh3[ton]',delim,...
            'redCLEPM10[ton]',delim,'redCLEPM25[ton]',delim,'redCLESO2[ton]',delim,'redwrtCLENox[ton]',delim,...
            'redwrtCLEVoc[ton]',delim,'redwrtCLENh3[ton]',delim,'redwrtCLEPM10[ton]',delim,...
            'redwrtCLEPM25[ton]',delim,'redwrtCLESO2[ton]',delim,'costTotal[Meuro]',...
            delim,'costCLE[Meuro]',delim,'costOverCLE[Meuro]');
        %'remNox',delim,'remVoc',...
        %delim,'remNH3',delim,'remPm10',delim,'remPM25',delim,...
        %'remSo2',delim,
        fclose(fid);
        dlmwrite(nomefile,emicosttab_msat_tot,'-append','roffset', 0, 'precision', 5);
        
        %data per m, global data over domain
        emicosttab_m_tot=[data_m_cel' emicosttab_m_tot];
        nomefile=strcat(outdir,'/emi_cost/m/m','_point',int2str(index),'.csv');
        fid=fopen(nomefile,'wt');
        delim=',';
        fprintf(fid,'%s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s %c %s\n',...
            'ms', delim,...
            'iniNox[ton]',delim,'iniVoc',delim,'iniNh3',delim,...
            'iniPM10',delim,'iniPM25',delim,'iniSO2',delim,...
            'redCLENox[ton]',delim,'redCLEVoc',delim,'redCLENh3',delim,...
            'redCLEPM10',delim,'redCLEPM25',delim,'redCLESO2',delim,...
            'redwrtCLENox[ton]',delim,'redwrtCLEVoc',delim,'redwrtCLENh3',delim,...
            'redwrtCLEPM10',delim,'redwrtCLEPM25',delim,'redwrtCLESO2',delim,...
            'remNox[ton]',delim,'remVoc',delim,'remNh3',delim,...
            'remPM10',delim,'remPM25',delim,'remSO2',delim,...
            'costTotal[Meuro]',delim,'costCLE[Meuro]',delim,'costOverCLE[Meuro]');
        fclose(fid);
        dlmwrite(nomefile,emicosttab_m_tot,'-append','roffset', 0, 'precision', 5);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FUNCTION TO COMPUTE EMISSIONS AND AQIS VALUES PER CELL

    function [emi_low,emi_high,aqi_val2]=POST_FG_compute_aqi(x,SR,DD,indexSea, aqiIndex, emi, commonDataInfo)
        
        %compute aqi related to an optimal solution
        %define low emissions matrix (to be used for emission computation)
        global_data=commonDataInfo.special.global_data;
        flag_region_dom=commonDataInfo.domainInfo.flag_region_dom;
        
        xtmplow=x(global_data(:,5)==1);
        xtmplow2=[xtmplow xtmplow xtmplow xtmplow xtmplow xtmplow];
        
        %define high emissions matrix (to be used for emission computation)
        xtmphigh=x(global_data(:,5)==2);
        xtmphigh2=[xtmphigh xtmphigh xtmphigh xtmphigh xtmphigh xtmphigh];
        
        %UNIBS(ET)20131001 - from here to line the end of the first for
        %flag_optim_dom has been changed with flag_region_dom
        
        %preallocate variables
        emi_low(length(flag_region_dom),6)=0;
        emi_high(length(flag_region_dom),6)=0;
        
        %loop over cells to create emissions file
        for i=1:length(flag_region_dom)
            
            if flag_region_dom(i)==0
                %case of cells without optimization domain
                %emission order: NOX, COV, NH3, PM10, PM25, SO2.
                emi_low(i,:)=flag_region_dom(i,:);
                emi_high(i,:)=flag_region_dom(i,:);
                
            elseif (flag_region_dom(i)==1 | flag_region_dom(i)==2)
                %for cells inside optimization domain
                %separate low and high emissions
                tmp_emi_low=emi{i,2}(global_data(:,5)==1,:);
                tmp_emi_high=emi{i,2}(global_data(:,5)==2,:);
                
                gdl=global_data(global_data(:,5)==1,:);
                gdl_sa=gdl(:,[2 3]);
                [a , ~, c]=unique(gdl_sa,'rows');
                for sa=1:size(a,1)
                    ind_sa_gdl=find(gdl(:,2)==a(sa,1) & gdl(:,3)==a(sa,2));
                    gdl_emibase=tmp_emi_low(ind_sa_gdl,:); %emissions for this cell, for all techs in sa
                    gdl_re=gdl(ind_sa_gdl,6:11)/100; %re for this cell, for all techs in sa
                    gdl_ar=xtmplow2(ind_sa_gdl,:); %ar for this cell, for all techs in sa
                    emiRed(sa,1:6)=sum(gdl_emibase(:,1:6).*gdl_re(:,:).*gdl_ar(:,:),1); %reduced emissions
                end
                emi_low(i,:)=base_emi_low(i,:)+base_emi_low_noc(i,:)-sum(emiRed);
                
                gdh=global_data(global_data(:,5)==2,:);
                gdh_sa=gdh(:,[2 3]);
                [a2 , ~, c]=unique(gdh_sa,'rows');
                emiRedHigh=[]; %20150304 ET - deal with DBs without punctual sources
                for sa=1:size(a2,1)
                    ind_sa_gdh=find(gdh(:,2)==a2(sa,1) & gdh(:,3)==a2(sa,2));
                    gdh_emibase=tmp_emi_high(ind_sa_gdh,:); %emissions for this cell, for all techs in sa
                    gdh_re=gdh(ind_sa_gdh,6:11)/100; %re for this cell, for all techs in sa
                    gdh_ar=xtmphigh2(ind_sa_gdh,:); %ar for this cell, for all techs in sa
                    emiRedHigh(sa,1:6)=sum(gdh_emibase(:,1:6).*gdh_re(:,:).*gdh_ar(:,:),1); %remaining emissions
                end
                emi_high(i,:)=base_emi_high(i,:)+base_emi_high_noc(i,:)-sum(emiRedHigh);
            end
        end
        
        %UNIBS(ET)20131001 - keep only cells in regional domain
        emi_low=emi_low(find(flag_region_dom==1 | flag_region_dom==2),:);
        emi_high=emi_high((flag_region_dom==1 | flag_region_dom==2),:);
        %emi_low=emi_low+emi_high;
        
        %in case of areal+point summed, D and d contain infos both or areal
        %and point, together
        %in case of areal and point separated, D and d are only related to
        %areal emissions (the point infos Dp and dp will be used later on)
        D=DD(aqiIndex).D;
        d=DD(aqiIndex).d;
        
        E = D * sparse(x);
        E = d - E;
        E_full = full(E);
        emis=E_full;
        
        %optimizerCondition=commonDataInfo.optimizerCondition;
        %optimizerValues=commonDataInfo.optimizerValues;

        precNames={'NOx';'NMVOC';'NH3';'PM10';'PM25';'SOx'};
        %execString=strcat('length(find(', optimizerCondition);
        %execString=strcat(execString, '))');
        %ncelopt=eval(execString);

        %restore first guess emissions
        s1_NOX = emis((ncelopt*0)+1:ncelopt*1);

        s1_VOC = emis((ncelopt*1)+1:ncelopt*2);

        s1_NH3 = emis((ncelopt*2)+1:ncelopt*3);

        s1_PM10 = emis((ncelopt*3)+1:ncelopt*4);

        s1_PM25 = emis((ncelopt*4)+1:ncelopt*5);

        s1_SO2 = emis((ncelopt*5)+1:ncelopt*6);
        
        %CASE OF AREAL+POINT, SUMMED
        %if areal_point==0
            
            %ET20140313
            %input selection based on ANNs features
            
            %read net and select precursors
            emissioni=[];
            for i_prec=1:size(precNames,1) %areal do always
                if strcmp(precNames(i_prec),'NH3')==1
                    emissioni=[emissioni, s1_NH3];
                elseif strcmp(precNames(i_prec),'NOx')==1
                    emissioni=[emissioni, s1_NOX];
                elseif strcmp(precNames(i_prec),'PM10')==1
                    emissioni=[emissioni, s1_PM10];
                elseif strcmp(precNames(i_prec),'PM25')==1
                    emissioni=[emissioni, s1_PM25];
                elseif strcmp(precNames(i_prec),'SOx')==1
                    emissioni=[emissioni, s1_SO2];
                elseif strcmp(precNames(i_prec),'NMVOC')==1
                    emissioni=[emissioni, s1_VOC];
                end
            end
        
        %input_rete2=emissioni';
        alpha=aggregationInfo.firstguess.alpha;
        omega=aggregationInfo.firstguess.omega;
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
        thisAlpha(commonDataInfo.flag_region_dom==0,:)=[];
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
        if ((aqiIndex == 1) || (aqiIndex == 2)) aqi_per_cell=sum(emissioni(:,[1 3 5 6]).*thisAlpha(:,[1 3 5 6]),2); end
        if (aqiIndex == 6) aqi_per_cell=sum(emissioni(:,[1]).*thisAlpha(:,[1]),2); end
        
        %change name to a better one!!! (too similar to caller...)
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint8_new','input_rete2');
        % 20160421 MM / EP SR / First guess Version
        %aqi_per_cell=interface_get_aqipercell(input_rete2, NN, aggregationInfo.mathIntermediateData, commonDataInfo, periodIndex, aqiIndex );
        periodIndex=1;
        aqi_val2=firstguess_get_aqipercell(aqi_per_cell, aggregationInfo, commonDataInfo, periodIndex, aqiIndex );
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FUNCTION TO COMPUTE EMISSIONS AND AQIS MAPS
    function []=POSTcomputeMaps(livelli,map,nomefile,aqilab,indexPareto,outdir,...
            pathPRB,pathPCI,flag_dom2,detailResults)
        %pathPRB,pathPCI)xutm,yutm,flag_dom2,detailResultsnx,ny2
        
        %create maps with optimization results, for aqi and emissions
        
        %map represents only the optimization domain cells, while flag_dom2
        %contains 0 outside optimization domain, and 1 inside optimization domain.
        %So I fill flag_dom2==1 with the values in map
        indrem=find(flag_dom2==0);
        flag_dom2(find(flag_dom2==1 | flag_dom2==2))=map;
        flag_dom2(indrem)=nan;
        aqi=[coordinate flag_dom2];
        
        %in case of detailed results, create maps
        if detailResults==1
            figure;
            %UNIBS(ET)20131002 - to transform coordinates are in m if they
            %are in Km
            if aqi(1,1)>1000
                unit=1;
            else
                unit=1000;
            end
            
            POSTcomputeMapsRead(livelli,aqi,unit,nomefile)
            hold on
            
            %draw regional boundaries
            POSTcomputeMapsBln(pathPRB);
            
            %draw cities
            POSTcomputeMapsPost(pathPCI);
            set(gcf,'Position',[1 31 1258 694]);
            set(gcf,'PaperPositionMode','auto');
            
            %save map of results
            file=strcat(outdir,'/',nomefile,'_point',int2str(indexPareto),'.png');
            print(gcf,'-dpng',file);
            close
        end
        
        %save ascii of results
        file=strcat(outdir,'/',nomefile,'_point',int2str(indexPareto),'.csv');
        
        fid=fopen(file,'wt');
        delim=',';
        
        fprintf(fid,'%s %c %s %c %s \n','xutm[km]', delim, 'yutm[km]', delim, cell2mat(aqilab));
        fclose(fid);
        
        aqi(isnan(aqi))=-999;
        dlmwrite(file,aqi,'-append','roffset', 0, 'precision', '%.6f');
        %         dlmwrite(file,aqi,'-append','roffset', 0);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FUNCTION ANCILLARY FOR MAPS CREATION
    function POSTcomputeMapsRead(livelli,dati,units,~)
        %contours of optimal aqi and emissions
        
        %create coordinates
        dx=dati(1,1);
        diversi=find(dati(:,1)-dx>0.001);
        %         diversi=find(dati(:,1)~=dx); %OLD RIAT
        dx=round(abs(dx-dati(diversi(1),1)));
        %         dx=abs(dx-dati(diversi(1),1));
        
        dy=dati(1,2);
        diversi=find(dati(:,2)-dy>0.001);
        %         diversi=find(dati(:,2)~=dy);
        dy=round(abs(dy-dati(diversi(1),2)));
        %         dy=abs(dy-dati(diversi(1),2));
        
        xorg=min(dati(:,1));
        yorg=min(dati(:,2));
        
        [nr nc]=size(dati);
        
        %create map
        for i=1:nr
            ii=round((dati(i,1)-xorg)/dx+1);
            jj=round((dati(i,2)-yorg)/dy+1);
            %             ii=(dati(i,1)-xorg)/dx+1; %OLD RIAT
            %             jj=(dati(i,2)-yorg)/dy+1;
            xcoord(ii)=(dati(i,1)+dx/2)*units;
            ycoord(jj)=(dati(i,2)+dy/2)*units;
            mappa(jj,ii)=dati(i,3);
        end
        
        %create contour
        contourf(xcoord,ycoord,mappa,livelli);
        set(gca,'fontsize',18)
        shading flat
        
        h=contourcmap(livelli,'jet','colorbar','on');
        set(h,'fontsize',18)
        
        y_tick=get(h,'ytick');
        y_tick=y_tick(1:2:length(y_tick));
        y_tick_label=get(h,'yticklabel');
        y_tick_label=y_tick_label(1:2:length(y_tick_label),:);
        set(h,'ytick',y_tick,'yticklabel',y_tick_label);
        
        minimo=min(dati(:,1:2));
        massimo=max(dati(:,1:2));
        
        set(gca,'XLim',[minimo(1)*units massimo(1)*units]);
        set(gca,'YLim',[minimo(2)*units massimo(2)*units]);
        
        asp_ratio=get(gca,'DataAspectRatio');
        asp_ratio(2)=asp_ratio(1);
        set(gca,'DataAspectRatio',asp_ratio);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FUNCTION ANCILLARY FOR MAPS CREATION
    function POSTcomputeMapsBln(infile)
        %depict domain and regional boundaries
        
        % read and plot bnl file
        hold on
        dati=load (infile);
        blocchi_dati=find(dati(:,2)==0);
        inizio=0;
        
        for i=1:length(blocchi_dati)
            dati_plot=dati(2+inizio:2+inizio+dati(blocchi_dati(i),1)-1,:);
            inizio=inizio+dati(blocchi_dati(i),1)+1;
            plot(dati_plot(:,1),dati_plot(:,2),'Color',[0.5 0.5 0.5])
        end
        
        dati(blocchi_dati,:)=[];
        
        minimo=min(dati);
        massimo=max(dati);
        
        asp_ratio=get(gca,'DataAspectRatio');
        asp_ratio(2)=asp_ratio(1);
        set(gca,'DataAspectRatio',asp_ratio);
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%FUNCTION ANCILLARY FOR MAPS CREATION
    function POSTcomputeMapsPost(infile)
        %depict cities
        
        %open file
        fid=fopen(infile);
        
        while (feof(fid)==0)
            riga=fgetl(fid);
            if (feof(fid)==0)
                coord_xy(1,:)=sscanf(riga,'%f');
                blank=regexpi(riga,'\s');
                nome_label{1}=strtrim(riga(blank(length(blank)):length(riga)));
                
                h=text(coord_xy(:,1),coord_xy(:,2),strcat('+ ',nome_label{1}));
                set(h,'Fontsize',12);
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [extCost_tot_morb,extCost_tot_yoll,resp_prob,card_prob,yoll_prob]=POSTexternalcosts(aqi_val,...
            pathextpop,pathextoth,pop,flag_optim_dom,flag_region_dom)
        %UNIBS(ET)20131001 - added  flag_region_dom in input
        
        % 	- "external_cost_pop_data"contains
        % 		% pop per age, in % (9 columns: 0-5 6-10 11-14 15-19 20-24 25-29 30-59 60-64 <65)
        % 		% mortality rate per age (10-18 columns: 0-5 6-10 11-14 15-19 20-24 25-29 30-59 60-64 <65)
        % 		% asmathics rate in % (19 column)
        % 		% fraction of population > 30 years (20 column)
        % 	- "external_cost_oth_data": contains impact and cost coefficients (2 columns)
        % 		ASMATICI
        % 		Adulti
        % 		 	1  uso di broncodilatatore
        % 			2  tosse
        % 			3  affaticamento nella respirazione (sintomi)
        % 		Bambini
        % 			4  uso di broncodilatatore
        % 			5  tosse
        % 			6  affaticamento nella respirazione (sintomi)
        % 		OLTRE I 65 ANNI
        % 			7  infarto
        % 		BAMBINI
        % 			8  tosse cronica
        % 		ADULTI
        % 			9  giorni di attività ridotta
        % 			10 bronchite cronica
        % 		POPOLAZ. TOT
        % 			11 mortalità cronica
        % 			12 ricoveri per problemi respiratori
        % 			13 ricoveri per problemi cardiocircolatori
        % 		SOPRA I 30 ANNI
        % 			14 anni di vita persi
        
        %reconstruct the full_aqi vector (with 0 outside optimization domain)
        full_aqi(1:length(flag_region_dom))=0;
        %         full_aqi(flag_region_dom==1)=aqi_val;
        full_aqi(flag_region_dom==1 | flag_region_dom==2)=aqi_val;
        full_aqi=full_aqi';
        
        %cut extpop data, if necessary (as in RIAT).
        extpop(length(flag_region_dom)+1:end,:)=[];
        
        %compute some interim values
        pop0_14=extpop(:,1)+extpop(:,2)+extpop(:,3);
        pop15_64=extpop(:,4)+extpop(:,5)+extpop(:,6)+extpop(:,7)+extpop(:,8);
        pop65_90=extpop(:,9);
        pop30_90=extpop(:,22);
        
        %asmatic 0-14
        pop_asm0_14=extpop(:,19);
        pop_nasm0_14=1-extpop(:,19);
        %asmatic-not asmatic 15-64
        pop_asm15_64=extpop(:,20);
        pop_nasm15_64=1-extpop(:,20);
        %asmatic-not asmatic 65-90
        pop_asm65_90=extpop(:,21);
        pop_nasm65_90=1-extpop(:,21);
        
        %compute impacts
        %         adults-asmatich: broncodilitatore, cough, respiratory problems
        impacts(:,1)=pop.*pop15_64.*pop_asm15_64.*full_aqi.*extoth(1,1);
        impacts(:,2)=pop.*pop15_64.*pop_asm15_64.*full_aqi.*extoth(2,1);
        impacts(:,3)=pop.*pop15_64.*pop_asm15_64.*full_aqi.*extoth(3,1);
        %         children-asmatich: broncodilitatore, cough, respiratory
        %         problems
        impacts(:,4)=pop.*pop0_14.*pop_asm0_14.*full_aqi.*extoth(4,1);
        impacts(:,5)=pop.*pop0_14.*pop_asm0_14.*full_aqi.*extoth(5,1);
        impacts(:,6)=pop.*pop0_14.*pop_asm0_14.*full_aqi.*extoth(6,1);
        %over 65-not asmatic: heart attack
        impacts(:,7)=pop.*pop65_90.*pop_nasm65_90.*full_aqi.*extoth(7,1);
        %children-not asmatic: chronic cough
        impacts(:,8)=pop.*pop0_14.*pop_nasm0_14.*full_aqi.*extoth(8,1);
        %adults-not asmatic: reduced activity days
        impacts(:,9)=pop.*pop15_64.*pop_nasm15_64.*full_aqi.*extoth(9,1);
        %adults-not asmatic: chronic cough
        impacts(:,10)=pop.*pop15_64.*pop_nasm15_64.*full_aqi.*extoth(10,1);
        %impact 11 is chronic death...to be computed later on
        
        %respiratory problems
        impacts(:,12)=pop.*full_aqi.*extoth(12,1);
        %define cells of optimization domain
        %UNIBS(ET)20131001 - create optfilter to filter external costs for PAD
        regdom=find(flag_region_dom==1 | flag_region_dom==2);
        optfilter=find(flag_optim_dom==1 | flag_optim_dom==2);
        %UNIBS(ET)20131001 - resp_prob data only for PAD
        resp_prob=impacts(optfilter,12);
        %problemi cardiovascolari
        impacts(:,13)=pop.*full_aqi.*extoth(13,1);
        %UNIBS(ET)20131001 - resp_prob data only for PAD
        card_prob=impacts(optfilter,13);
        %yoll
        impacts(:,14)=pop.*pop30_90.*full_aqi.*extoth(14,1);
        %         impacts(:,14)=full_aqi.*extoth(14,1)*12*60; %in monhts, over 60 years of life
        %         yoll_tmp=full_aqi.*extoth(14,1)*12; %over the 60 years of life, in months
        %20130213 - computing months of lost life (considering average
        %MOLL, you need to multiply the YOLL * 12, and also you need to
        %consider 60 years (to comulate yearly impact over the life of
        %people, from 30 to 90)
        yoll_prob=impacts(regdom,14)./(pop(regdom).*pop30_90(regdom))*12*60; %in monhts, over 60 years of life
        
        %computing chronic mortality - index 11
        for ageCat=1:9
            impacts_tmp(:,ageCat)=pop.*extpop(:,ageCat).*extpop(:,ageCat+9)...
                /100.*extoth(11,1)/100.*full_aqi;
        end
        impacts(:,11)=sum(impacts_tmp,2);
        
        %20130228 - separate yoll and morbidity
        %computing external costs
        for indImp=1:13
            costs_morb(:,indImp)=impacts(:,indImp).*extoth(indImp,2);
        end
        
        for indImp=14:14
            costs_yoll(:,indImp)=impacts(:,indImp).*extoth(indImp,2);
        end
        
        %UNIBS(ET)20131001 - ext costs data only for PAD
        costs_morb=costs_morb(optfilter,:);
        costs_yoll=costs_yoll(optfilter,14);
        
        extCost_tot_morb=sum(sum(costs_morb,2));
        extCost_tot_yoll=sum(sum(costs_yoll,2));
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [output_net_norm]=sim_exe(NN,input_rete_norm2)
        %20141121-ET new function to substitute matlab function sim
        sz=size(input_rete_norm2,2);
        %neural network simulation
        s_a=NN.net.IW{1,1}*input_rete_norm2+repmat(NN.net.b{1},1,sz);
        s_a1=eval(strcat(NN.net.layers{1}.transferFcn,'(s_a)'));
        s_b=NN.net.LW{2}*s_a1+repmat(NN.net.b{2},1,sz);
        output_net_norm=eval(strcat(NN.net.layers{2}.transferFcn,'(s_b)'));
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
