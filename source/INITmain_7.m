%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%20130830-EP
% - added the words "PROGRESSION: " in the beginning of command window messages

%20120517 - updates in v1.4.2
% INIT is modified to manage the "border cells". Now the tool read emissions
% emissions directory (usually "./input/emissions/") considering two subdirectories:
% 1) detailed (for emissions detailed as sec-act-tec)
% 2) aggregated (for emissions aggregated for macrosector)
% in case the cells are partly inside and partly outside the optimization
% domain, the emissions are read in both the directories

%read emissions are saved in emi{:,1} (aggregated emissions), and emi{:,2}
%(detailed emissions). For border cells you have both columns of data

%modified routines: MAINload_emi, MAINbase_emi (if cell is completely or
%partly in optimization domainm you sum up detailed and, if available,
%aggregated data); MAINmain and POSTmain now consider that optimization
%domain is defined by flag_region_dom==1 and flag_region_dom==2

%for border cells, in "MAINbase_emi" the contribution of the part of the
%cell that is outside the optimization domain, is added to the base_emi_low
%and base_emi_low

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NB: the INIT is able to create new D and d matrices, processing the emission
%DB (first line read in file flag_optim_path, using icells that is read from ANNs files
% (line 10-15 file flag_optim_path), and creating the
% D,d file using names from first 6 lines of Dd matrices (line 28-33
% flag_optim_path)
% Even if you process a different database (not yearly, but seasonal i.e.,
% you need to use these lines to define ANNs paths and D,d filenames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This version of the code is running both on linux and matlab
%function []=INITmain_(f1, aggregationInfoFile, prepareDataInfoFile, finalizeJobInfoFile)
function []=INITmain_7(f1, f2)

%to run the routine, write
% INITmain('inputOPERAaspa/flag_optim_path.txt');

%20111206 generalization of quadrants: loading the ANNs, the tool loads the
%dimension of quadrants (icells) contained in the ANNs
%So D and d are computed, depending on the dimension of quadrants (icells)
%, that can be different for each considered pollutant

%load configuration files
% f1='inputOPERAaspa/flag_optim_path.txt';
% f2='input/POSTflag_optim_path.txt';

%to keep track of computing time
tic

%remove warning
%needed otherwise the compile toolbox is not working
warning off MCR:ClassInUseAtExit
warning off MATLAB:ClassInstanceExists
warning('off','all')
%read file input/flag_optim_path


%new MM

%commonDataInfo=load_CommonDataInfo(f1, 'init');
aggregationInfo=load_AggregationInfo(f2);
commonDataInfo=load_CommonDataInfo(f1, aggregationInfo.type);

%new MM end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       aqi_definition.txt

% preparing the data structure to store the AQIs configuration
% 7 columns are the 7 AQIs
% 3 structures are: 1)annual, 2)winter, 3)summer

%commonDataInfo=interface_prepareCommonData(aggregationInfo, commonDataInfo, f1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       flag_optim_oth.txt
%                  (optimization configuration)
%output file name
% nomeOUTPUT=pathOUF;

%create directory name
mkdir('output')
%load gains DB
%load file "out_clemfr.txt", containing informations all quadrupes corinair
%macrosector-sectorIiasa-activityIiasa-tech, and related informations
% global_data=load('input/out_clemfr.txt');
global_data=load(commonDataInfo.dirs.pathOCM);

%extract part of gains DB related to emissions for which no reduction
%technologies are available
global_data_noc=global_data(find(global_data(:,4)==commonDataInfo.info.nocID),:);
global_data(find(global_data(:,4)==commonDataInfo.info.nocID),:)=[];
%global_data_noc1=global_data_noc;
%save('C:\data\work\projects\riat\RiatPlus-v3beta\global_data_noc1', 'global_data_noc1');

%extract arADSdet and arREF variable...
% 20130402 - arADSdet is a variable used for the detailed scenario mode...it
% has to contain AR values between 0 and 1
arADSdet=global_data(:,24);
arREF=global_data(:,25);
global_data(:,25)=[];
global_data(:,24)=[];

%COLUMNS OF OUT_CLEMFR FILE
%macrosector sector activity tecnology lowHighEmissionFlag
%removalEfficiencies (NOx,VOC,NH3,PM10,PM25,SO2)
%unabatedEmissionFactores (NOx,VOC,NH3,PM10,PM25,SO2)
%unitcost arCLE_year arPOT_year flagReplace flagOptim flagNTM_TM arADSdet
%atREFyear

%COLUMNS OF OUT_CLEMFR FILE
%macrosector sector activity tecnology lowHighEmissionFlag
%removalEfficiencies (NOx,VOC,NH3,PM10,PM25,SO2)
%unabatedEmissionFactores (NOx,VOC,NH3,PM10,PM25,SO2)
%unitcost arCLE_year arPOT_year flagReplace flagOptim ar2005

%hl_flag=1 (low), =2(high);
%replace=1 (not replaceable)
%optim=1 (to be used in the optimization)

%NOTE: AR and RE are in %
[commonDataInfo]=init_Domain(commonDataInfo);
% commonDataInfo.dirs.pathFOD
%FROM HERE INSERT NEW LOOP
%new pathEMI updated
pathEMIyea=strcat(commonDataInfo.dirs.pathEMI,'TP1/');
pathEMIwin=strcat(commonDataInfo.dirs.pathEMI,'TP2/');
pathEMIsum=strcat(commonDataInfo.dirs.pathEMI,'TP3/');

pathVEC={pathEMIyea,pathEMIwin,pathEMIsum};

%if summer and winter directory do not exist, only create yearly matrices
%of precomputed emissions
if (exist(pathVEC{2},'dir')==0 & (exist(pathVEC{3},'dir')==0))
    indini=1;
    indfini=1;
    %if both summer and winter directory exist, create summer and winter
    %precomputed emissions preprocessing data, and yearl preprocessed
    %emissions as the sum of winter and summer
elseif (exist(pathVEC{2},'dir')==7 & (exist(pathVEC{3},'dir')==7))
    indini=1;
    indfini=3;
else
    error('No data available for emission preprocessing')
end
coordinate=commonDataInfo.domainData.data(:,1:2);
commonDataInfo.coordinate=coordinate;
for k=indini:indfini %year, winter, summmer
    
    emi={};
    pathEMI=pathVEC{k};
    
    if size(commonDataInfo.domainData.data,2)==6
        flag_region_dom=commonDataInfo.domainData.data(:,3); %regional domain
        flag_optim_dom=commonDataInfo.domainData.data(:,4);%optimization domain
        flag_aqi_dom=commonDataInfo.domainData.data(:,5); %aqi computation domain
        pop=commonDataInfo.domainData.data(:,6); %population
    else
        flag_region_dom=commonDataInfo.domainData.data(:,3);
        flag_optim_dom=commonDataInfo.domainData.data(:,3);%optimization domain
        flag_aqi_dom=commonDataInfo.domainData.data(:,4); %aqi computation domain
        
        pop=commonDataInfo.domainData.data(:,5); %population
    end
    %     stepsize=max(coordinate(2,1)-coordinate(1,1),coordinate(2,2)-coordinate(1,2));
    %
    if isequal(aggregationInfo.type, 'FIRSTGUESS')
        [ncel, nx, ny]=calcCellNo('latlon', commonDataInfo);
    else
        [ncel, nx, ny]=calcCellNo('utm', commonDataInfo);
    end
    
    %load cells emissions
    %emi: emissions of the starting case (2 rows in the case of cells
    %outside domain, sector-activity-techs rows if inside domain)
    %if the cell is inside the domain, the correspondent file contains both
    %emissions (first 6 columns) and activity level (seventh column)
    %emi contains two columnes: column 1 is the aggregated data, column 2 the
    %detailed data
    [emi,~]=MAINload_emi(ncel,pathEMI,flag_optim_dom,commonDataInfo.optim_flags.mode_ce_mo);
    
    %load PM10 to PM25 relationship
    %     pm10aveToExceed=load(pathPM);
    
    %in case of RIAT, matrices and input are cut and flipped (due to an initial
    %problem in the data). Otherwise no change to the data is performed
    %newEmi=emi;
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\newEmi', 'newEmi');
    if commonDataInfo.optim_flags.coordinateFlip==1
        [newDomainData, newDomainInfo, newEmi]=gridCorrectionCase1(domainData,domainInfo);
        domainData=newDomainData;
        domainInfo=newDomainInfo;
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\newEmi', 'newEmi');
        %emi=newEmi;
    end
    
    %imagesc(flipud(reshape(flag_optim_dom,ny,nx))) % spatial map
    if commonDataInfo.optim_flags.coordinateFlip==2 %ASPA CASE
        [newDomainData, newDomainInfo, newEmi]=gridCorrectionCase2(domainData,domainInfo);
        domainData=newDomainData;
        domainInfo=newDomainInfo;
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\newEmi', 'newEmi');
        emi=newEmi;
    end
    
    %define number of optimization cells: to be used to apply D and d in the
    %routine MAINcompute_aqi (1 is for cells inside optim domain, 2 for border
    %cells)
    %ncelopt=length(find(flag_optim_dom==1 | flag_optim_dom==2));
    
    %load or compute base emissions for each cell sum of emissions, for each cell,
    %separated for low and high rows extracte from emi (one row per each sec-act) and then summed
    %emissions of the part of the border cell outside optimziation domain, is
    %summed up to base_emi_low and base_emi_high
    %actlev_final=activity level of all sector-activity-techs
    %actlev_final_sum=optimization domain wide sum of activity levels
    
    %global_data1=global_data;
    %save('C:\data\work\projects\riat\RiatPlus-v3beta\global_data1', 'global_data1')
    [actlev_final,actlev_final_sum,base_emi_low,base_emi_high,base_emi_low_noc,...
        base_emi_high_noc]=MAINbase_emi(emi,global_data,flag_optim_dom,flag_region_dom,global_data_noc);
    
    %create matrices D and d, for quadrant computation (matrix are
    %4times greater than the starting one
    %in case of areal+point emissions, only D and d are created
    %in case ot areal and point emissions separated, D, d, and also Dp and dp are created
    %in this second cass, the matrices contain only the relative areal or point
    %contribution (D and d areal, Dp and dp point)
    %CREATE D and d only if ANN is available for the considered AQI
    for indaqi = 1:size(commonDataInfo.pathANN(1).ANNs,1), %aqis
        if (isequal(strtrim(commonDataInfo.pathANN(k).ANNs(indaqi,:)),'-999')==0) %check if AQI is defined
            if exist(strcat((strtrim(commonDataInfo.pathDd(k).Dd(indaqi,:))),'.mat'), 'file')==0 %check if Dd is already existing
                
                disp(strcat('PROGRESSION: Creating preprocessed emissions, for AQI -->',strtrim(commonDataInfo.pathANN(k).ANNs(indaqi,:))))
                %[D,d,Dp,dp]=INITaggregation(emi,global_data,...
                %    base_emi_low, base_emi_high, base_emi_low_noc,...
                %    base_emi_high_noc,flag_optim_dom,nx,ny,strtrim(pathANN(k).ANNs(indaqi,:)),optim_flags.areal_point);
                %load('C:\data\work\projects\riat\RiatPlus-v6.1\inputVars');
                %, 'emi1', ...
                %'global_data1', ...
                %'base_emi_low1', ...
                %'base_emi_high1', ...
                %'base_emi_low_noc1', ...
                %'base_emi_high_noc1', ...
                %'flag_region_dom1', ...
                %'nx1', ...
                %'ny1', 'areal_point1');
                
                %[D,d,Dp,dp]=INITaggregation(emi,global_data,...
                %    base_emi_low, base_emi_high, base_emi_low_noc,...
                %    base_emi_high_noc,flag_optim_dom,nx,ny,strtrim(pathANN(k).ANNs(indaqi,:)),optim_flags.areal_point,...
                %    aggregationInfo);
                %aggregationInfo=interface_setAggregationInfo(aggregationInfo);
                commonDataInfo=interface_new_fill_common_data_info(commonDataInfo, k, indaqi, aggregationInfo);
                
                %                 if isequal(aggregationInfo.type, 'FIRSTGUESS')
                %                     commonDataInfo=firstguess_fill_commonDataInfo(commonDataInfo, k, indaqi);
                %                 else
                %                     commonDataInfo=quadrant_fill_commonDataInfo(commonDataInfo, k, indaqi);
                %                 end
                
                [D,d,Dp,dp]=INITaggregation(emi,global_data,...
                    base_emi_low, base_emi_high, base_emi_low_noc,...
                    base_emi_high_noc,flag_region_dom,nx,ny,commonDataInfo.optim_flags.areal_point,...
                    commonDataInfo, aggregationInfo, k, indaqi);
                %base_emi_high_noc,flag_region_dom,nx,ny, strtrim(pathANN(k).ANNs(indaqi,:)),optim_flags.areal_point,...
                save(strtrim(commonDataInfo.pathDd(k).Dd(indaqi,:)),'D','d','Dp','dp');
            else
                disp(strcat('PROGRESSION: Preprocessed emissions already created, for AQI -->',strtrim(commonDataInfo.pathANN(k).ANNs(indaqi,:))))
            end
        end
    end
end

toc
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DDsparse,ddsparse,DDsparsep,ddsparsep]=INITaggregation(emi,global_data,...
    base_emi_low, base_emi_high, base_emi_low_noc,...
    base_emi_high_noc,flag_region_dom,nx,ny,areal_point, ...
    commonDataInfo, aggregationInfo, index1, index2)

%in case of areal+point emissions, only DDsparse and ddsparse contain
%useful informations (the DDsparsep and ddsparsep are empty matrices)

%routine to create matrices D and d, to be used to compute
%quadrant emissions, applying:
% d - D * x

DDsparsep=[];
ddsparsep=[];
DDsparse=[];
ddsparse=[];
%to define icells, load pathANNpm10
%#function network
%[net]=net_read(pathANNpm10); %20140617 ET - read model from txt file

%icells=net.icells;

%defininf variables
np=size(base_emi_low,2);    %number of pollutants
nc=length(flag_region_dom);  %number of cells
ncv=size(global_data,1);    %number of control variables
re=global_data(:,6:11)/100; %removal efficiencies

% 20160418: MM
% set to 1 with SR/Sherpa/First Guess mode

quad=interface_new_get_quadrant(aggregationInfo);
% if isequal(aggregationInfo.type, 'FIRSTGUESS')
%     quad =1;
% else
%     % set to 4 with Quadrant
%     quad=4;
% end

%info to cut ddsparse and ddsparsep, considering only optimization cells
%this clipping allows to speed up the procedure

indopt=flag_region_dom; %cells in optimization domain
indoptrep=repmat(indopt,quad*np,1); %index repeated for pollutants and quadrants
indiciMAT=reshape(indopt,ny,nx);

%initialize matrices
emi_red(1:nc,1:np,1:ncv)=0;
emi_rem(1:nc,1:np)=0;

%loop over cells to create emission matrices

if isequal(aggregationInfo.type, 'FIRSTGUESS')
    % load
    areas=importdata(commonDataInfo.dirs.pathArea);
    areaIn=areas.data(:,5);
    areaOut=areas.data(:,6);
    %put 1 in the two variables, to avoid division per zero
    %     areaIn(areaIn==0)=1;
    %     areaOut(areaOut==0)=1;
else
    areaIn=repmat(1,1,length(flag_region_dom))';%repmat(1,length(flag_region_dom));
    areaOut=repmat(1,1,length(flag_region_dom))';%repmat(1,length(flag_region_dom));
end
for i=1:length(flag_region_dom)
    
    %case of cells without optimization domain
    %emission order: NOX, COV, NH3, PM10, PM25, SO2.
    if flag_region_dom(i)==0
        % remember to cancel...
        %         i
        %INITIAL EMISSIONS
        %case areal and point summed:
        if areal_point==0
            %old version
            %emi_rem(i,:)=base_emi_low(i,:)+base_emi_high(i,:);
            % new Version
            nArea=repmat(areaOut(i), 1, np)+repmat(areaIn(i), 1, np);
            if sum(nArea)==0
                nArea=repmat(1,1,np);
            end
            emi_rem(i,:)=(base_emi_low(i,:)+base_emi_high(i,:))./nArea;
            % only for FIRST GUESS
            %case areal and point separated
        elseif areal_point==1
            emi_rem_low(i,:)=base_emi_low(i,:);
            emi_rem_high(i,:)=base_emi_high(i,:);
        end
        
    else     %for cells inside optimization domain
        
        %EMISSIONS REDUCED (BOTH "AREAL+POINT" AND "AREAL AND POINT")
        %define emission matrix
        emi_red_tmp{i}=emi{i,2}(1:ncv,1:6);
        %compute emissions*RE
        % working version
        %emi_red(i,1:np,1:ncv)=(emi_red_tmp{i}.*re)';
        % new version
        nAreaIn=repmat(areaIn(i), ncv, np);
        emi_red(i,1:np,1:ncv)=((emi_red_tmp{i}.*re)./nAreaIn)';
        
        %INITIAL EMISSIONS
        %compute final component of D matrix
        if areal_point==0
            % working version
            %emi_rem(i,:)=base_emi_low(i,:)+base_emi_high(i,:)+...
            %    base_emi_low_noc(i,:)+base_emi_high_noc(i,:);
            % new version
            nArea=repmat(areaIn(i), 1, np)+repmat(areaOut(i), 1, np);
            emi_rem(i,:)=(base_emi_low(i,:)+base_emi_high(i,:)+...
                base_emi_low_noc(i,:)+base_emi_high_noc(i,:))./nArea;
            %case areal and point separated
        elseif areal_point==1
            emi_rem_low(i,:)=base_emi_low(i,:)+...
                base_emi_low_noc(i,:);
            emi_rem_high(i,:)=base_emi_high(i,:)+...
                base_emi_high_noc(i,:);
        end
    end
end

%COMPUTE d: INITIAL EMISSIONS
%grid emissions are as NOX, COV, NH3, PM10, PM25, SO2
%areal+point emissions
%optimizerCondition='(optimizerValues==1 | optimizerValues==2)';
%unused but to fill with a value
indicators=0;
x=0;
y=0;

if areal_point==0
    
    % 20160418 MM Quadrant case
    % res=do_FullJob(aggregationInfo.geometryDataInfo, commonDataInfo, emi_rem, orderedPolls, indicators, x, y, nx, ny, indiciMAT, indoptrep, 1, 1);
    % orderedPolls = {'NOX';'VOC';'NH3';'PM10';'PM25';'SO2'};
    % 20160418 MM First Guess case
    
    orderedPolls=interface_new_get_pollutant_list(commonDataInfo);
    
    [jobIntermediate jobRes]=do_FullJob(aggregationInfo.geometryDataInfo, commonDataInfo, emi_rem, orderedPolls, indicators, x, y, nx, ny, indiciMAT, indoptrep, 1, 1);
    ddsparse=jobRes.finalGrid;
    % END
    clearvars dd
    %areal and point emissions separated
elseif areal_point==1
    
    %areal
    orderedPolls={'NOX_d';'VOC_d';'NH3_d';'PM10_d';'PM25_d';'SO2_d'};
    [jobIntermediate, jobRes]=do_FullJob(aggregationInfo.geometryDataInfo, commonDataInfo, emi_rem_low, orderedPolls, indicators, x, y, nx, ny, indiciMAT, indoptrep, 1, 1);
    ddsparse=jobRes.finalGrid;
    %point
    orderedPolls = {'NOX_dp';'VOC_dp';'NH3_dp';'PM10_dp';'PM25_dp';'SO2_dp'};
    [jobIntermediate, jobRes]=do_FullJob(aggregationInfo.geometryDataInfo, commonDataInfo, emi_rem_high, orderedPolls, indicators, x, y, nx, ny, indiciMAT, indoptrep, 1, 1);
    ddsparsep=jobRes.finalGrid;
    
    %keep optimization cells (inside optimization==1 and border cells==2)
    clearvars dd ddp
end


%NOW START TO CONSIDER EMISSIONS THAT ARE REMOVED DUE TO APPLICATION OF TECHNOLOGIES

%define D matrix as empty matrix (only considering PAD cells)
%indoptrep%%indopt repmat(indopt,quad*np,1); %index repeated for pollutants and quadrants
%indoptrep=repmat(indopt,intermediateResult.geomFactor*np,1); %index repeated for pollutants and quadrants

DD=zeros(size(find(indoptrep==1 | indoptrep==2),1),size(global_data,1));
if areal_point==1
    DDp=zeros(size(find(indoptrep==1 | indoptrep==2),1),size(global_data,1));
end

%loop to create matrix D
%areal+point emissions
if areal_point==0
    orderedPolls = {'NOX_D';'VOC_D';'NH3_D';'PM10_D';'PM25_D';'SO2_D'};
    for i=1:size(global_data,1)
        %new version        
        [jobIntermediate, jobRes]=do_FullJob(aggregationInfo.geometryDataInfo, commonDataInfo, emi_red(:,:,i), orderedPolls, indicators, x, y, nx, ny, indiciMAT, indoptrep, 0, 1);
        tmp=jobRes.finalGrid;
        
        DD(:,i)=tmp;%((indoptrep==1 | indoptrep==2));
    end
    DDsparse=sparse(DD);
    clearvars DD
    
    %areal and point emissions separated
elseif areal_point==1
    
    %     DDp=zeros(size(ddp,1),size(global_data,1));
    for i=1:size(global_data,1)
        %disp(strcat('PROGRESSION: step', num2str(i), '/',num2str(size(global_data,1))));
        orderedPolls = {'NOX_D';'VOC_D';'NH3_D';'PM10_D';'PM25_D';'SO2_D'};
        [jobIntermediate, jobRes]=do_FullJob(aggregationInfo.geometryDataInfo, commonDataInfo, emi_red(:,:,i), orderedPolls, indicators, x, y, nx, ny, indiciMAT, indoptrep, 0, 1);
        temp=jobRes.finalGrid;
        
        if global_data(i,5)==1 %areal emissions
            DD(:,i)=temp;
            
        elseif global_data(i,5)==2 %point emissions
            DDp(:,i)=temp;
        end
        
    end
    %create sparse matrix
    DDsparse=sparse(DD);
    DDsparsep=sparse(DDp);
end

%cut ddsparse and ddsparsep
disp('PROGRESSION: matrices d and D created');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%