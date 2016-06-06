function [emi,flag_ADS_tp]=MAINload_emi(ncel,pathEMI,flag_optim_dom,flag_mode_ce_mo)
%load_emi loop over domain cells, create file_name to be loaded, and load
%in cell array the emissions, without checking if the cell is in or out the
%optimization domain (so depending on the cell, emi{i} contains 2 rows and
%6 columns (outside optim domain), or n rows and 6 columns (inside optim
%domain, n is quadruple numbers, both low and high)

%20120517 - now "emi" is a structure that contains in the first column the
%aggregated data, in the second column the detailed data (only for "border
%cells" the two dataset are loaded)

%20130424 - different data load depending on the type of use
%if MO, CE, detSA

%in this case 
% - emi{i,1} contains 6 columns (nox, voc, nh3, p10, p25, so2) and 2 rows (low and high emissions at CLE) 
% - emi{i,2} contains 6 columns (nox, voc, nh3, p10, p25, so2) and as many rows as many available technologies
%flag_optim_dom=commonDataInfo.domainInfo.flag_optim_dom;
%ncelopt=length(find(flag_optim_dom==1 | flag_optim_dom==2));

if (flag_mode_ce_mo==0 | flag_mode_ce_mo==1 | flag_mode_ce_mo==2)
    %emi=cell(ncel,2);
    for i=1:ncel
        if flag_optim_dom(i)==0 %cells outside optimization cell
            ci=int2str(i);
            %      tmp1=sprintf('%4s',ci);
            fileemi1=strcat(pathEMI,'aggregated/',ci,'.txt');
            %     fileemi=strcat('./input/emissions/',tmp1,'.txt');
            emi{i,1}=load(fileemi1);
        elseif flag_optim_dom(i)==1
            ci=int2str(i);
            fileemi2=strcat(pathEMI,'detailed/',ci,'.txt');
            emi{i,2}=load(fileemi2);
        elseif flag_optim_dom(i)==2
            ci=int2str(i);
            fileemi1=strcat(pathEMI,'aggregated/',ci,'.txt');
            emi{i,1}=load(fileemi1);
            fileemi2=strcat(pathEMI,'detailed/',ci,'.txt');
            emi{i,2}=load(fileemi2);
        end
    end
    flag_ADS_tp=[];
%aggregated SA
%in this case 
% - emi{i,1} contains 12 columns (6 nox, voc, nh3, p10, p25, so2 areal emissions, and then point emissions)...
%   and 11 row (emissions at CLE, per macrosector) 
% - emi{i,2} = 0
elseif flag_mode_ce_mo==3
    
    %periods yea, win, sum
    period={'TP1/','TP2/','TP3/'};
    
    %UNIBS(ET)20130927 
    %Check the existence of the emission files for each period
    %keep track of existent-nonExistent files in a flag
    flag_ADS_tp(1)=exist(strcat(pathEMI,period{1},'areal/NOx.txt'), 'file');
    flag_ADS_tp(2)=exist(strcat(pathEMI,period{2},'areal/NOx.txt'), 'file');
    flag_ADS_tp(3)=exist(strcat(pathEMI,period{3},'areal/NOx.txt'), 'file');
    
    for per=1:3
        
        %UNIBS(ET)20130927
        %if the files fot that period don'exist use existent files to fill the matrix
        if flag_ADS_tp(per)==2
        per2=per;
        else
          ExistEmi=find(flag_ADS_tp==2);
          per2=ExistEmi(1);  
        end
        
        %AREAL
        nox{per}=importdata(strcat(pathEMI,period{per2},'areal/NOx.txt'));
        voc{per}=importdata(strcat(pathEMI,period{per2},'areal/VOC.txt'));
        nh3{per}=importdata(strcat(pathEMI,period{per2},'areal/NH3.txt'));
        p10{per}=importdata(strcat(pathEMI,period{per2},'areal/PM10.txt'));
        p25{per}=importdata(strcat(pathEMI,period{per2},'areal/PM25.txt'));
        so2{per}=importdata(strcat(pathEMI,period{per2},'areal/SO2.txt'));
        
        %POINT
        noxp{per}=importdata(strcat(pathEMI,period{per2},'point/NOx.txt'));
        vocp{per}=importdata(strcat(pathEMI,period{per2},'point/VOC.txt'));
        nh3p{per}=importdata(strcat(pathEMI,period{per2},'point/NH3.txt'));
        p10p{per}=importdata(strcat(pathEMI,period{per2},'point/PM10.txt'));
        p25p{per}=importdata(strcat(pathEMI,period{per2},'point/PM25.txt'));
        so2p{per}=importdata(strcat(pathEMI,period{per2},'point/SO2.txt'));
    end
    
    %data structure - yea (nox voc nh3 p10 p25 so2, at first areal and then
    %point), win (nox voc nh3 p10 p25 so2, at first areal and then
    %point), sum(nox voc nh3 p10 p25 so2, at first areal and then
    %point)
    %
    for i=1:ncel
        emi{i,1}=[nox{1}.data(i,3:end)' voc{1}.data(i,3:end)' nh3{1}.data(i,3:end)' p10{1}.data(i,3:end)' p25{1}.data(i,3:end)' so2{1}.data(i,3:end)'...
                  noxp{1}.data(i,3:end)' vocp{1}.data(i,3:end)' nh3p{1}.data(i,3:end)' p10p{1}.data(i,3:end)' p25p{1}.data(i,3:end)' so2p{1}.data(i,3:end)'...
                  nox{2}.data(i,3:end)' voc{2}.data(i,3:end)' nh3{2}.data(i,3:end)' p10{2}.data(i,3:end)' p25{2}.data(i,3:end)' so2{2}.data(i,3:end)'...
                  noxp{2}.data(i,3:end)' vocp{2}.data(i,3:end)' nh3p{2}.data(i,3:end)' p10p{2}.data(i,3:end)' p25p{2}.data(i,3:end)' so2p{2}.data(i,3:end)'...
                  nox{3}.data(i,3:end)' voc{3}.data(i,3:end)' nh3{3}.data(i,3:end)' p10{3}.data(i,3:end)' p25{3}.data(i,3:end)' so2{3}.data(i,3:end)'...
                  noxp{3}.data(i,3:end)' vocp{3}.data(i,3:end)' nh3p{3}.data(i,3:end)' p10p{3}.data(i,3:end)' p25p{3}.data(i,3:end)' so2p{3}.data(i,3:end)'];
        emi{i,2}=0;
    end
end