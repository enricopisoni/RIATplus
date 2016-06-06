function [actlev_final,actlev_final_sum,base_emi_low,base_emi_high,...
    base_emi_low_noc,base_emi_high_noc]=MAINbase_emi(emi,global_data,...
    flag_optim_dom,flag_region_dom,global_data_noc)
%compute base emissions
%UNIBS(ET)20131001 - added flag_region_dom in input
%for border cells, in "MAINbase_emi" the contribution of the part of the
%cell that is outside the optimization domain, is added to the base_emi_low
%and base_emi_low

%emi{:,1} are emissions aggregated for macrosector (for cells completely or 
%partly outside optimization domain) 
%emi{:,2} are detailed emissions (completely or partly inside optimization domain) 

%when you compute "base_emi", you sum up aggregated emissions (completely outside
%optimization domain), detailed emission (inside optimization domain), 
%aggregated and detailed emissions (for border cells)

%or aggregated emissions, or part of emissions for border cells)

%looking at columns 2 and 3 (sector-activity), and extract rows
%representing each s-a (so remove double lines in the tmp_sec_act_### file)
ind1=find(global_data(:,5)==1);
tmp_sec_act_low=global_data(ind1,[2 3]);
[a blow c]=unique(tmp_sec_act_low,'rows');

%looking at columns 2 and 3 (sector-activity), and extract rows
%representing each s-a (so remove double lines in the tmp_sec_act_### file)
ind2=find(global_data(:,5)==2);
tmp_sec_act_high=global_data(ind2,[2 3]);
[a bhigh c]=unique(tmp_sec_act_high,'rows');

actlev_final_sum=0;

%loop over cells to create emissions file
for i=1:length(emi)
    

    %UNIBS(ET)20131001 - if cell outside regional domain,
                       % only rewrite file content
    if flag_region_dom(i)==0
        base_emi_low(i,:)=emi{i,1}(1,:);
        base_emi_high(i,:)=emi{i,1}(2,:);
        
        % 20160418: First Guess
        %base_emi_low_noc(i,1:5)=zeros(1,5);
        %ase_emi_high_noc(i,1:5)=zeros(1,5);
        % 20160418: Quadrant
        base_emi_low_noc(i,1:6)=zeros(1,6);
        base_emi_high_noc(i,1:6)=zeros(1,6);
        % END

        actlev_final{i}=0;
        
            %UNIBS(ET)20131001 - inside regional domain, or border cells
    elseif (flag_region_dom(i)==1 | flag_region_dom(i)==2)
        %extract emissions where no technologies can be applied
        emi_noc=emi{i,2}(size(global_data,1)+1:end,:);
        emi{i,2}(size(global_data,1)+1:end,:)=[];
        
        %emission that cannot be reduced %20120722 Nicola found bug
        % 20160418: First Guess
        %base_emi_low_noc(i,1:5)=sum(emi_noc(find(global_data_noc(:,5)==1),1:5),1);
        %base_emi_high_noc(i,1:5)=sum(emi_noc(find(global_data_noc(:,5)==2),1:5),1);
        % 20160418: Quadrant
        base_emi_low_noc(i,1:6)=sum(emi_noc(find(global_data_noc(:,5)==1),1:6),1);
        base_emi_high_noc(i,1:6)=sum(emi_noc(find(global_data_noc(:,5)==2),1:6),1);
        % END
        
        
        %separate low and high emissions
        tmp_emi_low=emi{i,2}(ind1,1:6);
        tmp_emi_high=emi{i,2}(ind2,1:6);
        
        %low emissions treatment
        %base emissions: select rows to get list of sector-activity emis
        tmp_base_emi_low=tmp_emi_low(sort(blow),:);
        %order of emissions NOX, COV, NH3, PM10, PM25, SO2
        if isempty(emi{i,1})==0
            base_emi_low(i,:)=sum(tmp_base_emi_low,1)+emi{i,1}(1,:); %border cell: sum emi cell part outside optimization domain;
        else
            base_emi_low(i,:)=sum(tmp_base_emi_low,1); %not border cell
        end
        
        %high emission treatment
        tmp_base_emi_high=tmp_emi_high(sort(bhigh),:);
        %order of emissions NOX, COV, NH3, PM10, PM25, SO2
        if isempty(emi{i,1})==0
            base_emi_high(i,:)=sum(tmp_base_emi_high,1)+emi{i,1}(2,:); %sum emi cell part outside optimization domain;
        else
            base_emi_high(i,:)=sum(tmp_base_emi_high,1); %sum emi cell part outside optimization domain;
        end
        
        %UNIBS(ET)20131001 - actlev considers only PAD cells
        if (flag_optim_dom(i)==1 | flag_optim_dom(i)==2)
            %activity level NOC
            actlev_noc{i}=emi_noc(:,7);
            
            %activity level
            actlev{i}=emi{i,2}(:,7);
            
            %activity level of all sector-activity-techs
            actlev_final{i}=[actlev{i}; actlev_noc{i}];
            %optimization domain wide sum of activity levels
            actlev_final_sum=actlev_final_sum+actlev_final{i};
        end
        
    end
end



