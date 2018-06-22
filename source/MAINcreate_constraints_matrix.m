function [A,B] = ...
    MAINcreate_constraints_matrix(global_data,CLE,flag_constraints,pathLST)

[technologiesNum, foo] = size(global_data);

% Select sec-act-Low|High
sec_act = unique(global_data(:,[2 3 5]),'rows');
[sec_act_num,foo] = size(sec_act);

% Removal efficiencies with respect to the six pollutants
eff = global_data(:, 6:11);

% distinguish between technical and non-technical measures:
% technical measures: flag_tech_nontech = 1
% non-technical measures: flag_tech_nontech = 0
flag_tech_nontech=global_data(:,23);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         TECHNICAL MEASURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           NH3 constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The NH3 constraints are created on-the-fly in the code, depending on
% the technologies codes defined in the file "./input/techsList.csv"

% List of NH3 technologies to be involved in NH3 constraints creation
NH3technologies = ...
    {'SA','LNF','LNA','CS','BF'};
[foo, NH3technologiesNum] = size(NH3technologies);

% Mapping among technologies' ids and technologies' names
[technologiesNameCodeMapping]=importdata(pathLST);
[technologiesNameCodeMappingNum,foo] = ...
    size(technologiesNameCodeMapping.textdata);

[Code I] = sort(technologiesNameCodeMapping.data);
technologiesNameCodeMappingFull(Code) = ...
    technologiesNameCodeMapping.textdata(I);

technologiesNameCodeMappingFull = technologiesNameCodeMappingFull';

%fill missing data with NAN
%for k=1:size(technologiesNameCodeMappingFull,1)
%    if isempty(technologiesNameCodeMappingFull{k})==1
%        technologiesNameCodeMappingFull{k}='NAN';
%    end
%end

A_NH3 = [];
B_NH3 = [];

for sa=1:sec_act_num,
    
    for NH3tt=1:NH3technologiesNum,
        
        NH3ttName = NH3technologies(NH3tt);
        
        constraintLHS = [];
        constraintRHS = 1;
        
        for tt=1:technologiesNum
            
            technologyId = global_data(tt,4);
            
            
            if((sec_act(sa,1) == global_data(tt,2)) && ...
                    (sec_act(sa,2) == global_data(tt,3)) && ...
                    (sec_act(sa,3) == global_data(tt,5)) && ...
                    (flag_tech_nontech(tt) == 1) )

                ttName = technologiesNameCodeMappingFull{technologyId};
                
                pat = sprintf('(^%s_\\w*|\\w*_%s$|\\w*_%s_\\w*|^%s$)',...
                    NH3ttName{1},NH3ttName{1},NH3ttName{1},NH3ttName{1});
                m = regexp(ttName, pat, 'match');
                
                
                %                 if((NH3tt==3) || (NH3tt==6))
                %                     internalNH3ttName = NH3technologies(NH3tt+1);
                %                     pat = sprintf('\\w*%s\\w*',internalNH3ttName{1});
                %                     m1 = regexp(ttName, pat, 'match');
                %
                %                     internalNH3ttName = NH3technologies(NH3tt+2);
                %                     pat = sprintf('\\w*%s\\w*',internalNH3ttName{1});
                %                     m2 = regexp(ttName, pat, 'match');
                %
                %                     if((isempty(m1)) && (isempty(m2)))
                %                         if(isempty(m))
                %                             constraintLHS = [constraintLHS, 0];
                %                         else
                %                             constraintLHS = [constraintLHS, 1];
                %                             if(isequal(m(1), NH3ttName))
                %                                 constraintRHS = global_data(tt,20)/100;
                %                             end
                %                         end
                %                     else
                %                         constraintLHS = [constraintLHS, 0];
                %                     end;
                %                 else
                if(isempty(m))
                    constraintLHS = [constraintLHS, 0];
                else
                    constraintLHS = [constraintLHS, 1];
                    if(isequal(m(1), NH3ttName))
                        constraintRHS = global_data(tt,20)/100;
                    end
                end
                %                end
            else
                constraintLHS = [constraintLHS, 0];
            end
        end
        
        if(sum(constraintLHS) > 1)
            A_NH3 = [A_NH3 ; constraintLHS];
            B_NH3 = [B_NH3 ; constraintRHS];
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTRAINTS CONCERNING THE MUTUAL EXCLUSION OF TECHNOLOGIES APPLICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sum AR < 1

A_T_MUTEX = [];
B_T_MUTEX = [];

for sa=1:sec_act_num,
    for p=1:6,
        
        constraintLHS = [];
        constraintRHS = 1;
        
        for tt=1:technologiesNum
            
            if((sec_act(sa,1) == global_data(tt,2)) && ...
                    (sec_act(sa,2) == global_data(tt,3)) && ...
                    (sec_act(sa,3) == global_data(tt,5)) && ...
                    ((eff(tt,p) > 0) || (eff(tt,p) < 0)) && ...
                    (flag_tech_nontech(tt) == 1) )
                
                constraintLHS = [constraintLHS, 1];
            else
                constraintLHS = [constraintLHS, 0];
            end
            
        end
        
        if(sum(constraintLHS) > 0)
            A_T_MUTEX = [A_T_MUTEX ; constraintLHS];
            B_T_MUTEX = [B_T_MUTEX ; constraintRHS];
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        CONTRAINTS TO IMPOSE IN CASE OF REPLACEABLE TECHNOLOGIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sum eff*AR > sum eff*ARcle; i.e. -sum eff*AR < -sum eff*ARcle
A_EMISSIONS_REDUCTION = [];
B_EMISSIONS_REDUCTION = [];

% sum AR > sum ARcle; i.e. -sum AR < -sum ARcle
A_EMISSIONS_CONTROLLED = [];
B_EMISSIONS_CONTROLLED = [];


if flag_constraints==1
    
    % -sum eff*AR < -sum eff*ARcle
%    for sa=1:sec_act_num,
%        for p=1:6,
            
%            constraintLHS = [];
%            constraintRHS = 0;
            
%            insert = 0;
            
%            for tt=1:technologiesNum
                
%                if((sec_act(sa,1) == global_data(tt,2)) && ...
%                        (sec_act(sa,2) == global_data(tt,3)) && ...
%                        (sec_act(sa,3) == global_data(tt,5)) && ...
%                       ((eff(tt,p) > 0) || (eff(tt,p) < 0)) && ...
%                       (flag_tech_nontech(tt) == 1) )
                    
                    
%                    constraintLHS = [constraintLHS, -eff(tt,p)/100];
                    
%                   constraintRHS = constraintRHS - ...
%                        CLE(tt) * eff(tt,p)/100;
                    
%                    insert = 1;
%                else
%                    constraintLHS = [constraintLHS, 0];
%                end
                
%            end
            
%            if(insert == 1)
%                A_EMISSIONS_REDUCTION = [A_EMISSIONS_REDUCTION ; constraintLHS];
%                B_EMISSIONS_REDUCTION = [B_EMISSIONS_REDUCTION ; constraintRHS];
%            end
%       end
%    end
    
    
    % -sum AR < -sum ARcle
    for sa=1:sec_act_num,
        for p=1:6,
            
            constraintLHS = [];
            constraintRHS = 0;
            
            for tt=1:technologiesNum
                
                if((sec_act(sa,1) == global_data(tt,2)) && ...
                        (sec_act(sa,2) == global_data(tt,3)) && ...
                        (sec_act(sa,3) == global_data(tt,5)) && ...
                        ((eff(tt,p) > 0) || (eff(tt,p) < 0)) && ...
                        (flag_tech_nontech(tt) == 1) )
                    
                    constraintLHS = [constraintLHS, -1];
                    
                    constraintRHS = constraintRHS - ...
                        CLE(tt);
                else
                    constraintLHS = [constraintLHS, 0];
                end
                
            end
            
            if(sum(constraintLHS) < 0)
                A_EMISSIONS_CONTROLLED = [A_EMISSIONS_CONTROLLED ; constraintLHS];
                B_EMISSIONS_CONTROLLED = [B_EMISSIONS_CONTROLLED ; constraintRHS];
            end
            
        end
    end
    
end

A = [A_NH3; A_T_MUTEX; A_EMISSIONS_REDUCTION; A_EMISSIONS_CONTROLLED];
B = [B_NH3; B_T_MUTEX; B_EMISSIONS_REDUCTION; B_EMISSIONS_CONTROLLED];

% A = [A_T_MUTEX; A_EMISSIONS_REDUCTION; A_EMISSIONS_CONTROLLED];
% B = [B_T_MUTEX; B_EMISSIONS_REDUCTION; B_EMISSIONS_CONTROLLED];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         NON-TECHNICAL MEASURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONTRAINTS CONCERNING THE MUTUAL EXCLUSION OF TECHNOLOGIES APPLICATION
% REMOVED !!! ALL THE NTM CAN BE SUPERIMPOSED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sum AR < 1

A_NT_MUTEX = [];
B_NT_MUTEX = [];

% for sa=1:sec_act_num,
%     for p=1:6,
%         
%         constraintLHS = [];
%         constraintRHS = 1;
%         
%         for tt=1:technologiesNum
%             
%             if((sec_act(sa,1) == global_data(tt,2)) && ...
%                     (sec_act(sa,2) == global_data(tt,3)) && ...
%                     (sec_act(sa,3) == global_data(tt,5)) && ...
%                     ((eff(tt,p) > 0) || (eff(tt,p) < 0)) && ...
%                     (flag_tech_nontech(tt) == 0) )
%                 
%                 constraintLHS = [constraintLHS, 1];
%             else
%                 constraintLHS = [constraintLHS, 0];
%             end
%             
%         end
%         
%         if(sum(constraintLHS) > 0)
%             A_NT_MUTEX = [A_NT_MUTEX ; constraintLHS];
%             B_NT_MUTEX = [B_NT_MUTEX ; constraintRHS];
%         end
%         
%     end
% end

A = [A; A_NT_MUTEX];
B = [B; B_NT_MUTEX];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    TECHNICAL/NON-TECHNICAL MEASURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    CONSERVATION OF MASS CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


A_CM = [];
B_CM = [];

%constraint is activated only if there are some non technical measures:
if isempty(find(flag_tech_nontech==0))==0
    for sa=1:sec_act_num,
        for p=1:6,
            
            constraintLHS = [];
            constraintRHS = 1;
            
            insert = 0;
            
            for tt=1:technologiesNum
                
                if((sec_act(sa,1) == global_data(tt,2)) && ...
                        (sec_act(sa,2) == global_data(tt,3)) && ...
                        (sec_act(sa,3) == global_data(tt,5)) && ...
                        ((eff(tt,p) > 0) || (eff(tt,p) < 0)) )
                    
                    constraintLHS = [constraintLHS, eff(tt,p)/100];
                    insert = 1;
                else
                    constraintLHS = [constraintLHS, 0];
                end
                
            end
            
            if(insert == 1)
                A_CM = [A_CM ; constraintLHS];
                B_CM = [B_CM ; constraintRHS];
            end
            
        end
    end
end

A = [A; A_CM];
B = [B; B_CM];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            POST-PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove redundant constraints
[A I J]=unique(A,'rows');
B = B(I);



