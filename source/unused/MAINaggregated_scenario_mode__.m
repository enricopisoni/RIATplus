    function [aqi_per_cell]=MAINaggregated_scenario_mode__(emiTMP,NN,ii,jj, commonDataInfo, aggregationInfo)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %ANNS CONSIDERED
        %NN=nnSuperSet(ii).nnSet(jj);
        %save('C:\data\work\projects\riat\RiatPlus-v6.1b\datasave\sa_start_fix', 'emiTMP','NN');
        icells=NN.icells;
        coordinate=commonDataInfo.domainData.data(:,1:2);
        stepsize=round(max(abs(coordinate(2,1)-coordinate(1,1)),abs(coordinate(2,2)-coordinate(1,2))));
        nxny=round((max(coordinate)-min(coordinate))/stepsize+1);
        areal_point=commonDataInfo.optim_flags.areal_point;
        flag_optim_dom=commonDataInfo.domainInfo.flag_optim_dom;
        pathAR=commonDataInfo.dirs.pathAR;
        nx=nxny(1,1);
        ny=nxny(1,2);
        ncel=nx*ny;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %AREAL CASE
        [s1_NOX,s2_NOX,s3_NOX,s4_NOX]=MAINquadrant_emi__(emiTMP(:,1),icells,nx,ny);
        [s1_VOC,s2_VOC,s3_VOC,s4_VOC]=MAINquadrant_emi__(emiTMP(:,2),icells,nx,ny);
        [s1_NH3,s2_NH3,s3_NH3,s4_NH3]=MAINquadrant_emi__(emiTMP(:,3),icells,nx,ny);
        [s1_PM10,s2_PM10,s3_PM10,s4_PM10]=MAINquadrant_emi__(emiTMP(:,4),icells,nx,ny);
        [s1_PM25,s2_PM25,s3_PM25,s4_PM25]=MAINquadrant_emi__(emiTMP(:,5),icells,nx,ny);
        [s1_SO2,s2_SO2,s3_SO2,s4_SO2]=MAINquadrant_emi__(emiTMP(:,6),icells,nx,ny);
        NH3_all=[s1_NH3,s2_NH3,s3_NH3,s4_NH3];
        %save('C:\data\work\projects\riat\RiatPlus-v6.1b\datasave\sa_NH3_all_fix', 'NH3_all');
        NOX_all=[s1_NOX,s2_NOX,s3_NOX,s4_NOX];
        PM10_all=[s1_PM10,s2_PM10,s3_PM10,s4_PM10];
        PM25_all=[s1_PM25,s2_PM25,s3_PM25,s4_PM25];
        SO2_all=[s1_SO2,s2_SO2,s3_SO2,s4_SO2];
        VOC_all=[s1_VOC,s2_VOC,s3_VOC,s4_VOC];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CONSIDER CASE WITH 6 OR 2 PERCURSOR EMISSIONS AS INPUT,
        %AREAL CASE. ALWAYS THERE ARE 4 QUADRANTS TO CONSIDER WIND
        %DIRETCIONS
        if areal_point==0
            %read net and select precursors
            emissioni=[];
            for i_prec=1:size(NN.PRECs,2) %areal
                if strcmp(NN.PRECs(1,i_prec),'NH3')==1
                    emissioni=[emissioni, NH3_all];
                elseif strcmp(NN.PRECs(1,i_prec),'NOX')==1
                    emissioni=[emissioni, NOX_all];
                elseif strcmp(NN.PRECs(1,i_prec),'PM10')==1
                    emissioni=[emissioni, PM10_all];
                elseif strcmp(NN.PRECs(1,i_prec),'PM25')==1
                    emissioni=[emissioni, PM25_all];
                elseif strcmp(NN.PRECs(1,i_prec),'SO2')==1
                    emissioni=[emissioni, SO2_all];
                elseif strcmp(NN.PRECs(1,i_prec),'VOC')==1
                    emissioni=[emissioni, VOC_all];
                end
            end
            
            
            %             if (size(NN.net,1)==24 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==24))% ANNlinear, 6 input
            %                 emissioni=[NH3_all,NOX_all,PM10_all,PM25_all,SO2_all,VOC_all];
            %             elseif (size(NN.net,1)==8 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==8))%ANNlinear, 2 input
            %                 emissioni=[NOX_all,VOC_all];
            %             end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %IF ALSO POINT SOURCES
        elseif areal_point==1
            %point
            [s1_NOXp,s2_NOXp,s3_NOXp,s4_NOXp]=MAINquadrant_emi__(emiTMP(:,7),icells,nx,ny);
            [s1_VOCp,s2_VOCp,s3_VOCp,s4_VOCp]=MAINquadrant_emi__(emiTMP(:,8),icells,nx,ny);
            [s1_NH3p,s2_NH3p,s3_NH3p,s4_NH3p]=MAINquadrant_emi__(emiTMP(:,9),icells,nx,ny);
            [s1_PM10p,s2_PM10p,s3_PM10p,s4_PM10p]=MAINquadrant_emi__(emiTMP(:,10),icells,nx,ny);
            [s1_PM25p,s2_PM25p,s3_PM25p,s4_PM25p]=MAINquadrant_emi__(emiTMP(:,11),icells,nx,ny);
            [s1_SO2p,s2_SO2p,s3_SO2p,s4_SO2p]=MAINquadrant_emi__(emiTMP(:,12),icells,nx,ny);
            NH3_allp=[s1_NH3p,s2_NH3p,s3_NH3p,s4_NH3p];
            %save('C:\data\work\projects\riat\RiatPlus-v6.1b\datasave\sa_NH3_allp_fix', 'NH3_allp');
            NOX_allp=[s1_NOXp,s2_NOXp,s3_NOXp,s4_NOXp];
            PM10_allp=[s1_PM10p,s2_PM10p,s3_PM10p,s4_PM10p];
            PM25_allp=[s1_PM25p,s2_PM25p,s3_PM25p,s4_PM25p];
            SO2_allp=[s1_SO2p,s2_SO2p,s3_SO2p,s4_SO2p];
            VOC_allp=[s1_VOCp,s2_VOCp,s3_VOCp,s4_VOCp];
            
            
            
            %create input structure, to be used in the ANNs
            %read net and select precursors
            emissioni=[];
            for i_prec=1:size(NN.PRECs,2) %areal
                if strcmp(NN.PRECs(1,i_prec),'NH3')==1
                    emissioni=[emissioni, NH3_all];
                elseif strcmp(NN.PRECs(1,i_prec),'NOX')==1
                    emissioni=[emissioni, NOX_all];
                elseif strcmp(NN.PRECs(1,i_prec),'PM10')==1
                    emissioni=[emissioni, PM10_all];
                elseif strcmp(NN.PRECs(1,i_prec),'PM25')==1
                    emissioni=[emissioni, PM25_all];
                elseif strcmp(NN.PRECs(1,i_prec),'SO2')==1
                    emissioni=[emissioni, SO2_all];
                elseif strcmp(NN.PRECs(1,i_prec),'VOC')==1
                    emissioni=[emissioni, VOC_all];
                end
            end
            for i_prec=1:size(NN.PRECs,2) %point
                
                if strcmp(NN.PRECs(1,i_prec),'NH3')==1
                    emissioni=[emissioni, NH3_allp];
                elseif strcmp(NN.PRECs(1,i_prec),'NOX')==1
                    emissioni=[emissioni, NOX_allp];
                elseif strcmp(NN.PRECs(1,i_prec),'PM10')==1
                    emissioni=[emissioni, PM10_allp];
                elseif strcmp(NN.PRECs(1,i_prec),'PM25')==1
                    emissioni=[emissioni, PM25_allp];
                elseif strcmp(NN.PRECs(1,i_prec),'SO2')==1
                    emissioni=[emissioni, SO2_allp];
                elseif strcmp(NN.PRECs(1,i_prec),'VOC')==1
                    emissioni=[emissioni, VOC_allp];
                end
            end
            
            %             if (size(NN.net,1)==48 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==48))% ANN/linear, 6 input
            %                 emissioni=[NH3_all,NOX_all,PM10_all,PM25_all,SO2_all,VOC_all,...
            %                     NH3_allp,NOX_allp,PM10_allp,PM25_allp,SO2_allp,VOC_allp];
            %             elseif (size(NN.net,1)==16 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==16))%ANNlinear, 2 input
            %                 emissioni=[NOX_all,VOC_all,NOX_allp,VOC_allp];
            %             end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %ANN INPUT DATA
        %keep only optim domain
        emissioni(find(flag_optim_dom==0),:)=[];
        %20130820 - consider only if cell completerly in PAD
        %         emissioni(find(flag_optim_dom==0 | flag_optim_dom==2),:)=[];
        input_rete2=emissioni';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %IN CASE QUADRANT EMISSIONS TOO CLOSE TO CTM BOUNDARY
        if strcmp(pathAR,'-1')==0
            AR=load(pathAR);                             %load not in optim
            AR=AR.Ratio;                                 %rename variable
            AR(flag_optim_dom==0,:)=[];                  %remove cells outside flag_optim_dom
            
            if (size(NN.net.IW{1,1},2)==48)% ANN/linear, 6 input
                ARreordrepm=repmat(AR,1,12);            %repmat
                
            elseif (size(NN.net.IW{1,1},2)==16)%ANNlinear, 2 input
                ARreordrepm=repmat(AR,1,4);            %repmat
            end
            
            %             if (size(NN.net,1)==48 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==48))% ANN/linear, 6 input
            %                 ARreordrepm=repmat(AR,1,12);            %repmat
            %             elseif (size(NN.net,1)==16 || (size(NN.net,1)==1 && NN.net.inputs{1}.size==16))%ANNlinear, 2 input
            %                 ARreordrepm=repmat(AR,1,4);            %repmat
            %             end
            
            input_rete2=emissioni'./ARreordrepm';        %rewrite input_rete2 dividing emissions by area_ratio
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %CHECK IF EMISSIONS ARE INSIDE ANNS IDENTIFICATION
        %BOUNDS...IF NOT, STOP RIAT+
        %compare ANNs input and ps_input, to see if we are in the
        %ANNs bounds
        %check in the PAD
        %         for indcel=1:ncel
        %                 if sum(input_rete2(:,indcel)<NN.ps_input.xmin)>0 | sum(input_rete2(:,indcel)>NN.ps_input.xmax)>0
        %                     find(input_rete2(:,indcel)<NN.ps_input.xmin)
        %                     find(input_rete2(:,indcel)>NN.ps_input.xmax)
        %                     figure;plot([input_rete2(:,indcel)-NN.ps_input.xmin  NN.ps_input.xmax-input_rete2(:,indcel)])
        %                     error('SR model identification bounds not respected...RIAT+ is terminated')
        %                 end
        %         end
        %matrix of min and max ANNS values
        minVal=repmat(NN.ps_input.xmin,1,size(emissioni,1));
        maxVal=repmat(NN.ps_input.xmax,1,size(emissioni,1));
        %save('C:\data\work\projects\riat\RiatPlus-v6.1b\datasave\sa_minmax_fix', 'minVal', 'maxVal');
        %check if my scenario is out of bounds
        %         indmin=find(input_rete2<minVal);  OK
        %         indmax=find(input_rete2>maxVal);  OK
        %20130828 - input_rete2>0 added to manage issue of the Alsace cells
        %adjacent to the CTM grid
        indmin=find(input_rete2<minVal & input_rete2>0);
        indmax=find(input_rete2>maxVal & input_rete2>0);
        %check bounds problem.....if limited problems, simulate the
        %scenario, otherwise exit
        %         ii
        %         jj
        %         max((input_rete2(indmax)-maxVal(indmax))./(maxVal(indmax))*100)
        %         min((input_rete2(indmin)-minVal(indmin))./minVal(indmin)*100)
        
        if length(indmax)>0
            %             if max(input_rete2(indmax)-maxVal(indmax))>1
            if max((input_rete2(indmax)-maxVal(indmax))./(maxVal(indmax))*100)>30
                
                Error=3;
                
                fprintf(commonDataInfo.fidExit, int2str(Error));
                strStatus='PROGRESSION: Aggregated scenario analysis finished.';
                disp(strStatus);
                fprintf(commonDataInfo.fidStatus, '%s\n',strStatus);
                error('SR model identification bounds not respected...RIAT+ is terminated')
            end
        end
        
        if length(indmin)>0
            %             if max(input_rete2(indmin)-minVal(indmin))<-1
            if min((input_rete2(indmin)-minVal(indmin))./minVal(indmin)*100)<-30
                
                Error=3;
                fprintf(commonDataInfo.fidExit, int2str(Error));
                strStatus='PROGRESSION: Aggregated scenario analysis finished.';
                disp(strStatus);
                fprintf(commonDataInfo.fidStatus, '%s\n',strStatus);
                error('SR model identification bounds not respected...RIAT+ is terminated')
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %IF BOUNDS RESPECTED, RUN SR MODEL
        %if constraints respected for this AQI, create results
        %define if lin o net
        if size(NN.net,1)==1
            flag_reg_net=1;
        elseif size(NN.net,1)>1
            flag_reg_net=0;
        end
        
        switch flag_reg_net
            case 0 %linear case
                aqi_per_cell=(input_rete2'*NN.net)';
            case 1 %nonlinear case
                %20111215 - run the network only on the optimization cells
                
                %20141205-ET - mapminmax replaced with a handwritten normalization
                %to deal with the problems of standalone windows executables
                ymin=NN.ps_input.ymin;
                ymax=NN.ps_input.ymax;
                xmin=NN.ps_input.xmin;
                xmax=NN.ps_input.xmax;
                G1=(ymax-ymin)./(xmax-xmin);
                
                input_rete_norm2=((input_rete2-repmat(xmin,1,size(input_rete2,2))).*repmat(G1,1,size(input_rete2,2))+repmat(ymin,size(input_rete2,1),size(input_rete2,2)));
                %                 NN.ps_input.no_change=0; %to be compatible among different matlab versions
                %                 [input_rete_norm2]=mapminmax('apply',input_rete2,NN.ps_input);
                
                %check for results outside bounds - force to be inside bounds
                input_rete_norm2(find(input_rete_norm2<-1))=-1;
                input_rete_norm2(find(input_rete_norm2>1))=1;
                
                %20141121-ET - Use new sim_exe function instead of matlab
                %sim function
                output_rete_norm2=sim_exe(NN,input_rete_norm2);
                
                %20141205-ET - mapminmax replaced with a handwritten normalization
                %to deal with the problems of standalone windows
                %executables, matlab transformation included
                G=NN.ps_target.gain;
                offst=NN.ps_target.offset;
                ymin=NN.ps_target.ymin;
                ymax=NN.ps_target.ymax;
                xmin=NN.ps_target.xmin;
                xmax=NN.ps_target.xmax;
                
                aqi_per_cell=((output_rete_norm2-ymin)*G*(xmax-xmin))/(ymax-ymin)+xmin+((xmax-xmin)*(-ymin+offst)/(ymax-ymin));
                %                 NN.ps_target.no_change=0; %to be compatible among different matlab versions
                %                 aqi_per_cell=mapminmax('reverse',output_rete_norm2,NN.ps_target);
                
                %%%
                %check for results outside bounds - force to be inside bounds
                aqi_per_cell(find(aqi_per_cell<xmin))=xmin;
                aqi_per_cell(find(aqi_per_cell>xmax))=xmax;
                %%%
                
                
        end
    end