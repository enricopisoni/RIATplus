%function [aqi_per_cell]= neuralNet_do_aqi_per_cell( emissioni, aggregationInfo, commonDataInfo, periodIndex, aqiIndex)
function [aqi_per_cell]= neuralNet_do_aqi_per_cell( emissioni, NN, aggregationInfo, commonDataInfo, periodIndex, aqiIndex)

        input_rete2=emissioni';
        
        %in case it is necessary to process quadrant emissions (if too close
        %to domain boundary, it is necessary to increment emissions with
        %the assumptions that part of the quadrant in which emissions are
        %not available, still contain same emission average)
        dirs=commonDataInfo.dirs;
        domainInfo=commonDataInfo.domainInfo;

        if strcmp(dirs.pathAR,'-1')==0
            %load path for Area_ratio variable...file contain ratios  with
            %the order as NORTH SOUTH NWEST EAST
            
            %remove cells not in flag_optim_dom
            
            %repmat to get same "emissioni'"
            %dimensions
            AR=load(dirs.pathAR);                             %load not in optim
            AR=AR.Ratio;                                 %rename variable
            AR(domainInfo.flag_optim_dom==0,:)=[];                  %remove cells outside flag_optim_dom
            
            %ENR20130313
            %manage case with 2 input for ozone, or 6 input for ozone and
            %PM ANNs
            %             if NN.net.inputs{1}.size==48
            %                 ARreordrepm=repmat(AR,1,12);            %repmat
            %             elseif NN.net.inputs{1}.size==16
            %                 ARreordrepm=repmat(AR,1,4);            %repmat
            %             end
            
            if (size(NN.net.IW{1,1},2)==48)% ANN/linear, 6 input
                ARreordrepm=repmat(AR,1,12);            %repmat
                
            elseif (size(NN.net.IW{1,1},2)==16)%ANNlinear, 2 input
                ARreordrepm=repmat(AR,1,4);            %repmat
                
            end
            
            %             if size(NN.net,1)==48 % linear, 6 input
            %                 ARreordrepm=repmat(AR,1,12);            %repmat
            %
            %             elseif size(NN.net,1)==16 %linear, 2 input
            %                 ARreordrepm=repmat(AR,1,4);            %repmat
            %
            %             elseif size(NN.net,1)==1 %ANNs
            % %                 if NN.net.inputs{1}.size==48 %ANNs, 6 input
            %                     ARreordrepm=repmat(AR,1,12);            %repmat
            %
            %                 elseif NN.net.inputs{1}.size==16 %ANNs, 2 input
            %                     ARreordrepm=repmat(AR,1,4);            %repmat
            %                 end
            %             end
            
            
            %             ARreord=[AR(:,2) AR(:,3) AR(:,1) AR(:,4)];   %reorder
            %             ARreordrepm=repmat(AR,1,12);            %repmat
            input_rete2=emissioni'./ARreordrepm';        %rewrite input_rete2 dividing emissions by area_ratio
        end
        %change name to a better one!!! (too similar to caller...)
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint8_new','input_rete2');
        aqi_per_cell=interface_get_aqipercell(input_rete2, NN, aggregationInfo.mathIntermediateData, commonDataInfo, periodIndex, aqiIndex );
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\varPoint9_new','aqi_per_cell');       
        %define if lin o net

end
