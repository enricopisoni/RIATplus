%function [ aqi_per_cell ]= neuralNet_get_aqipercell(input_rete2, refData, commonDataInfo, periodIndex, aqiIndex)
function [ aqi_per_cell ]= neuralNet_get_aqipercell(input_rete2, NN, refData, commonDataInfo, periodIndex, aqiIndex)

%NN=refData.nnSuperSet(periodIndex).nnSet(aqiIndex);

if size(NN.net,1)==1
    flag_reg_net=1;
else
    flag_reg_net=0;
end

switch flag_reg_net
    case 0 %linear case
        icel=0;
        input_rete=input_rete2';
        for dimen1=1:size(input_rete,2)
            for dimen2=dimen1:size(input_rete,2)
                icel=icel+1;
                TMPnet(:,icel)=input_rete(:,dimen1).*input_rete(:,dimen2);
            end
        end
        aqi_per_cell=NN.net'*[input_rete TMPnet]';
        %if delta net transform net output in abs values
        if strcmp(NN.Class,'Delta')==1
            aqi_per_cell=BCconc'+(aqi_per_cell.*BCconc');
        end
        myNet=NN.net;
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\MAINcompute_aqi_per_cell_3_new','myNet','aqi_per_cell');
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
        myNetNew=NN;
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\MAINcompute_aqipercell_3_new','myNetNew', 'input_rete_norm2');
        %NN.ps_input.no_change=0; %to be compatible among different matlab versions
        %[input_rete_norm2]=mapminmax('apply',input_rete2,NN.ps_input);
        
        %%%
        %check for results outside bounds - force to be inside bounds
        input_rete_norm2(find(input_rete_norm2<-1))=-1;
        input_rete_norm2(find(input_rete_norm2>1))=1;
        %%%
        
        %20141121-ET - Use new sim_exe function instead of matlab
        %sim function
        output_rete_norm2=sim_exe(NN,input_rete_norm2);
        %%%
        
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
        
        %check for results outside bounds - force to be inside bounds
        %check for results outside bounds - force to be inside bounds
        aqi_per_cell(find(aqi_per_cell<xmin))=xmin;
        aqi_per_cell(find(aqi_per_cell>xmax))=xmax;
        %%%
        
        %30130402  ET -  if delta net transform net output in abs values
        if strcmp(NN.Class,'Delta')==1
            aqi_per_cell=BCconc'+(aqi_per_cell.*BCconc');
        end
        %save('C:\data\work\projects\riat\RiatPlus-v3beta\datasave\MAINcompute_aqi_per_cell_nl_new','aqi_per_cell','output_rete_norm2');
end

end
