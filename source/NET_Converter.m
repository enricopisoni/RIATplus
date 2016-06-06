function NET_Converter(f1)
%NET CONVERTER
%v 2.0
%#function network
file_name=f1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long
%load net
load(strcat(f1,'.mat'))
fileID=fopen(strcat(file_name,'.txt'),'wt');
fprintf('\n Lettura rete eseguita');
%old networks
if exist('Model','var')==0
    %add new variables
    Model='net';
    ArPt='Separated';
    Class='EmiConc';
    if cell2mat(net.inputs.size)==48
    PRECs={'NH3','NOX','PM10','PM25','SO2','VOC'};
    elseif cell2mat(net.inputs.size)==16
    PRECs={'NOX','VOC'};   
    elseif cell2mat(net.inputs.size)==40
    PRECs={'NH3','NOX','PM10','SO2','VOC'};
    elseif cell2mat(net.inputs.size)==8
    PRECs={'NOX','VOC'}; 
    ArPt='Aggregated';
    else
        error
    end
    ps_pca=[];
end

if strcmp(Model,'net')==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%net

%write input number
ninput=net.inputs{1}.size;
fprintf(fileID,'%u\n',ninput);
%write hidden neurons
nhidden=net.layers{1}.size;
fprintf(fileID,'%u\n',nhidden);
%write transfer function 1
tf1=net.layers{1}.transferFcn;
fprintf(fileID,'%s\n',tf1);
%write IW
IW_write=net.IW{1};
for i=1:nhidden
    fprintf(fileID,strcat(repmat('%.15f ',1,ninput),'\n'),IW_write(i,:));
end
%write b1
b1_write=net.b{1};
fprintf(fileID,strcat(repmat('%.15f ',1,nhidden),'\n'),b1_write');
%write transfer function 2
tf2=net.layers{2}.transferFcn;
fprintf(fileID,'%s\n',tf2);
%write LW
LW_write=net.LW{2};
fprintf(fileID,strcat(repmat('%.15f ',1,nhidden),'\n'),LW_write);
%write b2
b2_write=net.b{2};
fprintf(fileID,'%.15f\n',b2_write);
fprintf('\n Scrittura Struttura e Pesi eseguita');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ps variables

% write ps_input
%fprintf(fileID,'%s\n',ps_input.name);
%fprintf(fileID,'%u\n',ps_input.xrows);
fprintf(fileID,strcat(repmat('%.15f ',1,ps_input.xrows),'\n'),ps_input.xmax);
fprintf(fileID,strcat(repmat('%.15f ',1,ps_input.xrows),'\n'),ps_input.xmin);
%fprintf(fileID,strcat(repmat('%.15f ',1,ps_input.xrows),'\n'),ps_input.xrange);
%fprintf(fileID,'%u\n',ps_input.yrows);
fprintf(fileID,'%.15f\n',ps_input.ymax);
fprintf(fileID,'%.15f\n',ps_input.ymin);
%fprintf(fileID,'%.15f\n',ps_input.yrange);
%fprintf(fileID,'%u\n',ps_input.no_change);
% write ps_target
%fprintf(fileID,'%s\n',ps_target.name);
%fprintf(fileID,'%u\n',ps_target.xrows);
fprintf(fileID,strcat(repmat('%.15f ',1,ps_target.xrows),'\n'),ps_target.xmax);
fprintf(fileID,strcat(repmat('%.15f ',1,ps_target.xrows),'\n'),ps_target.xmin);
%fprintf(fileID,strcat(repmat('%.15f ',1,ps_target.xrows),'\n'),ps_target.xrange);
%fprintf(fileID,'%u\n',ps_target.yrows);
fprintf(fileID,'%.15f\n',ps_target.ymax);
fprintf(fileID,'%.15f\n',ps_target.ymin);
%fprintf(fileID,'%.15f\n',ps_target.yrange);
%fprintf(fileID,'%u\n',ps_target.no_change);
zmin=net.outputs{1,2}.processSettings{1,end}.xmin;
zmax=net.outputs{1,2}.processSettings{1,end}.xmax;
GAIN=(zmax-zmin)/2;
OFFSET=zmin;
fprintf(fileID,'%.15f\n',GAIN);
fprintf(fileID,'%.15f\n',OFFSET);

% write ps_pca
if isempty(ps_pca)==0
fprintf(fileID,'%u\n',1);    
fprintf(fileID,'%s\n',ps_pca.name);
fprintf(fileID,'%u\n',ps_pca.xrows);
fprintf(fileID,'%.15f\n',ps_pca.maxfrac);
fprintf(fileID,'%u\n',ps_pca.yrows);
for i_pca=1:ps_pca.xrows
fprintf(fileID,strcat(repmat('%.15f ',1,ps_pca.yrows),'\n'),ps_pca.transform(i_pca,:));
end
fprintf(fileID,'%u\n',ps_pca.no_change);
else
fprintf(fileID,'%u\n',0); %if 0 use isempty [] in read
end



elseif strcmp(Model,'lin')==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%linear model
 fprintf(fileID,strcat(repmat('%.15f ',1,size(net,1)),'\n'),net);     
% ps_input  
fprintf(fileID,'%d\n',-999);  
%ps_target
fprintf(fileID,'%d\n',-999);
%ps_pca
fprintf(fileID,'%d\n',-999);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%flags

% write ArPt
fprintf(fileID,'%s\n',ArPt);
% write Class
fprintf(fileID,'%s\n',Class);
% write PRECs
fprintf(fileID,'%s\n',num2str(size(PRECs,2)));
for i_prec=1:size(PRECs,2)
fprintf(fileID,'%s\n',PRECs{1,i_prec});
end
% write icells
fprintf(fileID,'%u\n',icells);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%print fileEND
fprintf(fileID,'%d\n',-9999);
fclose(fileID);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if delta write basecase files
if strcmp(Class,'Delta')==1
    
    %load concBC
    file_concBC=strcat(PATHfolder,'dataset/',TARGET,'_concBC_',period);
    load(strcat(file_concBC,'.mat'));
    fileID=fopen(strcat(file_concBC,'.txt'),'wt');
    for i=1:size(BC_conc,1)
    fprintf(fileID,strcat(repmat('%.15f ',1,size(BC_conc,2)),'\n'),BC_conc(i,:));
    end
    fclose(fileID);
    
    %load emiBC
    file_emiBC=strcat(PATHfolder,'dataset/',TARGET,'_emiBC_',period);
    load(strcat(file_emiBC,'.mat'));
    fileID=fopen(strcat(file_emiBC,'.txt'),'wt');
    for j=1:size(BC_emi,1)
    fprintf(fileID,strcat(repmat('%.15f ',1,size(BC_emi,2)),'\n'),BC_emi(j,:));
    end
    
    fclose(fileID);
    
    
end
fprintf('\n Fine Scrittura \n');
end