%20140617 ET - read model from txt file
function [nn]=net_read(file_in)
%file_in
fid=fopen(strcat(file_in,'.txt'));

%linear or non linear model
%read input number
riga=fgetl(fid);
ninput=sscanf(riga,'%d');

if size(ninput)==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%net

%read hidden neurons
riga=fgetl(fid);
nhidden=sscanf(riga,'%d');
riga=fgetl(fid);
%read transfer function 1
tf1=strtrim(riga);
for i=1:nhidden
    riga=fgetl(fid);
    appo=sscanf(riga,'%f');
    IW_read(i,:)=appo;
end
%read b1
b1_read=sscanf(fgetl(fid),'%f');
%read transfer function 2
riga=fgetl(fid);
tf2=strtrim(riga);
%read LW
riga=fgetl(fid);
LW_read=sscanf(riga,'%f');
%read b2
b2_read=sscanf(fgetl(fid),'%f');

%20141121-ET - No more networks, only structures from now on!
%create net
net=struct;
%set net parameters
net.layers{1}.transferFcn=tf1;
net.IW{1}=IW_read;
net.layers{2}.transferFcn=tf2;
net.LW{2}=LW_read';
net.b{1}=b1_read;
net.b{2}=b2_read;
net.inputs{1}.range=repmat([-1 1],ninput,1); % che impatto ha????
net.outputs{2}.range=[-1 1]; % che impatto ha????

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ps variables
%20141205-ET - removed mapminmax lines,not used with new normalization
%and added new lines to deal with matlab transformation inside sim function
% read ps_input
%nn.ps_input.name=strtrim(fgetl(fid));
%nn.ps_input.xrows=sscanf(fgetl(fid),'%d');
nn.ps_input.xmax=sscanf(fgetl(fid),'%f');
nn.ps_input.xmin=sscanf(fgetl(fid),'%f');
%nn.ps_input.xrange=sscanf(fgetl(fid),'%f');
%nn.ps_input.yrows=sscanf(fgetl(fid),'%d');
nn.ps_input.ymax=sscanf(fgetl(fid),'%f');
nn.ps_input.ymin=sscanf(fgetl(fid),'%f');
%nn.ps_input.yrange=sscanf(fgetl(fid),'%f');
%nn.ps_input.no_change=sscanf(fgetl(fid),'%d');
% read ps_target
%nn.ps_target.name=strtrim(fgetl(fid));
%nn.ps_target.xrows=sscanf(fgetl(fid),'%d');
nn.ps_target.xmax=sscanf(fgetl(fid),'%f');
nn.ps_target.xmin=sscanf(fgetl(fid),'%f');
%nn.ps_target.xrange=sscanf(fgetl(fid),'%f');
%nn.ps_target.yrows=sscanf(fgetl(fid),'%d');
nn.ps_target.ymax=sscanf(fgetl(fid),'%f');
nn.ps_target.ymin=sscanf(fgetl(fid),'%f');
%nn.ps_target.yrange=sscanf(fgetl(fid),'%f');
%nn.ps_target.no_change=sscanf(fgetl(fid),'%d');
nn.ps_target.gain=sscanf(fgetl(fid),'%f');
nn.ps_target.offset=sscanf(fgetl(fid),'%f');
% read ps_pca
flag_pca=sscanf(fgetl(fid),'%d');
if flag_pca==1
nn.ps_pca.name=strtrim(fgetl(fid));
nn.ps_pca.xrows=sscanf(fgetl(fid),'%d');
nn.ps_pca.maxfrac=sscanf(fgetl(fid),'%f');
nn.ps_pca.yrows=sscanf(fgetl(fid),'%d');
for i_pca=1:ps_pca.yrows
nn.ps_pca.transform(i_pca,:)=sscanf(fgetl(fid),'%f');
end
nn.ps_pca.no_change=sscanf(fgetl(fid),'%d');
elseif flag_pca==0
nn.ps_pca=[];
end



elseif size(ninput)~=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%linear model
net=ninput;
nn.ps_input=sscanf(fgetl(fid),'%d');
nn.ps_target=sscanf(fgetl(fid),'%d');
nn.ps_pca=sscanf(fgetl(fid),'%d');

end
nn.net=net;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%flags

% read ArPt
riga=fgetl(fid);
nn.ArPt=strtrim(riga);
% read Class
riga=fgetl(fid);
nn.Class=strtrim(riga);
% read PRECs
riga=fgetl(fid);
n_precs=sscanf(riga,'%f');
% read PRECs
for i_prec=1:n_precs
riga=fgetl(fid);
nn.PRECs{1,i_prec}=strtrim(riga);
end
% read icells
riga=fgetl(fid);
nn.icells=sscanf(riga,'%f');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check if readed correctly
fileEND=sscanf(fgetl(fid),'%d');
if fileEND~=-9999
    error 'Problems reading the net file'
else
    disp('Net loaded correctly')
end

fclose(fid);



end
