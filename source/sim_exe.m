% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output_net_norm]=sim_exe(NN,input_rete_norm2)  
%20141121-ET new function to substitute matlab function sim
sz=size(input_rete_norm2,2);
%neural network simulation
        s_a=NN.net.IW{1,1}*input_rete_norm2+repmat(NN.net.b{1},1,sz);
        s_a1=eval(strcat(NN.net.layers{1}.transferFcn,'(s_a)'));
        s_b=NN.net.LW{2}*s_a1+repmat(NN.net.b{2},1,sz);
        output_net_norm=eval(strcat(NN.net.layers{2}.transferFcn,'(s_b)'));     
        
end
