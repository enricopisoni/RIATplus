function [output XMin XMax yMin yMax bint] = Regression_MinMax(PrecPatch,indicReg,pcaFlag)

X=PrecPatch;
y=indicReg;

%normalize data
[n, m]=size(X);
XMax=max(X);
XMin=min(X);
XNorm=(X - repmat(XMin,[n 1])) ./ (repmat(XMax,[n 1])-repmat(XMin,[n 1]));
XNorm(isnan(XNorm))=0; %remove nan
XNorm(isfinite(XNorm)==0)=0; %remove nan
yMax=max(y);
yMin=min(y);
yNorm=(y - repmat(yMin,[n 1])) ./ (repmat(yMax,[n 1]) - repmat(yMin,[n 1]));

% mdl=fitlm(PCAScores(:,1:ind1), yNorm,'linear','Intercept',false,'RobustOpts','on');
% output=PCALoadings(:,1:ind1)*mdl.Coefficients.Estimate;

if pcaFlag==0
    [output bint]=regress(y,X);
elseif pcaFlag==1
    %GIVE LESS IMPORTANCE TO INPUT CLOSE TO 0, IN PCR
%     [PCALoadings,PCAScores,PCAVar] = pca(XNorm,'Centered',false);
%     percent_explained = 100*PCAVar/sum(PCAVar);
%     percent_explained_cum=cumsum(percent_explained);
%     ind2=find(percent_explained_cum>95);
%     ind1=ind2(1);
%     [b bint]=regress(yNorm,PCAScores(:,1:ind1));
%     output=PCALoadings(:,1:ind1)*b;
    [output bint]=regress(yNorm,XNorm);
end






