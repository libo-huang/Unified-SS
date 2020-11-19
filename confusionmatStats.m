function stats = confusionmatStats(group,grouphat)
% INPUT
% group = true class labels
% grouphat = predicted class labels
%
% OUTPUT
% stats is a structure array
% stats.confusionMat
%               Predicted Classes
%                    p'    n'
%              ___|_____|_____| 
%       Actual  p |     |     |
%      Classes  n |     |     |
% TP: true positive, TN: true negative, 
% FP: false positive, FN: false negative
% 
% ref: https://en.wikipedia.org/wiki/Confusion_matrix
% 

fCMAT = 'confusionMat';
if nargin < 2
    vCMAT = group;
else
    vCMAT = confusionmat(group,grouphat);
end

numOfClasses = size(vCMAT,1);
totalSamples = sum(sum(vCMAT));
    
fTACC = 'TOTAL_ACC';
vTACC = trace(vCMAT) / totalSamples;

[TP,TN,FP,FN,vTPR,vTNR,vPPV,vNPV,vFNR,vFPR,vFDR,vFOR,vTS,vACC,vF1,vMCC,vBM,vMK]...
    = deal(zeros(numOfClasses,1));
for i_class = 1:numOfClasses
    TP(i_class) = vCMAT(i_class,i_class);
    auxCMat = vCMAT;
    auxCMat(:,i_class) = []; % remove column
    auxCMat(i_class,:) = []; % remove row
    TN(i_class) = sum(sum(auxCMat));
    FP(i_class) = sum(vCMAT(:,i_class))-TP(i_class);
    FN(i_class) = sum(vCMAT(i_class,:))-TP(i_class);
    
    vTPR(i_class) = TP(i_class) / (TP(i_class)+FN(i_class)); % True Positive Rate
    vTNR(i_class) = TN(i_class) / (TN(i_class)+FP(i_class)); % True Negative Rate
    vPPV(i_class) = TP(i_class) / (TP(i_class)+FP(i_class)); % Positive Predictive Value
    vNPV(i_class) = TN(i_class) / (TN(i_class)+FN(i_class)); % Negative Predictive Value
    vFNR(i_class) = FN(i_class) / (FN(i_class)+TP(i_class)); % False Negative Rate
    vFPR(i_class) = FP(i_class) / (FP(i_class)+TN(i_class)); % False Pasitive Rate
    vFDR(i_class) = FP(i_class) / (FP(i_class)+TP(i_class)); % False Discovery Rate
    vFOR(i_class) = FN(i_class) / (FN(i_class)+TN(i_class)); % False Omission Rate
    vTS(i_class) = TP(i_class) / (TP(i_class)+FN(i_class)+FP(i_class)); % Threat Score
    vACC(i_class) = (TP(i_class)+TN(i_class)) / (TP(i_class)+TN(i_class)+FP(i_class)+FN(i_class)); % Accuracy
    vF1(i_class) = 2*TP(i_class) / (2*TP(i_class)+FP(i_class)+FN(i_class)); % F1 score
    vMCC(i_class) = (TP(i_class)*TN(i_class)-FP(i_class)*FN(i_class)) / sqrt((TP(i_class)+FP(i_class))*(TP(i_class)+FN(i_class))*(TN(i_class)+FP(i_class))*(TN(i_class)+FN(i_class))); % Matthews Correlation Coefficient
    vBM(i_class) = vTPR(i_class) + vTNR(i_class) - 1; % Bookaker Informedness
    vMK(i_class) = vPPV(i_class) + vNPV(i_class) - 1; % Markedness
end
[fTPR,fTNR,fPPV,fNPV,fFNR,fFPR,fFDR,fFOR,fTS,fACC,fF1,fMCC,fBM,fMK] = ...
deal('TPR','TNR','PPV','NPV','FNR','FPR','FDR','FOR','TS','ACC','F1','MCC','BM','MK');
stats = struct(fTACC,vTACC,fCMAT,vCMAT,fTPR,vTPR,fTNR,vTNR,fPPV,vPPV,fNPV,vNPV,fFNR,vFNR,fFPR,vFPR,...
    fFDR,vFDR,fFOR,vFOR,fTS,vTS,fACC,vACC,fF1,vF1,fMCC,vMCC,fBM,vBM,fMK,vMK);

if any(group<=0)
    [otl_raw,otl_dtc,otl_hit] = outlinerCount(group(:),grouphat(:));
    stats.outlierRaw = otl_raw;
    stats.outlierDetect = otl_dtc;
    stats.outlierHit = otl_hit;
end
end

% outlinerCount(group(:),grouphat(:))
function [otl_raw,otl_dtc,otl_hit] = outlinerCount(group,grouphat)
g = group<=0;
gh = grouphat<=0;
otl_hit = sum(and(g,gh));
otl_raw = sum(g);
otl_dtc = sum(gh);
end