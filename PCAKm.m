function [vecG,W,objV_buff] = PCAKm(X, k, m)
% function [vecG,W,objV_buff,proOlp] = PCAKm(X, k)
%PCAKm
% X: data matrix, d_variables * n_observations
% k: Number of clusters
% 
% vecG_buff: label vector buffer, x_variables * n_observations
% W: transforms matrix, m * k
% objV_buff: optimal value buffer, x_variables * 1
warning off;
if ~exist('m','var')
    m = k-1; 
end
nSmp = size(X,2);
iterMax = 100;
% centralization ( it is a litte time-consuming, instead by X=X-repmat(mean(X,2),1,nSmp) or better by X=bsxfun(@minus,X,mean(X,2)) )
X=bsxfun(@minus,X,mean(X,2)); %H = eye(nSmp)-ones(nSmp)/nSmp; X = X*H;
St = X*X';

% initialization: using PCA to initialize W
[W,~,~] = svds(St,m,'largest');
% random indexing to initialize G, then get M
WX = W'*X; % equal to W0'*X*H
%startIdx = randi(k,1,nSmp); [vecG,~,~] = kmeansTR(WX,startIdx); 
%G=idxvec2mtx(vecG); M=WX*G*diag(1./diag(G'*G));
vecG = kmeans(WX',k,'MaxIter',200,'Replicates',5,'Start','plus');

% main run
vecG_buff = zeros(1,nSmp);
objV_buff = zeros(1,1);
for iter = 1:iterMax
    % fix G and M to update W
    G = full(sparse(1:nSmp,vecG,1,nSmp,k,nSmp)); 
    Sw = X*(eye(nSmp)-G*diag(1./diag(G'*G))*G')*X';
    [W0,objV_Temp] = eigs(St,Sw,m); %[objV_temp,D_idx]=sort(diag(D_raw),'descend'); W0 = W_raw(:,D_idx(1:m));
    W = W0*mpower(W0'*St*W0,-0.5);
    if any(~isreal(W(:)))
        1; %need more mathmatically analyze
    end
    objV_buff(iter) = sum(diag(objV_Temp));
    % fix W to update G (and M) with Comparison Rule
	vecG_buff(iter,:) = vecG;
    %[vecG,objVal,~] = kmeansTR(W'*X,idxmtx2vec(G));
    X_in = X'*W;
    [vecG,~,objVal] = kmeans(X_in,k,'MaxIter',200,'Replicates',1,'Start','plus');
    objVal = sum(objVal(:));
    for iTime=1:10 %10->4   
        %iIdx = randi(k,1,nSmp); [vecGi,objVi,~] = kmeansTR(W'*X,iIdx);
        [vecGi,~,objVi] = kmeans(X_in,k,'MaxIter',200,'Replicates',1,'Start','plus');
        objVi = sum(objVi(:));
        if objVi<objVal
            vecG = vecGi;
            break;
        end
    end
    vecG = bestMap(vecG_buff(1,:),vecG);
    if vecG_buff(iter,:)==vecG'
        break;
    end
end
% post analysis
%{
if iter > 1
    cntG = histc(vecG_buff,1:k);
    proG = bsxfun(@rdivide,cntG,sum(cntG,1));
    [clsIdx, smpIdx] = find(proG==1);
    vecG = zeros(1,nSmp); vecG(smpIdx)=clsIdx;
    olpIdx = setdiff(1:nSmp,smpIdx);
    proOlp = proG(:,olpIdx); proOlp(end+1,:) = olpIdx; %probability of overlapped spikes
    purX = X(:,smpIdx);
    purX = bsxfun(@minus,purX,mean(purX,2));
    purSt = purX*purX';
    purG = full(sparse(1:nSmp,clsIdx,1,nSmp,k,nSmp)); % purG = idxvec2mtx(clsIdx);    
    purSw = purX*(eye(length(clsIdx))-purG*diag(1./diag(purG'*purG))*purG')*purX';
    [W,~] = eigs(purSt,purSw,m);
end
%}
end

function [postIdx,objVal,M] = kmeansTR(X,preIdx)
%K-means for RTLDA-Km
% X: data matrix, m_variables * n_observations
% preIdx: pre indexes, 1 * n_observations
% 
% postIdx: post indexes, 1 * n_observations
% objVal: objective values
% C: k centres, m_variables * k_centres
[mVar, nObs] = size(X);
k = max(preIdx);
M = zeros(mVar,k);
XX = repmat(sum(X.*X,1),k,1); %each sample stands x'*x
X2 = 2*X;
for iter = 1:200
    for i=1:k %G=idxvec2mtx(preIdx); M=X*G*diag(1./diag(G'*G));
        M(:,i) = mean(X(:,preIdx==i),2);
    end
    CC = repmat(sum(M.*M,1)',1,nObs); %each sample stands c'*c
    difM = abs(XX+CC-M'*X2); %(difference matrix)each sample stands (x-c)'*(x-c)
    [distVec, postIdx] = min( difM,[],1 ); %no need match since it is one time k-means
    if preIdx==postIdx
        break;
    end
    preIdx = postIdx;
end
distSum = zeros(k,1); %each sample stands the distance of within class
for j=1:k
    distTem = distVec(postIdx==j);
    distSum(j) = sum(distTem);
end
objVal = sum(distSum);
end

function idxMtx = idxvec2mtx(idxVec) %E = sparse(1:n,label,1,n,k,n); 
%transform indexes from vector to matrix
% idxVec: indexes vector, 1 * n_observations
% 
% idxMtx: indexes matrix, n_observations * k_classes
k = max(idxVec);
n = numel(idxVec);
idxMtx = zeros(n,k);
for i = 1:k
    idxMtx(idxVec==i,i) = 1;
end
end

function idxVec = idxmtx2vec(idxMtx)
%transform indexes from matrix to vector
% idxMtx: indexes matrix, n_observations * k_classes
% 
% idxVex: indexes vector, 1 * n_observations
[n,k] = size(idxMtx);
idxVec = zeros(1,n);
for i = 1:k
    idxVec(idxMtx(:,i)==1) = i;    
end
end


%{

H = eye(nSmp)-ones(nSmp)/nSmp; 
W = pca((X*H)'); W = W(:,1:m);
objV0 = -1;
objVal = 0;
iter = 1;
epsilon = 1e-1;%1e-10; %eps

while abs(objVal-objV0) > 1e-2
    M1 = W'*X*H;
    M2 = X*H*X'+epsilon*eye(dFea);
    M3 = W'*M2*W;
    T = mpower(M3,-0.5)*M1;
    if iter > 1
        oldG = vecG;
    end
    vecG = kmeans(T',k,'start','plus','EmptyAction','singleton','Replicates',5,'Distance','sqEuclidean','maxiter',200); %k-means
	vecG = idxvec2mtx(vecG);
    if iter > 1
        oldG_vec = idxmtx2vec(oldG);
        G_vec = idxmtx2vec(vecG);
        remap_Gvec = bestMap(oldG_vec,G_vec);
        if all(oldG_vec == remap_Gvec')
            break;
        end
    end
    M4 = vecG*inv(vecG'*vecG)*vecG'*H*X';
    [W,~] = eig(X*H*M4,M2);
    
%     [eigVec,eigVal] = eig(Sb,Sw);
%     [~,inx] = sort(diag(eigVal),'descend');
%     eigVec = eigVec(:,inx);
%     W = eigVec(:,1:d); 
    
    
    objV0 = objVal;
    objVal = sum( diag(M3\M1*M4*W) );
    iter = iter + 1;
end
vecG = idxmtx2vec(vecG);
% gscatter3(Y,H);
end
%}