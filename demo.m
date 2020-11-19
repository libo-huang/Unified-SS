%% fix k=3
clc; clear;
% Quiroga, R. Quian
aveTim = [];
for i_data = 13
    [spikeRaw,spikeClean,spikeOutlier] = loadSpks_QQ(i_data);
    Data = spikeRaw;
    X = Data.spikes;
    cntAcc=[]; cntTim=[];
    for i_time = 1:20
        tic;
        reaIdx = Data.spike_class;
        [estIdx,W] = PCAKm(X,3);
        remIdx = bestMap(estIdx, reaIdx);
        evaStats = confusionmatStats(estIdx, remIdx);
        cntAcc = [cntAcc, evaStats.TOTAL_ACC];
        cntTim = [cntTim, toc];
    end
	aveTim = [aveTim, mean(cntTim)];
    fprintf('\t PCAKm: \t %4.2f±%4.2f \t %4.2f±%4.2f \n',...
      mean(cntAcc)*100,std(cntAcc)*100, mean(cntTim),std(cntTim));
end
fprintf('\t average time: \t %4.2f±%4.2f \n',mean(aveTim),std(aveTim));  
  
%% estimate class number with 4 methods
clc;clear;
% Quiroga, R. Quian
for i_data = 13 %1:20
    [spikeRaw,spikeClean,spikeOutlier] = loadSpks_QQ(i_data);
    Data = spikeClean;
    X = Data.spikes;
    reaIdx = Data.spike_class;
    [estIdx,W] = PCAKm_auto(X,3);
    if max(estIdx)==3
        remIdx = bestMap(estIdx, reaIdx);
        evaSta = confusionmatStats(estIdx, remIdx);
        fprintf('\t PCAKm_auto: \t %d \t %6.2f \n',...
            max(estIdx), evaSta.TOTAL_ACC*100);
    else
        fprintf('\t PCAKm_auto: \t %d \n',max(estIdx));
    end
end