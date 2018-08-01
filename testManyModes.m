asavePluginStuff = load('savePluginStuff.mat');
savePluginStuff = savePluginStuff.savePluginStuff;
i = 1;
prevModes = savePluginStuff(1,2);
meanToc = zeros(1,10);
meanTocDiv = 1;
meanNext = 1;
modes = zeros(1,10);
while i < length(savePluginStuff(:,1))
    if savePluginStuff(i,1) > 0.5 || savePluginStuff(i,1) == 0
        savePluginStuff(i,:) = [];
    else
        if prevModes == savePluginStuff(i,2)
            meanToc(1,meanNext) = meanToc(1,meanNext) + savePluginStuff(i,1);
            meanTocDiv = meanTocDiv + 1;
        else
            modes(1,meanNext) = prevModes;
            meanToc(1,meanNext) = meanToc(1,meanNext)./meanTocDiv;
            meanNext = meanNext + 1;
            meanTocDiv = 1;
        end
        prevModes = savePluginStuff(i,2);
        i = i+1;
    end
end
 modes(1,meanNext) = prevModes;
meanToc(1,meanNext) = meanToc(1,meanNext)./meanTocDiv;
meanNext = meanNext + 1;
meanTocDiv = 1;
i = 1;
while i < length(meanToc)
    if meanToc == 0
        meanToc(i) = [];
        modes(i) = [];
    else
        i = i + 1;
    end
end
corr(meanToc', modes')