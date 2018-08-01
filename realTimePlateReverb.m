clear all;
close all;
clc;

frameLength = 256;
deviceReader = audioDeviceReader('SamplesPerFrame', frameLength);
deviceWriter = audioDeviceWriter('SampleRate',deviceReader.SampleRate);
scope = dsp.TimeScope(...
    'SampleRate',deviceReader.SampleRate,...
    'TimeSpan',1,...
    'BufferLength',deviceReader.SampleRate,...
    'YLimits',[-1 1]);

disp('Begin Signal Input...')
output = [];
samp = 0;
options = [1 1 1 0 0.5]; %[delModes calcCent phasing stretching cents]
phasing = 0;
stretching = 0;
calcCent = true;
[factorBdA, factorCdA, factorIndA, omega, phiOutL, ...
    phiOutR, phiOutLPre, phiOutRPre, circXLength]...
    = initPlate(2,1,options);

Lspeed = 50;
Rspeed = 30;
qNext = zeros(length(omega(:,1)),1);
qPre = zeros(length(omega(:,1)),1);
qNow = zeros(length(omega(:,1)),1);
qPrev = zeros(length(omega(:,1)),1);
%%
tic
while toc < 5
    
    mySignal = deviceReader();
    output = zeros(1,length(mySignal));
    for t = 1:length(mySignal)
        if t == 1
            qNext =(factorBdA.*qNow+factorCdA.*qPrev+factorIndA.*samp);
        else
            qNext =(factorBdA.*qNow+factorCdA.*qPrev+factorIndA.*mySignal(t-1));
        end
        if phasing == true
            output(1,t) = 50000*qNext'*phiOutLPre(:,floor(mod(t/Lspeed,circXLength)+1));
            output(2,t) = 50000*qNext'*phiOutRPre(:,floor(mod(t/Rspeed,circXLength)+1));
        else
            output(1,t) = 50000*qNext'*phiOutL;
            output(2,t) = 50000*qNext'*phiOutR;
        end
        qPrev = qNow;
        qNow = qNext;
    end
    samp = mySignal(end);
    deviceWriter(output');
    
    %scope([mySignal output']);
end
disp('End Signal Input') 
release(deviceReader);
release(deviceWriter);
%release(scope);