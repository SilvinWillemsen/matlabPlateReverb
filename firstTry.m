clear all;
close all;
clc;
frameLength = 512;
deviceReader = audioDeviceReader;
deviceWriter = audioDeviceWriter('SampleRate',deviceReader.SampleRate);
scope = dsp.TimeScope(...
    'SampleRate',deviceReader.SampleRate,...
    'TimeSpan',2,...
    'BufferLength',2*deviceReader.SampleRate,...
    'YLimits',[-1 1]);

disp('Begin Signal Input...')
output = [];
samp = 0;
[factorBdA, factorCdA, factorIndA, omega, phiOutL, phiOutR] = initPlate(2,1,0,0);

qNext = zeros(length(omega(:,1)),1);
qPre = zeros(length(omega(:,1)),1);
qNow = zeros(length(omega(:,1)),1);
qPrev = zeros(length(omega(:,1)),1);

tic
while toc < 10
    mySignal = deviceReader();
    output = zeros(1,length(mySignal));
    for t = 1:length(mySignal)
        if t == 1
            qNext =(factorBdA.*qNow+factorCdA.*qPrev+factorIndA.*samp);
        else
            qNext =(factorBdA.*qNow+factorCdA.*qPrev+factorIndA.*mySignal(t-1));
        end
        
        output(1,t) = 100000*qNext'*phiOutL;
        output(2,t) = 100000*qNext'*phiOutR;
        qPrev = qNow;
        qNow = qNext;
    end
    samp = mySignal(end);
    deviceWriter(output');
    scope([mySignal output']);
end
disp('End Signal Input') 
release(deviceReader);
release(deviceWriter);
%release(scope);