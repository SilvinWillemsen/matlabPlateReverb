frameLength = 256;
deviceReader = audioDeviceReader('SamplesPerFrame', frameLength);
sampleRate = fileReader.SampleRate;
hostedPlugin = loadAudioPlugin('realTimePlateReverbPlugin.vst');
deviceWriter = audioDeviceWriter('SampleRate',sampleRate);
setSampleRate(hostedPlugin,sampleRate);
tic
while toc<25
    audioIn = deviceReader();

    % The hosted plugin requires a stereo input.
    stereoAudioIn = [audioIn,audioIn];

    x = process(hostedPlugin,stereoAudioIn);

    deviceWriter(x);
end

release(fileReader)
release(deviceWriter)