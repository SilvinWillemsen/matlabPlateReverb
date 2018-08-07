%{
    Dynamic Plate Reverb Plugin
    
    The Dynamic Plate Reverb Plugin is based on the paper 'Plate reverberation:
        Towards the development of a real-time physical model for the working musician'
        by Michele Ducceschi and Craig Webb (http://mdphys.org/PDF/icaPLA_2016.pdf).

    In their implementation it is possible to change the dimensions of the plate,
        but this does not alter the sound in real-time. This plugin extends on their 
        implementation by dynamically altering the sound as the user
        interacts with these parameters. 

    In addition to this, this plugin add the option to move the microphones
    in Lissajous-figures over the plate, creating a flanging-like effect.

%}

classdef realTimePlateReverbPlugin < audioPlugin
    properties
        % Parameters that the user can interact with
        
        % Comment legend: Description, (unit), [minimum value - maximum value]
        
        wetness = 75;       % Dry/Wetness, (%), [0 - 100] (dry signal - only effect)
        Lx = 2;             % Plate Width, (m),  [1 - 3]
        Ly = 1;             % Plate Length, (m), [0.5 - 2]
        cents = 10;          % Amount of cents between eigenmodes (Quality vs. Performance), (Cents), [0.01 - 10]
                      
        physDamp = false;   % Physical Damping [off/on]
        flanging = true;   % Flanging [off/on]
        LFO = false;        % LFO [off/on]
        init = true;        % (Re)Initialise variables
    end
    properties (Access = private)
        
        % Parameters that the user can't interact with
        
        currentSample = 0; % Used for time-variant processes (LFO and flanging)
        circXLength = 0;
        LxSmooth = false;
        LySmooth = false;
        Lxpre = 2;
        Lypre = 1;
        T60 = 4;    %
        rho = 7850; % Material density 
        h = 0.0005; % Plate thickness (m)
        
        stretching = true;  % Plate stretching [off/on]
        lengthOmega = 832; % Number of eigenfrequencies
        omega = zeros (832, 3); % Eigenfrequencies. Col 1: frequency, 
                                                  % Col 2, horizontal modal index,
                                                  % Col 3, vertical modal index
        
        % Coefficients used in the update equation                                     
        coeffBdA = zeros (832, 1);
        coeffCdA = zeros (832, 1);
        coeffIndA = zeros (832, 1)
        
        % Vectors to save the output coefficients for the 
        phiOutL = zeros (832, 1);
        phiOutR = zeros (832, 1);
        
        % Matrix to save all possible output coefficients for the flanging effect
        phiOutMat = zeros (832, 44100 / 4);
        
        % X and Y positions for the dynamic outputs
        circX = zeros (44100 / 4, 1);
        circY = zeros (44100 / 4, 1);

        % QVectors used in the update equation
        qNext = zeros (832, 1);
        qNow =  zeros (832, 1);
        qPrev = zeros (832, 1);
        
        % Save the last sample of a buffer to use in the next buffer
        samp = 0;
        
        p =  [0.4  0.415]; % Input position (x, y), [0 - 1]
        qL = [0.1  0.45];  % Left output position (x, y), [0 - 1]
        qR = [0.84 0.45];  % Right output position (x, y), [0 - 1]
        
        % Mode reduction techniques
        square = false; 
        delModes = true;
        calcCents = true;
        
    end
    properties (Constant)
        
        % Initialise the audioPluginInterface
        
        PluginInterface = audioPluginInterface (...
        audioPluginParameter ('wetness', ...
            'DisplayName', 'Dry/Wet', ...
            'Label', '%', ...
            'Mapping', {'lin', 0, 100}), ...
        audioPluginParameter ('Lx', ...
            'DisplayName', 'Plate Width', ...
            'Label', 'm', ...
            'Mapping', {'lin', 1, 3}), ...
        audioPluginParameter ('Ly', ...
            'DisplayName', 'Plate Height', ...
            'Label', 'm', ...
            'Mapping', {'lin', 0.5, 2}), ...
        audioPluginParameter ('cents', ...
            'DisplayName', 'Cents', ...
            'Label', 'cents', ...
            'Mapping', {'lin', 0.01, 10}), ...
        audioPluginParameter ('flanging', ...
            'DisplayName', 'Flanging', ...
            'Label', 'off/on', ...
            'Mapping', {'enum', 'off', 'on'}), ...
        audioPluginParameter ('LFO', ...
            'DisplayName', 'LFO Stretch', ...
            'Label', 'off/on', ...
            'Mapping', {'enum', 'off', 'on'}), ...
        audioPluginParameter ('init', ...
            'DisplayName', 'Re-initialise', ...
            'Label', 'Click twice', ...
            'Mapping', {'enum', 'off', 'on'}))
%         audioPluginParameter ('stretching', ...
%             'DisplayName', 'Stretching', ...
%             'Label', 'off/on', ...
%             'Mapping', {'enum', 'off', 'on'}), ...
%         audioPluginParameter ('calcCents', ...
%             'DisplayName', 'Calculate Cents', ...
%             'Label', 'off/on', ...
%             'Mapping', {'enum', 'off', 'on'}), ...
%         audioPluginParameter ('cm', ...
%             'DisplayName', 'Physical Damping', ...
%             'Label', 'off/on', ...
%             'Mapping', {'enum', 'off', 'on'}), ...
    end
    methods
        function plugin = realTimePlateReverbPlugin        %<---
            
        end   
        function out = process (plugin, in)
            % Initialise the plugin
%             tic
            if plugin.init == true
                
                % Create an vector including settings and in- and outputs
                settings = [plugin.delModes plugin.square plugin.calcCents ...
                           plugin.flanging plugin.physDamp];
                inOutputs = [plugin.p; plugin.qL; plugin.qR];
                disp ('Initialising Plugin')
                
                % Obtain initial coefficients using the variables initialised above
                % Note: always initialise with Lx = 2 and Ly = 1 so that flanging will work correctly
                [plugin.coeffBdA, plugin.coeffCdA, plugin.coeffIndA, ...
                 plugin.omega, plugin.phiOutL, plugin.phiOutR, ...
                 plugin.phiOutMat, plugin.circXLength, plugin.rho, plugin.T60]...
                    = initPlate (2, 1, plugin.cents, ...
                                 inOutputs, settings);
                
                plugin.qNext = zeros (length (plugin.omega (:, 1)), 1);
                plugin.qNow = zeros (length (plugin.omega (:, 1)), 1);
                plugin.qPrev = zeros (length (plugin.omega (:, 1)), 1);
                
                % First sample is 0
                plugin.samp = 0;
                
                plugin.lengthOmega = length (plugin.omega (:, 1));
                disp (plugin.lengthOmega)
                
                % Set the init variable to false
                plugin.init = false;
            end
            
            % Initialise the stereo out-vector
            out = zeros (length (in), 2);
            
            % Initialise new variables 
            qNextPre =         plugin.qNext;
            qNowPre =          plugin.qNow;
            qPrevPre =         plugin.qPrev;
            coeffBdAPre =      plugin.coeffBdA;
            coeffCdAPre =      plugin.coeffCdA;
            coeffIndAPre =     plugin.coeffIndA;
            phiOutLPre =       plugin.phiOutL;
            phiOutRPre =       plugin.phiOutR;
            
            M = plugin.lengthOmega;
            pLoop = plugin.p;
            qLLoop = plugin.qL;
            qRLoop = plugin.qR;
            kSquared = (2e11 * plugin.h^2) / (12 * plugin.rho * (1 - 0.3^2));
            k = 1 / 44100;
            LxSmoothUse = plugin.Lxpre;
            LySmoothUse = plugin.Lypre;
            
            sS = round (plugin.lengthOmega / 20);

            if abs (plugin.Lx - plugin.Lxpre) > 1 / sS
                if round (LxSmoothUse * sS) / sS > round (plugin.Lx * sS) / sS
                    LxSmoothUse = LxSmoothUse - 1 / sS;
                else
                    if round (LxSmoothUse * sS) / sS < round (plugin.Lx * sS) / sS
                        LxSmoothUse = LxSmoothUse + 1 / sS;
                    end
                end
            else
                LxSmoothUse = plugin.Lx;
            end

            if abs (plugin.Ly - plugin.Lypre) > 1 / sS
                if round (LySmoothUse * sS) / sS > round (plugin.Ly * sS) / sS
                    LySmoothUse = LySmoothUse - 1 / sS;
                else
                    if round (LySmoothUse * sS) / sS < round (plugin.Ly * sS) / sS
                        LySmoothUse = LySmoothUse + 1 / sS;
                    end
                end
            else
                LySmoothUse = plugin.Ly;
            end
            plugin.Lxpre = LxSmoothUse;
            plugin.Lypre = LySmoothUse;

            if plugin.LFO == true
                LxLoop = LxSmoothUse + (sin (2 * pi * plugin.currentSample / 44100) / 4);
            else
                LxLoop = LxSmoothUse;
            end

            LyLoop = LySmoothUse;
            omegaLoop = plugin.omega;

            % Disable stretching if flanging is true (done for speed)
            if plugin.flanging == true
                LxLoop = 2;
                LyLoop = 1;
            end
            
            phiOutLPre = zeros (M, 1);
            phiOutRPre = zeros (M, 1);
            coeffAPre = (1 / k^2) + ((12 .* (log(10) ./ 4)) / (plugin.rho * plugin.h * k));
            coeffAPreAll = coeffAPre * plugin.rho * plugin.h;

            % Update all coefficients that depend on the stretching of the plate
            omegaLoop (:, 1)= ((omegaLoop (:, 2) / LxLoop).^2 + (omegaLoop (:, 3) / LyLoop).^2) * sqrt (kSquared) * pi^2;
            coeffBdAPre (:, 1) = ((2 / k^2) - (omegaLoop (:, 1)).^2) ./ coeffAPre;
            coeffIndAPre (:, 1) = ((4 / (LxLoop * LyLoop)) * sin (omegaLoop (:, 2) * pi * pLoop (1)) .* sin (omegaLoop (:, 3) * pi * pLoop (2))) ./ coeffAPreAll;
            phiOutLPre (:, 1) = (4 / (LxLoop * LyLoop)) * sin (omegaLoop (:, 2) * pi * qLLoop (1)) .* sin (omegaLoop (:, 3) * pi * qLLoop (2));
            phiOutRPre (:, 1) = (4 / (LxLoop * LyLoop)) * sin (omegaLoop (:, 2) * pi * qRLoop (1)) .* sin (omegaLoop (:, 3) * pi * qRLoop (2));


            % Make sure that eigenfrequencies that are higher than the stability condition (2*fs) are ignored
            index = zeros (M, 1);
            i = 1;
            for m = 1 : M
                if omegaLoop (m, 1) < 44100 * 2 
                    index (i) = m;
                    i = i + 1;
                end
            end
            index = index (1 : i - 1);
            plugin.omega = omegaLoop;

            % Update the loop variables
            qNextLoop =      qNextPre (index);
            qNowLoop =       qNowPre (index);
            qPrevLoop =      qPrevPre (index);
            coeffBdALoop =   coeffBdAPre (index);
            coeffCdALoop =   coeffCdAPre (index);
            coeffIndALoop =  coeffIndAPre (index);
            phiOutLLoop =    phiOutLPre (index);
            phiOutRLoop =    phiOutRPre (index);

%             if plugin.flanging == true
                curSamp = plugin.currentSample;
                lengthCircX = plugin.circXLength;
%                 phiOutFlange = plugin.phiOutMat; 
                
%             end
            
%             divX =  omegaLoop (:, 2) * pi / LxLoop;                    
%             divY =  omegaLoop (:, 3) * pi / LyLoop;
%             scalar = (4 / (LxLoop * LyLoop));
%             idx = floor (mod ((curSamp + (1:length(in))) / 2, lengthCircX) + 1);
%             phiOutLoopPre = scalar *  sin (circXLoop (idx) .* divX) .* sin (circYLoop (idx) .* divY);
            
            % Update Equation
            for t = 1 : length(in)
                if t == 1
                    qNextLoop = (coeffBdALoop .* qNowLoop + coeffCdALoop .* qPrevLoop + coeffIndALoop .* plugin.samp);
                else
                    qNextLoop = (coeffBdALoop .* qNowLoop + coeffCdALoop .* qPrevLoop + coeffIndALoop .* in (t - 1, 1));
                end

                
                % If the flanging is turned on, find the correct vector from the phiOutPre matrices              
                if plugin.flanging == true
                    phiOutLLoop = plugin.phiOutMat (:, floor (mod ((curSamp + t) / 128, lengthCircX) + 1));
                    phiOutRLoop = plugin.phiOutMat (:, floor (mod ((curSamp + t) / 64, lengthCircX) + 1));
% %                     phiOutLLoop = phiOutLoopPre (:, t);
% %                     phiOutRLoop = phiOutLLoop;
                end
                
                out(t, 1) = plugin.wetness / 100 * 25000 * sum (qNextLoop .* phiOutLLoop) ...
                    + (1 - plugin.wetness / 100) * in (t, 1);
                out(t, 2) = plugin.wetness / 100 * 25000 * sum (qNextLoop .* phiOutRLoop) ...
                    + (1 - plugin.wetness / 100) * in (t, 1);

                qPrevLoop = qNowLoop;
                qNowLoop = qNextLoop;
            end
            
            % Set the first sample in the next buffer to the last sample in the current buffer
            plugin.samp = in (end);
            
            % Save the QVectors for the next loop
            plugin.qNext (index) =  qNextLoop;
            plugin.qPrev (index) =  qPrevLoop;
            plugin.qNow (index) =   qNowLoop;
            
            % Update the currentSample
            plugin.currentSample = plugin.currentSample + length (in (:, 1));
%             toc
        end
        
        function reset (plugin)
        end
        
        
    end
end