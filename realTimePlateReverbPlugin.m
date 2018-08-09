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

    Additional instructions: 
    If many samples are underrun, close the audioTestBench and 'clear all'.
%}

classdef realTimePlateReverbPlugin < audioPlugin
    properties
        % Parameters that the user can interact with
        
        % Comment legend: Description, (unit), [minimum value - maximum value]
        
        wetness = 75;       % Dry/Wetness, (%), [0 - 100] (dry signal - only effect)
        Lx = 2;             % Plate Width, (m),  [1 - 3]
        Ly = 1;             % Plate Length, (m), [0.5 - 2]
        cents = 5;          % Amount of cents between eigenmodes (Quality vs. Performance), (Cents), [0.01 - 10]
           
        speed = 0.33;       % Speed of the flanging (Revolutions^-3 / s), [0.1 - 1]
        physDamp = false;   % Physical Damping [off/on]
        flanging = true;    % Flanging [off/on]
        init = false;       % (Re)Initialise variables
    end
    properties (Access = private)
        
        % Parameters that the user can't interact with
        
        circXLength = 0;
        curLx = 2;
        curLy = 1;
        T60 = 4;    %
        rho = 7850; % Material density 
        h = 0.0005; % Plate thickness (m)
        
        stretching = true;  % Plate stretching [off/on]
        lengthOmega = 2; % Number of eigenfrequencies
        omega = zeros (2, 3); % Eigenfrequencies. Col 1: frequency, 
                                                  % Col 2, horizontal modal index,
                                                  % Col 3, vertical modal index
        
        % Coefficients used in the update equation                                     
        coeffBdA = zeros (2, 1);
        coeffCdA = zeros (2, 1);
        coeffIndA = zeros (2, 1)
        
        % Vectors to save the output coefficients for the 
        phiOutL = zeros (2, 1);
        phiOutR = zeros (2, 1);
        
        % Matrix to save all possible output coefficients for the flanging effect
        flangeMatSize = 44100;
        phiOutMat = zeros (1206, 44100);
        
        % X and Y positions for the dynamic outputs
        circX = zeros (44100, 1);
        circY = zeros (44100, 1);

        % QVectors used in the update equation
        qNext = zeros (2, 1);
        qNow =  zeros (2, 1);
        qPrev = zeros (2, 1);
        
        % Save the last sample of a buffer to use in the next buffer
        samp = 0;
        
        p =  [0.4  0.415]; % Input position (x, y), [0 - 1]
        qL = [0.1  0.45];  % Left output position (x, y), [0 - 1]
        qR = [0.84 0.45];  % Right output position (x, y), [0 - 1]
        
        % Mode reduction techniques
        square = false; 
        delModes = true;
        calcCents = true;
        
        % Scaling the output
        scaling = (2 * 1)^-0.25 * 150000;
        
        flangeIdx = 0;
        
        index = zeros (2, 1);
        
        wetnessLoop = 75;
        
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
            'DisplayName', 'Quality', ...
            'Mapping', {'lin', 0.01, 10}), ...
        audioPluginParameter ('flanging', ...
            'DisplayName', 'Flanging', ...
            'Label', 'off/on', ...
            'Mapping', {'enum', 'off', 'on'}), ...
        audioPluginParameter ('speed', ...
            'DisplayName', 'Flanging Speed', ...
            'Label', 'Rev / s', ...
            'Mapping', {'lin', 0.1, 1}), ...
        audioPluginParameter ('init', ...
            'DisplayName', 'Re-initialise', ...
            'Label', 'Click twice', ...
            'Mapping', {'enum', 'off', 'on'}))
    end
    methods
        function plugin = realTimePlateReverbPlugin
            initialise (plugin);
        end 
        
        function initialise (plugin)
            % Create an vector including settings and in- and outputs
                inOutputs = [plugin.p; plugin.qL; plugin.qR];
                disp ('Initialising Plugin')
                
                % Obtain initial coefficients
                % Note: always initialise with Lx = 2 and Ly = 1 so that flanging will work correctly.
                % If plugin.Lx or plugin.Ly is not 2 or 1 respectively, update loopvectors  
                [plugin.coeffBdA, plugin.coeffCdA, plugin.coeffIndA, ...
                 plugin.omega, plugin.phiOutL, plugin.phiOutR, ...
                 plugin.phiOutMat, plugin.circXLength, plugin.rho, plugin.T60, ...
                 plugin.circX, plugin.circY]...
                    = initPlate (2, 1, 10.01 - plugin.cents, plugin.flangeMatSize, inOutputs);
                plugin.qNext = zeros (length (plugin.omega (:, 1)), 1);
                plugin.qNow = zeros (length (plugin.omega (:, 1)), 1);
                plugin.qPrev = zeros (length (plugin.omega (:, 1)), 1);
                
                % First sample is 0
                plugin.samp = 0;
                
                plugin.lengthOmega = length (plugin.omega (:, 1));
                disp (plugin.lengthOmega)
                
                % If the sliders are not what the coefficients were initialised with (Lx = 2, Ly = 1),
                % wait with setting the init variable to false as the loopvectors need to be updated 
                if (plugin.Lx == 2 && plugin.Ly == 1)
                    plugin.init = false;
                end
                
                % Set the init variable to false
                plugin.index = (1:plugin.lengthOmega)';
                
                plugin.wetnessLoop = plugin.wetness;
                
        end
        function out = process (plugin, in)
            % Initialise the plugin
%             tic
            if plugin.init == true
                initialise (plugin);
            end
            
            % Initialise the stereo out-vector
            out = zeros (length (in), 2);
          
            % If there is no change between the width & length sliders and the current width and length AND
            % the init variable is false, retrieve the loopvectors directly from the plugin.
            if plugin.Lx == plugin.curLx && plugin.Ly == plugin.curLy && plugin.init == false
                qNextLoop =         plugin.qNext (plugin.index);
                qNowLoop =          plugin.qNow (plugin.index);
                qPrevLoop =         plugin.qPrev (plugin.index);
                coeffBdALoop =      plugin.coeffBdA (plugin.index);
                coeffCdALoop =      plugin.coeffCdA (plugin.index);
                coeffIndALoop =     plugin.coeffIndA (plugin.index);
                phiOutLLoop =       plugin.phiOutL (plugin.index);
                phiOutRLoop =       plugin.phiOutR (plugin.index);
            else
            % Otherwise recalculate the loopvariables
                % If the init variable is still true, this means that the plugin was re-initialised with Lx ~= 2 and Ly ~= 1
                if plugin.init == true
                    plugin.init = false;
                end
                % Copy the entire vectors
                qNextPre =         plugin.qNext;
                qNowPre =          plugin.qNow;
                qPrevPre =         plugin.qPrev;
                coeffBdAPre =      plugin.coeffBdA;
                coeffCdAPre =      plugin.coeffCdA;
                coeffIndAPre =     plugin.coeffIndA;

                M = plugin.lengthOmega;
                pLoop = plugin.p;
                qLLoop = plugin.qL;
                qRLoop = plugin.qR;
                kSquared = (2e11 * plugin.h^2) / (12 * plugin.rho * (1 - 0.3^2));
                k = 1 / 44100;

                sS = round (plugin.lengthOmega / 20);
                if round (plugin.curLx * sS) / sS > round (plugin.Lx * sS) / sS
                    plugin.curLx = plugin.curLx - 1 / sS;
                else
                    if round (plugin.curLx * sS) / sS < round (plugin.Lx * sS) / sS
                        plugin.curLx = plugin.curLx + 1 / sS;
                    else
                        plugin.curLx = plugin.Lx;
                    end
                end

                if round (plugin.curLy * sS) / sS > round (plugin.Ly * sS) / sS
                    plugin.curLy = plugin.curLy - 1 / sS;
                else
                    if round (plugin.curLy * sS) / sS < round (plugin.Ly * sS) / sS
                        plugin.curLy = plugin.curLy + 1 / sS;
                    else
                        plugin.curLy = plugin.Ly;
                    end
                end
                LxLoop = plugin.curLx;
                LyLoop = plugin.curLy;
                omegaLoop = plugin.omega;

                phiOutLPre = zeros (M, 1);
                phiOutRPre = zeros (M, 1);
                coeffAPre = (1 / k^2) + ((12 .* (log(10) ./ 4)) / (plugin.rho * plugin.h * k));
                coeffAPreAll = coeffAPre * plugin.rho * plugin.h;

                % Update all coefficients that depend on the stretching of the plate
                omegaLoop (:, 1)= ((omegaLoop (:, 2) / LxLoop).^2 + (omegaLoop (:, 3) / LyLoop).^2) * sqrt (kSquared) * pi^2;
                coeffBdAPre (:, 1) = ((2 / k^2) - (omegaLoop (:, 1)).^2) ./ coeffAPre;
                coeffIndAPre (:, 1) = (sin (omegaLoop (:, 2) * pi * pLoop (1)) .* sin (omegaLoop (:, 3) * pi * pLoop (2))) ./ coeffAPreAll;
                phiOutLPre (:, 1) = sin (omegaLoop (:, 2) * pi * qLLoop (1)) .* sin (omegaLoop (:, 3) * pi * qLLoop (2));
                phiOutRPre (:, 1) = sin (omegaLoop (:, 2) * pi * qRLoop (1)) .* sin (omegaLoop (:, 3) * pi * qRLoop (2));


                % Make sure that eigenfrequencies that are higher than the stability condition (2*fs) are ignored
                plugin.index = zeros (M, 1);
                i = 1;
                for m = 1 : M
                    if omegaLoop (m, 1) < 44100 * 2 
                        plugin.index (i) = m;
                        i = i + 1;
                    end
                end
                plugin.index = plugin.index (1 : i - 1);
                plugin.omega = omegaLoop;

                % Update the loopvectors
                qNextLoop =      qNextPre (plugin.index);
                qNowLoop =       qNowPre (plugin.index);
                qPrevLoop =      qPrevPre (plugin.index);
                coeffBdALoop =   coeffBdAPre (plugin.index);
                coeffCdALoop =   coeffCdAPre (plugin.index);
                coeffIndALoop =  coeffIndAPre (plugin.index);
                phiOutLLoop =    phiOutLPre (plugin.index);
                phiOutRLoop =    phiOutRPre (plugin.index);
            end
            
            plugin.scaling = (plugin.curLx * plugin.curLy)^-0.25 * 150000; 

            % The value by which the flangeIndex will increment. Depends on the speed variable.   
%             flangeIncVal = (plugin.flangeMatSize / 44100) / (128 * (1 - log10 (plugin.speed * 10) + (1 / 64)));
            flangeIncVal = (plugin.flangeMatSize / 44100) * plugin.speed^3;
            
            % if the difference between the previous wetness value and
            % the current value, smooth out 
            if round (plugin.wetnessLoop - plugin.wetness) > 1
                plugin.wetnessLoop = plugin.wetnessLoop - 1;
            else
                if round (plugin.wetnessLoop - plugin.wetness) < -1
                    plugin.wetnessLoop = plugin.wetnessLoop + 1;
                else 
                    plugin.wetnessLoop = plugin.wetness;
                end
            end
            
            % Update Equation
            for t = 1 : length(in)
                if t == 1
                    qNextLoop = (coeffBdALoop .* qNowLoop + coeffCdALoop .* qPrevLoop + coeffIndALoop .* plugin.samp);
                else
                    qNextLoop = (coeffBdALoop .* qNowLoop + coeffCdALoop .* qPrevLoop + coeffIndALoop .* in (t - 1, 1));
                end

                % If the flanging is turned on, find the correct vector from the phiOutPre matrices              
                if plugin.flanging == true
                    plugin.flangeIdx = plugin.flangeIdx + flangeIncVal;
                    % Let the left microphone move twice as fast as the right microphone
                    phiOutLLoop = plugin.phiOutMat (plugin.index, floor (mod (plugin.flangeIdx, plugin.circXLength) + 1));
                    phiOutRLoop = plugin.phiOutMat (plugin.index, floor (mod (plugin.flangeIdx / 2, plugin.circXLength) + 1));
                end
                
                out(t, 1) = plugin.wetnessLoop / 100 * plugin.scaling * sum (qNextLoop .* phiOutLLoop) ...
                    + (1 - plugin.wetnessLoop / 100) * in (t, 1);
                out(t, 2) = plugin.wetnessLoop / 100 * plugin.scaling * sum (qNextLoop .* phiOutRLoop) ...
                    + (1 - plugin.wetnessLoop / 100) * in (t, 1);

                qPrevLoop = qNowLoop;
                qNowLoop = qNextLoop;
            end

            % Set the first sample in the next buffer to the last sample in the current buffer
            plugin.samp = in (end);
            
            % Save the vectors for the next loop
            plugin.qNext (plugin.index) =       qNextLoop;
            plugin.qPrev (plugin.index) =       qPrevLoop;
            plugin.qNow (plugin.index) =        qNowLoop;
            plugin.coeffBdA (plugin.index) =    coeffBdALoop;
            plugin.coeffCdA (plugin.index) =    coeffCdALoop;
            plugin.coeffIndA (plugin.index) =   coeffIndALoop;
            plugin.phiOutL (plugin.index) =     phiOutLLoop;
            plugin.phiOutR (plugin.index) =     phiOutRLoop;

        end
        
        function reset (plugin)
        end
    end
end