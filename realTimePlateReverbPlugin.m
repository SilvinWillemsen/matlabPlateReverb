classdef realTimePlateReverbPlugin < audioPlugin
    properties
        
        % Comment legend: Description, (unit), [minimum value - maximum value]
        
        wetness = 50;       % Dry/Wetness, (%), [0 - 100] (dry signal - only effect)
        Lx = 2;             % Plate Width, (m),  [1 - 3]
        Ly = 1;             % Plate Length, (m), [0.5 - 2]
        cents = 4;          % Amount of cents between eigenmodes, Cents, [0.01 - 10]
        
        cm = false;         % Physical Damping [off/on]
        flanging = false;   % Flanging [off/on]
        square = true;      % 
        stretching = true;  % Plate stretching [off/on]
        LFO = false;        % LFO [off/on]
        init = true;        % (Re)Initialise variables
    end
    properties (Access = private)
        
        currentSample = 0;
        circXLength = 0;
        LxSmooth = 2;
        LySmooth = 1;
        Lxpre = 2;
        Lypre = 1;
        T60 = 4;
        rho = 7850;
        h = 0.0005;
        lengthOmega = 832;
        
        omega = zeros (832, 3);
        factorBdA = zeros (832, 1);
        factorCdA = zeros (832, 1);
        factorIndA = zeros (832, 1);
        phiOutL = zeros (832, 1);
        phiOutR = zeros (832, 1);
        phiOutLPre = zeros (832, 22050);
        phiOutRPre = zeros (832, 22050);
        qNext = zeros (832, 1);
        qNow = zeros (832, 1);
        qPrev = zeros (832, 1);
        qPre = zeros (832, 1);
        samp = 0;
        prevLengthOmega = 0;
        p = [0.4 0.415];
        qL = [0.1 0.45];
        qR = [0.84 0.45];
        saveMat = zeros (10000,2);
        initNum = 1;
        
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
        audioPluginParameter ('stretching', ...
            'DisplayName', 'Stretching', ...
            'Label', 'off/on', ...
            'Mapping', {'enum', 'off', 'on'}), ...
        audioPluginParameter ('LFO', ...
            'DisplayName', 'LFO Stretch', ...
            'Label', 'off/on', ...
            'Mapping', {'enum', 'off', 'on'}), ...
        audioPluginParameter ('cm', ...
            'DisplayName', 'Physical Damping', ...
            'Label', 'off/on', ...
            'Mapping', {'enum', 'off', 'on'}), ...
        audioPluginParameter ('init', ...
            'DisplayName', 'Re-initialise', ...
            'Label', 'Click twice', ...
            'Mapping', {'enum', 'off', 'on'}))
   
%         audioPluginParameter ('calcCents', ...
%             'DisplayName', 'Calculate Cents', ...
%             'Label', 'off/on', ...
%             'Mapping', {'enum', 'off', 'on'}), ...
    end
    methods
        function plugin = realTimePlateReverbPlugin        %<---
            
        end   
        function out = process (plugin, in)
            
            % Initialise the plugin
            if plugin.init == true
                options = [plugin.delModes plugin.square plugin.calcCents plugin.flanging plugin.stretching plugin.cents plugin.cm]; %[delModes calcCent flanging stretching cents]
                disp ('Initialising Plugin')
                
                [plugin.factorBdA, plugin.factorCdA, plugin.factorIndA, ...
                    plugin.omega, plugin.phiOutL, plugin.phiOutR, plugin.phiOutLPre, plugin.phiOutRPre, ...
                    plugin.p, plugin.qL, plugin.qR, plugin.circXLength]...
                    = initPlate (plugin.Lx, plugin.Ly, options);
                
                plugin.qNext = zeros (length (plugin.omega (:, 1)), 1);
                plugin.qPre = zeros (length (plugin.omega (:, 1)), 1);
                plugin.qNow = zeros (length (plugin.omega (:, 1)), 1);
                plugin.qPrev = zeros (length (plugin.omega (:, 1)), 1);
                plugin.samp = 0;
                
                plugin.lengthOmega = length (plugin.omega (:, 1));
                plugin.prevLengthOmega = length (plugin.omega (:, 1));
                plugin.initNum = plugin.initNum + 1;
                disp (plugin.lengthOmega)
                
                plugin.init = false;
            end
            out = zeros (length (in), 2);
            if plugin.prevLengthOmega > plugin.lengthOmega 
                qNextLoop =         plugin.qNext (1 : plugin.lengthOmega);
                qNowLoop =          plugin.qNow (1 : plugin.lengthOmega);
                qPrevLoop =         plugin.qPrev (1 : plugin.lengthOmega);
                factorBdALoop =     plugin.factorBdA (1 : plugin.lengthOmega);
                factorCdALoop =     plugin.factorCdA (1 : plugin.lengthOmega);
                factorIndALoop =    plugin.factorIndA (1 : plugin.lengthOmega);
                phiOutLLoop =       plugin.phiOutL (1 : plugin.lengthOmega);
                phiOutRLoop =       plugin.phiOutR (1 : plugin.lengthOmega);
            else
                if plugin.prevLengthOmega < plugin.lengthOmega
                    qNextLoop =         [plugin.qNext; zeros(plugin.lengthOmega - length (plugin.qNext), 1)];
                    qNowLoop =          [plugin.qNow; zeros(plugin.lengthOmega - length (plugin.qNow), 1)];
                    qPrevLoop =         [plugin.qPrev; zeros(plugin.lengthOmega - length (plugin.qPrev), 1)];
                    factorBdALoop =     [plugin.factorBdA; zeros(plugin.lengthOmega - length (plugin.factorBdA), 1)];
                    factorCdALoop =     [plugin.factorCdA;zeros(plugin.lengthOmega - length (plugin.factorCdA), 1)];
                    factorIndALoop =    [plugin.factorIndA;zeros(plugin.lengthOmega - length (plugin.factorIndA), 1)];
                    phiOutLLoop =       [plugin.phiOutL;zeros(plugin.lengthOmega - length (plugin.phiOutL), 1)];
                    phiOutRLoop =       [plugin.phiOutR;zeros(plugin.lengthOmega - length (plugin.phiOutR), 1)];
                else
                    qNextLoop =         plugin.qNext;
                    qNowLoop =          plugin.qNow;
                    qPrevLoop =         plugin.qPrev;
                    factorBdALoop =     plugin.factorBdA;
                    factorCdALoop =     plugin.factorCdA;
                    factorIndALoop =    plugin.factorIndA;
                    phiOutLLoop =       plugin.phiOutL;
                    phiOutRLoop =       plugin.phiOutR;
                end
            end
            
            if plugin.stretching == true
                M = plugin.lengthOmega;
                pLoop = plugin.p;
                qLLoop = plugin.qL;
                qRLoop = plugin.qR;
                kSquared = 0.5833;
                k = 1 / 44100;
                LxSmoothUse = plugin.LxSmooth;
                LySmoothUse = plugin.LySmooth;
                sS = round (M / 20);
                
                if abs (plugin.Lx - plugin.Lxpre) > 1 / sS
                    plugin.smoothLx = true;
                else
                    plugin.smoothLx = false;
                end
                
                if abs (plugin.Ly - plugin.Lypre) > 1 / sS
                    plugin.smoothLy = true;
                else
                    plugin.smoothLy = false;
                end
                
                if plugin.smoothLx == true
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
                
                if plugin.smoothLy == true
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
%                 disp (LxSmoothUse);
                
                if plugin.LFO == true
                    LxLoop = LxSmoothUse + (sin (2 * pi * plugin.currentSample / 44100) / 4);
                else
                    LxLoop = LxSmoothUse;
                end
                
                plugin.LxSmooth = LxSmoothUse;
                plugin.LySmooth = LySmoothUse;
                LyLoop = plugin.LySmooth;
                omegaLoop = plugin.omega;
                i = 0;
                phiOutLLoop = zeros (M, 1);
                phiOutRLoop = zeros (M, 1);
                factorALoop = (1 / k^2) + ((12 .* (log(10) ./ plugin.T60)) / (7850 * 0.0005 * k));
                rhoUse = plugin.rho;
                hUse = plugin.h;
                factorALoopAll = factorALoop * rhoUse * hUse;
                
                for m = 1 : M
                    omegaLoop (m, 1)= (((omegaLoop(m, 2) * pi) / LxLoop)^2 + ((omegaLoop (m, 3) * pi) / LyLoop)^2) * sqrt (kSquared);
                    
                    if omegaLoop(m, 1) < 44100 * 2
                        i = i + 1;
                    end
                    factorBdALoop (m, 1) = ((2 / k^2) - (omegaLoop (m, 1)).^2) ./ factorALoop;
                    factorIndALoop (m, 1) = ((4 / (LxLoop * LyLoop)) * sin ((omegaLoop (m, 2) * pi * pLoop (1) * LxLoop) / LxLoop) * sin ((omegaLoop (m, 3) * pi * pLoop (2) * LyLoop) / LyLoop)) ./ (factorALoopAll);
                    phiOutLLoop (m, 1) = (4 / (LxLoop * LyLoop)) * sin ((omegaLoop (m, 2) * pi * qLLoop (1) * LxLoop) / LxLoop) * sin ((omegaLoop (m, 3) * pi * qLLoop (2) * LyLoop) / LyLoop);
                    phiOutRLoop (m, 1) = (4 / (LxLoop * LyLoop)) * sin ((omegaLoop (m, 2) * pi * qRLoop (1) * LxLoop) / LxLoop) * sin ((omegaLoop (m, 3) * pi * qRLoop (2) * LyLoop) / LyLoop);
                end
                
                index = zeros (1, i);
                i = 1;
                
                for m = 1 : M
                    if omegaLoop (m, 1) < 44100 * 2 
                        index (1, i) = m;
                        i = i + 1;
                    end
                end
                
                plugin.omega = omegaLoop;
                
                qNextLoopInd =      qNextLoop(index);
                qNowLoopInd =       qNowLoop(index);
                qPrevLoopInd =      qPrevLoop(index);
                factorBdALoopInd =  factorBdALoop(index);
                factorCdALoopInd =  factorCdALoop(index);
                factorIndALoopInd = factorIndALoop(index);
                phiOutLLoopInd =    phiOutLLoop(index);
                phiOutRLoopInd =    phiOutRLoop(index);
            else
                qNextLoopInd =      qNextLoop;
                qNowLoopInd =       qNowLoop;
                qPrevLoopInd =      qPrevLoop;
                factorBdALoopInd =  factorBdALoop;
                factorCdALoopInd =  factorCdALoop;
                factorIndALoopInd = factorIndALoop;
                phiOutLLoopInd =    phiOutLLoop;
                phiOutRLoopInd =    phiOutRLoop;
                
                index = 1 : plugin.lengthOmega;
            end
            phiOutLPhase = plugin.phiOutLPre;
            phiOutRPhase = plugin.phiOutRPre;
            
            curSamp = plugin.currentSample;
            lengthCircX = plugin.circXLength;
            for t = 1 : length(in)
                if t == 1
                    qNextLoopInd = (factorBdALoopInd .* qNowLoopInd + factorCdALoopInd .* qPrevLoopInd + factorIndALoopInd .* plugin.samp);
                else
                    qNextLoopInd = (factorBdALoopInd .* qNowLoopInd + factorCdALoopInd .* qPrevLoopInd + factorIndALoopInd .* in (t - 1, 1));
                end
                
                if plugin.flanging == true
                    phiOutLLoopInd = phiOutLPhase (:, floor (mod ((curSamp + t) / 8, lengthCircX) + 1));
                    phiOutRLoopInd = phiOutRPhase (:, floor (mod ((curSamp + t) / 16, lengthCircX) + 1));
                end
%                 
                out(t, 1) = plugin.wetness / 100 * 25000 * sum (qNextLoopInd .* phiOutLLoopInd) ...
                    + (1 - plugin.wetness / 100) * in(t, 1);
                out(t, 2) = plugin.wetness / 100 * 25000 * sum (qNextLoopInd .* phiOutRLoopInd) ...
                    + (1 - plugin.wetness / 100) * in (t, 1);
                qPrevLoopInd = qNowLoopInd;
                qNowLoopInd = qNextLoopInd;
            end
%             if toc > 0.015
%                 plugin.prevLengthOmega = plugin.lengthOmega;
%                 plugin.lengthOmega = plugin.lengthOmega - 100;
%                 disp(plugin.lengthOmega)
%             else
%                 if toc < 0.01
%                     plugin.prevLengthOmega = plugin.lengthOmega;
%                     plugin.lengthOmega = plugin.lengthOmega + 100;
%                     disp(plugin.lengthOmega)
%                 end
%             end
            %disp(get(gca, 'CurrentPoint'));
            plugin.samp = in (end);
            plugin.qNext (index) =  qNextLoopInd;
            plugin.qPrev (index) =  qPrevLoopInd;
            plugin.qNow (index) =   qNowLoopInd;
            plugin.currentSample = plugin.currentSample + length (in (:, 1));
%             figure(2);
%             scatter(1:length(omegaLoop),sort(omegaLoop(:,1)))
%             drawnow;
        end
        function reset (plugin)
            % This section contains instructions to reset the plugin
            % between uses or if the environment sample rate changes.
        end
        
        
    end
end