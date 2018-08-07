function [coeffBdA, coeffCdA, coeffIndA, omega, phiOutL, phiOutR, phiOutFlange, circXLength, rho, cm] = initPlate (Lx, Ly, C, inOutputs, settings)
%% Set Global Variables
fs = 44100;
ca = 343; % Speed of sound in air
pa = 1.225; % Air density

%% Obtain Settings
delModes =      settings (1);
square =        settings (2);
calcCent =      settings (3);
flanging =      settings (4);
useCm =         settings (5);

decay = 4;

%% Obtain in- and outputs
p = inOutputs(1, :);
qL = inOutputs(2, :);
qR = inOutputs(3, :);

%% Set plate parameters (using values from the EMT 140)
h = 0.0005; % Plate thickness (m)
rho = 7850; % Material Density (kg/m^3)
E = 2e11;   % Young's modulus
v = 0.3;    % Poissons ratio
kSquared = (E * h^2) / (12 * rho * (1 - v^2));

%set to plate dimensions
in =    [p(1) * Lx      p(2) * Ly];
outL =  [qL (1) * Lx    qL(2) * Ly];
outR =  [qR (1) * Lx    qR(2) * Ly];

%% Get input
offset = 0; %set offset (in samples) for testing after stretch
len = 4.2; %length of sound in seconds

%% Create Eigenfrequencies
disp('Create Omega')
val = 0;
m = 1;
m1 = 1;
m2 = 1;
omega = zeros(100000,3);

while val < fs * 2
    val = ((m1 / Lx)^2 + (m2 / Ly)^2) * sqrt (kSquared) * pi^2;
    
    if val < fs*2
        omega (m, 1) = val;
        omega (m, 2) = m1;
        omega (m, 3) = m2;
        m2 = m2 + 1;
        m = m + 1;
    else
        if m1 == 1
            stablem2 = m2-1;
        end
        m2 = 1;
        m1 = m1 + 1;
        val = ((m1 / Lx)^2 + (m2 / Ly)^2) * sqrt (kSquared) * pi^2;  
        if val > fs * 2
            stablem1 = m1 - 1;
            break;
        end 
    end
end
omega = omega (1 : m - 1,:);

%% Make modes square
if square == true
    highMode = omega (:, 2) .* omega(:, 3);
    indexFound = find (highMode == max (highMode (:)));
    i = 1;
    maxM1 = omega(indexFound(1), 2);
    maxM2 = omega(indexFound(1), 3);
    while i <= length(omega(:,1))
        if omega(i,2) > maxM1 || omega(i,3) > maxM2
            omega(i,:) = [];
            %indexFound = indexFound - 1;
        else
            i = i + 1;
        end
    end
end

%% Delete Neglegible Modes From Input
if delModes == true
    disp('Delete neglegible modes from input')
    i1 = 1;
    answ1 = p(1);
    while rem(answ1,1) ~= 0 
        i1 = i1+1;
        answ1 = p(1)*i1;
    end

    i2 = 1;
    answ2 = p(2);
    while rem(answ2,1) ~= 0 
        i2 = i2+1;
        answ2 = p(2)*i2;
    end

    n = 1;
    while n < length(omega(:,1))
        if mod(omega(n,2), i1) == 0 || mod(omega(n,3), i2) == 0
            omega(n,:) = [];
        else
            n = n + 1;
        end
    end
end
%% Calculate Cents
%omega = sortrows(omega,1);
if calcCent == 1
    disp('Calculate Cents')
    omega = sortrows (omega, 1);
    n = 1;
    omegaPrev = 0;
    %C = 0.1;
    %omega(:,4) = ones(length(omega(:,1)),1);
    index = [1:length(omega(:, 1))]';
    i = 1;
    ncent = nthroot(2, 12)^(C / 100) * (omega(1, 1) / (2 * pi)) - omega(1, 1) / (2 * pi);
    while n < length(omega(:,1))
        if omega(n, 1) / (2 * pi) - omegaPrev / (2 * pi) < ncent
            omega(n,:) = [];
            index(n) = [];
        else
            omegaPrev = omega(n, 1);
            n = n + 1;
            ncent = nthroot(2, 12)^(C / 100) * (omega(n, 1) / (2 * pi)) - omega(n, 1) / (2 * pi);
            
        end 
        i = i + 1;
    end
    %omega = sortrows(omega,[2,3]);
end

%omega = omega(1:1000,:);
M = length(omega(:, 1));
phiIn = zeros(M, 1);
for m = 1:M
    phiIn(m, 1) = (4 / (Lx * Ly)) * sin((omega(m, 2) * pi * in(1)) / Lx) * sin((omega(m, 3) * pi * in(2)) / Ly);
    %phiIn2(m,1) = (4/(Lx*Ly))*sin((omega(m,2)*pi*in2(1))/Lx)*sin((omega(m,3)*pi*in2(2))/Ly);
end

%% Set up Moving Outputs

disp('Set up Moving Outputs')
outputPointsX = 0:1/(10000*Lx):1;
outputPointsY = 0:1/(10000*Ly):1;
outputPointsX = outputPointsX*Lx;
outputPointsY = outputPointsY*Ly;

%set shape extremes
Rx = 0.4;
Ry = 0.4;

%set x and y speeds
Sx = 2;
Sy = 3;

%Create possible output positions
circX = outputPointsX (ceil (length (outputPointsX) * (Rx * sin (Sx * 2 * pi * (1 : 2 / (max([Lx Ly])) : (fs / 4)) / (fs / 4)) + 0.5)));
circY = outputPointsY (ceil (length (outputPointsY) * (Ry * sin (Sy * 2 * pi * (1 : 2 / (max([Lx Ly])) : (fs / 4)) / (fs / 4) + 0.5 * pi) + 0.5)));
circXLength = length(circX);

%% Create the Output Vector

disp('Create PhiOut')


phiOutFlange = (4 / (Lx * Ly)) * sin ((omega (:, 2) * pi * circX) / Lx) .* sin ((omega (:, 3) * pi * circY) / Ly);    

phiOutL = (4 / (Lx * Ly)) * sin (omega (:, 2) * pi * qL(1)) .* sin (omega (:, 3) * pi * qL(2));
phiOutR = (4 / (Lx * Ly)) * sin (omega (:, 2) * pi * qR(1)) .* sin (omega (:, 3) * pi * qR(2));


%% Calculate thermoelastic damping
if useCm == true
    
    % Set thermal coefficients
    R1 = 4.94e-3;
    C1 = 2.98e-4;

    % Calculate damping coefficient
    n1 = (omega(:,1)*R1*C1) ./ ((omega(:,1).^2*h^2)+((C1^2)/(h^2)));
    alphaTH = (omega(:,1) / 2).*n1;

    %% Calculate Radiation Damping
    fc = (ca^2) / (2 * pi * sqrt (kSquared)); %calculate critical frecuency
    phiRad = sqrt ((omega (:, 1) / (2 * pi)) / fc); 
    g = ((1 - phiRad.^2) .* log ((1 + phiRad) ./ (1 - phiRad)) + 2 .* phiRad) ./ ((1 - phiRad.^2).^(3 / 2));
    alphaRadPre = (1 / (4 * pi^2)) * (ca * pa) / (rho * h) * ((2 * (Lx + Ly)) / (Lx * Ly)) * (ca / fc);
    alphaRad = alphaRadPre .* g;

    %% Calculate total damping
    alphaTot = alphaRad+alphaTH;
    T60 = 3 .* log(10) ./ alphaTot; % Reverberation time per mode
    cm = 12 .* (log(10) ./ T60); % Physical Damping
else
    decayVector = zeros (length (omega (:,1)), 1);
    decayVector (1 : length (omega (:,1)), 1) = decay;
    cm = 12 .* (log(10) ./ decayVector); % Damping using a reverberation time of <decay> for all modes 
end

%% Create coefficients In/A C/A and B/A used in the update equation
k = 1 / fs;
coeffA = (1 / k^2) + (cm / (rho * h * k));
coeffB = ((2 / k^2) - (omega (:, 1)).^2);
coeffC = ((cm / (rho * h * k)) - (1 / (k^2)));
coeffIn = ((phiIn) ./ (rho * h));
coeffIndA = coeffIn ./ coeffA;   % In divided by A
coeffBdA = coeffB ./ coeffA;     % B divided by A
coeffCdA = coeffC ./ coeffA;     % C divided by A
