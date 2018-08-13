function [coeffBdA, coeffCdA, coeffIndA, kSquared, omega, phiOutL, phiOutR, phiOutFlange, circXLength, circX, circY] = initPlate (Lx, Ly, C, rho, T60, h, flangeMatSize, inOutputs)

fs = 44100;

%% Obtain in- and outputs
p = inOutputs(1, :);
qL = inOutputs(2, :);
qR = inOutputs(3, :);

%% Set plate parameters (using values from the EMT 140 Plate Reverb)
E = 2e11;   % Young's modulus
v = 0.3;    % Poissons ratio
kSquared = (E * h^2) / (12 * rho * (1 - v^2));

%% Create Eigenfrequencies
disp('Create Omega')
val = 0;
m = 1;
m1 = 1;
m2 = 1;
omega = zeros(18218, 3); % 18218 is the maximum number of modes

% While the value of the eigenfrequency is smaller than the stability
% condition of twice the sample rate, fill the matrix
while val < fs * 2 
    val = ((m1 / Lx)^2 + (m2 / Ly)^2) * sqrt(kSquared) * pi^2;
    
    if val < fs * 2
        omega (m, 1) = val;
        omega (m, 2) = m1;
        omega (m, 3) = m2;
        m2 = m2 + 1;
        m = m + 1;
    else
        m2 = 1;
        m1 = m1 + 1;
        val = ((m1 / Lx)^2 + (m2 / Ly)^2) * sqrt(kSquared) * pi^2;  
    end
end
omega = omega (1 : m - 1,:);

%% Delete Neglegible Modes From Input
% Modes that have a node at the input position can be neglected.

disp('Delete neglegible modes from input')

% Find which modes can be neglectd
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

% Delete the modes from the vector
n = 1;
while n < length(omega(:,1))
    if mod(omega(n,2), i1) == 0 || mod(omega(n,3), i2) == 0
        omega(n,:) = [];
    else
        n = n + 1;
    end
end

%% Delete Perceptually Unimportant Modes
disp('Delete Perceptually Unimportant Modes')

% Sort the omega matrix according to the eigenfrequencies
omega = sortrows (omega, 1);
n = 1;
omegaPrev = 0;

% Set the first threshold
ncent = nthroot(2, 12)^(C / 100) * (omega(1, 1) / (2 * pi)) - omega(1, 1) / (2 * pi);
while n < length(omega(:, 1))
    % If the difference in cents between two eigenfrequencies is smaller
    % than the treshold, delete the mode.
    if omega(n, 1) / (2 * pi) - omegaPrev / (2 * pi) < ncent
        omega(n,:) = [];
    else
    % Otherwise set a new threshold
        omegaPrev = omega(n, 1);
        n = n + 1;
        ncent = nthroot(2, 12)^(C / 100) * (omega(n, 1) / (2 * pi)) - omega(n, 1) / (2 * pi);

    end 
end

%% Create the Input Vector
phiIn = (2 / (Lx * Ly)) * sin(omega(:, 2) * pi * p(1)) .* sin(omega(:, 3) * pi * p(2));

%% Set up Moving Outputs
disp('Set up Moving Outputs')

% The parameters below (extremes and speeds) have been preset and can't be changed by the user.
% These have been selected as having the best quality for the purpose of creating a flanging effect

% Set shape extremes [-0.5 - 0.5]
Rx = 0.4;
Ry = 0.4;

% Set x and y speeds (Sx should be even, Sy odd)
Sx = 2;
Sy = 3;

%Create possible output positions
circX = Rx * sin (Sx * 2 * pi * (1 : 2 / (max([Lx Ly])) : (flangeMatSize)) / (flangeMatSize)) + 0.5;
circY = Ry * sin (Sy * 2 * pi * (1 : 2 / (max([Lx Ly])) : (flangeMatSize)) / (flangeMatSize) + 0.5 * pi) + 0.5;
circXLength = length(circX);

phiOutFlange = sin (omega (:, 2) * pi * circX) .* sin (omega (:, 3) * pi * circY);    

%% Create the Output Vector
phiOutL = sin (omega (:, 2) * pi * qL(1)) .* sin (omega (:, 3) * pi * qL(2));
phiOutR = sin (omega (:, 2) * pi * qR(1)) .* sin (omega (:, 3) * pi * qR(2));

%% Decay
decayVector = zeros (length (omega (:,1)), 1);
decayVector (1 : length (omega (:,1)), 1) = T60;
cm = 12 .* (log(10) ./ decayVector); 

%% Create coefficients In/A C/A and B/A used in the update equation
k = 1 / fs;
coeffA = (1 / k^2) + (cm / (rho * h * k));
coeffB = ((2 / k^2) - (omega (:, 1)).^2);
coeffC = ((cm / (rho * h * k)) - (1 / (k^2)));
coeffIn = ((phiIn) ./ (rho * h));
coeffIndA = coeffIn ./ coeffA;   % In divided by A
coeffBdA = coeffB ./ coeffA;     % B divided by A
coeffCdA = coeffC ./ coeffA;     % C divided by A
