function [factorBdA, factorCdA, factorIndA, omega, phiOutL, phiOutR, phiOutLPre, phiOutRPre, p, qL, qR, circXLength, rho, h] = initPlate(Lx, Ly, options)
%% Set Global Variables
fs = 44100;
ca = 343; % Speed of sound in air
pa = 1.225; % Air density

 %% Set VA effects
% flanging = false; %make microphones move or not
% stretching = 0; %stretch or not (0 = false, 1 = true)

%% Obtain Options
delModes =      options (1);
square =        options (2);
calcCent =      options (3);
flanging =      options (4);
stretching =    options (5);
C =             options(6);
useCm =         options(7);

decay = 4;

%% Set plate parameters
% Lx = 2; %Plate width
% Ly = 1; %Plate height
h = 0.0005; %plate thickness (m)
material = 'steel';
if strcmp(material, 'steel') == 1
    rho = 7850; %Material Density (kg/m^3)
    E = 2e11; %Young's modulus
    v = 0.3; % Poissons ratio
else
    if strcmp(material, 'glass') == 1
        rho = 2600;
        E = 7e10;
        v = 0.25;
    end
end
kSquared = (E*h^2)/(12*rho*(1-v^2));

%% Set input/output points
p = [0.4 0.415]; %input position at (0-1)
p2 = [0.9 0.5]; %input position at (0-1)
%p = [0.1 0.5]; %input position at (0-1)
qL = [0.1 0.45]; %left output position at (0-1)
qR = [0.84 0.45]; %right output position at (0-1)

%set to plate dimensions
in =    [p (1) * Lx     p(2) * Ly];
in2 =   [p2 (1) * Lx    p2(2) * Ly];
outL =  [qL (1) * Lx    qL(2) * Ly];
outR =  [qR (1) * Lx    qR(2) * Ly];

%% Get input
offset = 0; %set offset (in samples) for testing after stretch
len = 4.2; %length of sound in seconds

%[sound, soundfs] = audioread('rhodes2.mp3');
% input = zeros(soundfs*len+soundfs*5,1);
% input2 = zeros(soundfs*len+soundfs*5,1);
%input(1+offset:length(sound(1:soundfs*len,1))+offset) = 1;
%input2(1+offset:length(sound(1:soundfs*len,1))+offset) = sin(2*pi*440*16*(1:soundfs*len)/fs)*4;
% input(1+offset:length(sound(1:soundfs*len,1))+offset) = sound(1:soundfs*len,1);
%% Create Eigenfrequencies
disp('Create Omega')
val = 0;
m = 1;
m1 = 1;
m2 = 1;
omega = zeros(100000,3);
while val < fs*2
    val = ((m1/Lx)^2 + (m2/Ly)^2)*sqrt(kSquared)*pi^2;
    %val = sqrt(((pi^4*E*h^2)/(12*rho*(1-v^2)))*(m1^2/Lx^2+m2^2/Ly^2)^2); %Cummings paper
    if val < fs*2
        omega(m,1) = val;
        omega(m,2) = m1;
        omega(m,3) = m2;
        m2 = m2 + 1;
        m = m+1;
    else
        if m1 == 1
            stablem2 = m2-1;
        end
        m2 = 1;
        m1 = m1 + 1;
        val = ((m1/Lx)^2 + (m2/Ly)^2)*sqrt(kSquared)*pi^2;  
        if val > fs*2
            stablem1 = m1 - 1;
            break;
        end 
    end
end
omega = omega(1:m-1,:);

%% Make modes square
if square == true
    highMode = omega(:,2).*omega(:,3);
    indexFound = find(highMode==max(highMode(:)));
    i = 1;
    maxM1 = omega(indexFound(1),2);
    maxM2 = omega(indexFound(1),3);
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
    omega = sortrows(omega,1);
    n = 1;
    omegaPrev = 0;
    %C = 0.1;
    %omega(:,4) = ones(length(omega(:,1)),1);
    ncent = nthroot(2,12)^(C/100)*(omega(1,1)/(2*pi))-omega(1,1)/(2*pi);
    while n < length(omega(:,1))
        if omega(n,1)/(2*pi) - omegaPrev/(2*pi) < ncent
            %omega(n-1,4) = omega(n-1,4) + 1;
            omega(n,:) = [];
        else
            omegaPrev = omega(n,1);
            n = n+1;
            ncent = nthroot(2,12)^(C/100)*(omega(n,1)/(2*pi))-omega(n,1)/(2*pi);
        end 
    end
    %omega = sortrows(omega,[2,3]);
end
%omega = omega(1:1000,:);
M = length(omega(:,1));
phiIn = zeros(M,1);
for m = 1:M
    phiIn(m,1) = (4/(Lx*Ly))*sin((omega(m,2)*pi*in(1))/Lx)*sin((omega(m,3)*pi*in(2))/Ly);
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
Sx = 4;
Sy = 3;

%Create possible output positions
circX = outputPointsX(ceil(length(outputPointsX)*(Rx*sin(Sx*2*pi*(1:2/(max([Lx Ly])):fs)/fs)+0.5)));
circY = outputPointsY(ceil(length(outputPointsY)*(Ry*sin(Sy*2*pi*(1:2/(max([Lx Ly])):fs)/fs + 0.5*pi)+0.5)));

%% Create the Output Vector
phiOutL = zeros(M,1);
phiOutR = zeros(M,1);
disp('Create PhiOut')


if flanging == 1
    phiOutLPre = zeros(M,length(circX));
    phiOutRPre = zeros(M,length(circX));
    for t = 1:length(circX)
        for m = 1:M
            phiOutLPre(m,t) = (4/(Lx*Ly))*sin((omega(m,2)*pi*circX(t))/Lx)*sin((omega(m,3)*pi*circY(t))/Ly);
            phiOutRPre(m,t) = (4/(Lx*Ly))*sin((omega(m,2)*pi*circX(t))/Lx)*sin((omega(m,3)*pi*circY(t))/Ly);
        end
    end
else
    phiOutLPre = zeros(M,length(circX));
    phiOutRPre = zeros(M,length(circX));
end
for m = 1:M
    phiOutL(m,1) = (4/(Lx*Ly))*sin((omega(m,2)*pi*qL(1)*Lx)/Lx)*sin((omega(m,3)*pi*qL(2)*Ly)/Ly);
    phiOutR(m,1) = (4/(Lx*Ly))*sin((omega(m,2)*pi*qR(1)*Lx)/Lx)*sin((omega(m,3)*pi*qR(2)*Ly)/Ly);
end


creaMesh = 0;
if creaMesh == 1
    disp('Create Mesh')
   % La = 3*Lx;
   % Lb = 3*Ly;
   % mass = h*Lx*Ly*rho;
    gridSize = 100;
   % mode = 1;
    meshFunc = zeros(gridSize*Ly,gridSize*Lx);
    for x = 1:gridSize*Lx
        for y = 1:gridSize*Ly
            for mode = 1:length(omega(:,1))
                meshFunc(y,x,mode) = (4/(Lx*Ly))*sin(omega(mode,2)*pi*(x/(gridSize))/Lx)*sin(omega(mode,3)*pi*(y/(gridSize))/Ly);
                %meshFunc(y,x,mode) = (2/sqrt(mass))*sin(((x/gridSize)+Lx/2-La/2)*omega(mode,2)*pi*(x/gridSize)/Lx)*sin(((y/gridSize)+Ly/2-Lb/2)*omega(mode,3)*pi*(y/gridSize)/Ly);
            end
        end
    end
end
%Setoutput
%output = zeros(2,length(input));

%% Calculate thermoelastic damping
%values for the figure
% omega = [];
% freq(:,1) = 0:20000;
% omega(:,1) = freq*2*pi;

%Set thermal coefficients
R1 = 4.94e-3;
C1 = 2.98e-4;
%Calculate damping coefficient and damping factor
n1 = (omega(:,1)*R1*C1)./((omega(:,1).^2*h^2)+((C1^2)/(h^2)));
alphaTH = (omega(:,1)/2).*n1;

%% Calculate Radiation Damping
pa = 1.225;
fc = (ca^2)/(2*pi*sqrt(kSquared)); %calculate critical frecuency
phiRad = sqrt((omega(:,1)/(2*pi))/fc); 
g = ((1-phiRad.^2).*log((1+phiRad)./(1-phiRad))+2.*phiRad)./((1-phiRad.^2).^(3/2));
alphaRadPre = (1/(4*pi^2))*(ca*pa)/(rho*h)*((2*(Lx+Ly))/(Lx*Ly))*(ca/fc);
alphaRad = alphaRadPre.*g;

%% Calculate Damping induced by Porous Medium
zeroCheck = false;
alphaPorousTest = zeros(length(omega(:,1)),1);
for i = 1:length(omega(:,1))
    if zeroCheck == false
        alphaPorousTest(i,1) = sin(i/100)+1;
    else
        alphaPorousTest(i,1) = 0;
    end
    if alphaPorousTest(i,1) < 0.001
        zeroCheck = true;
    end
end
waveNo = omega(:,1)./ca;
a = ca./(omega(:,1)*(2*pi)*2);

%sigma = 1./sqrt(1-(pi./(waveNo.*a)));

%% Porous Medium Test
% omega = sortrows(omega,1);
% EPor = 11e7;
% rhoPor = 39;
% vPor = 0.46;
% cPor = sqrt((EPor*(1-vPor))/(rhoPor*(1+vPor)*(1-2*vPor)));
% 
% waveNo = omega(:,1)./ca;
% waveNoPor = omega(:,1)./cPor;
% a = ca./(omega(:,1)*(2*pi)*2);
% D = (E*h^3)/(12*(1-v^2));
% a = pi./(sqrt(omega(:,1)).*nthroot((rho/D),4));
% aPor = cPor./(omega(:,1)*(2*pi)*2);
% kY1 = sqrt(waveNo.^2-(pi./a).^2);
% kY2 = sqrt(waveNoPor.^2-(pi./a).^2);
% 
% p1 = exp(i.*omega(:,1)*t).*sin((pi*x)./(a.*Lx)*sin((pi*y)./(a.*Ly));
% p2;
% p3;

%% Calculate total damping
alphaTot = alphaRad+alphaTH;
T60 = 3.*log(10)./alphaTot; %reverberation time
cm = zeros(length(omega),1);
if useCm == true
    cm(1:length(omega(:,1)),1) = 12.*(log(10)./T60); %aloss coefficients
else
    cm(1:length(omega(:,1)),1) = 12.*(log(10)./decay); %loss coefficients
end
%% Initialise update equation
k = 1/fs;
% qNext = zeros(length(omega(:,1)),1);
% qPre = zeros(length(omega(:,1)),1);
% qNow = zeros(length(omega(:,1)),1);
% qPrev = zeros(length(omega(:,1)),1);
factorA = (1/k^2)+(cm/(rho*h*k));
factorB = ((2/k^2)-(omega(:,1)).^2);
factorIn = ((phiIn)./(rho*h));
factorIndA = factorIn./factorA;
factorC = ((cm/(rho*h*k))-(1/(k^2)));
factorCdA = factorC./factorA;
factorBdA = factorB./factorA;

circXLength = length(circX);
