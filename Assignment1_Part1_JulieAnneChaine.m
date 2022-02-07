%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        ELEC 4700 A - Winter 2022                        %
%                              Assignment 1                               %
%        Author: Julie-Anne Chaine             Student #: 101104568       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% Constants:
mo = 9.1093837015e-31;  % Electron rest mass in kg
mn = 0.26*mo;           % Effective mass of electrons
k = 1.38e-23;           % Boltzmann's constant 
W = 200e-9;             % Nominal region width in m
L = 100e-9;             % Nominal region length in m
T = 300;                % Room temperature in K
tau_mn = 0.2e-12;       % Room temperature in K
elecpop = 1000;         % Number of particles to simulate
iter = 100;             % # of iterations 

% Variables
x = zeros(elecpop,1);               % Electron y positions
y = zeros(elecpop,1);               % Electron x positions
oldx = zeros(elecpop,1);            % Previous electron x position
oldy = zeros(elecpop,1);            % Previous electron y position
phi = zeros(elecpop,1);             % Electron direction
vx = zeros(elecpop,1);              % Electron x direction velocity
vy = zeros(elecpop,1);              % Electron y direction velocity
dt = 1e-14;                         % Time step in s
samplepop = 7;                      % # of electrons to observe
samp = randi(elecpop,samplepop,1);  % Random selection of particles
SampVx = zeros(7,1);                % Vector to calculate T later
SampVy = zeros(7,1);                % Vector to calculate T later
SampV = zeros(7,1);                 % Vector to calculate T later

% For Q1, thermal velocity
Vth = sqrt(2*k*T/mn);   % Thermal velocity in m/s

% For Q2, mean free path
mfp = Vth*tau_mn;   % Mean free path in m

% Code for Q3
% Random position generation in m
x(:,1) = rand(elecpop,1)*W;
y(:,1) = rand(elecpop,1)*L;
    
% Random direction generation in rad
phi(:,1) = rand(elecpop,1)*2*pi;
vx(:,1) = Vth*cos(phi(:,1)); 
vy(:,1) = Vth*sin(phi(:,1));
    
    for j = 1:iter 
        oldx = x;
        oldy = y;
        x = oldx + vx*dt;       % Previous x position + delta L
        y = oldy + vy*dt;       % Previous y position + delta L
        oldx(x<0) = W;          % Making the particles on the left boundary, appear on the right
        oldx(x>W) = 0;          % Making the particles on the right boundary, appear on the left
        x(x<0) = x(x<0) + W;    % All points passed the left get pushed to the other side
        x(x>W) = x(x>W) - W;    % All points passed the right get pushed to the other side
        vy(y>L) = -vy(y>L);     % Particle Y direction gets flipped if it hits the top
        vy(y<0) = -vy(y<0);     % Particle Y direction gets flipped if it hits the bottom
        for m = 1:samplepop     % For plotting
            figure(1)
            plot([oldx(samp(m)), x(samp(m))],[oldy(samp(m)),y(samp(m))],'SeriesIndex',m)
            axis([0 200e-9 0 100e-9])
            hold on
            SampVx(m) = vx(samp(m));    % To calculate T
            SampVy(m) = vy(samp(m));
            SampV(m) = sqrt( SampVx(m)^2+SampVy(m)^2 );
            NewT(m) = SampV(m)^2*mn/(2*k*T);            
        end
        pause(0.01)
    end
    figure(2)   % Plotting T
    plot(1:7,NewT)
    title('Temperature at random points')
    xlabel('random occurence')
    ylabel('Temperature')
    


       
