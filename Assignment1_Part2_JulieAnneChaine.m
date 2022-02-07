%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        ELEC 4700 A - Winter 2022                        %
%                              Assignment 1                               %
%        Author: Julie-Anne Chaine             Student #: 101104568       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% Constants:
mo = 9.1093837015e-31;              % Electron rest mass in kg
mn = 0.26*mo;                       % Effective mass of electrons
k = 1.38e-23;                       % Boltzmann's constant 
W = 200e-9;                         % Nominal region width in m
L = 100e-9;                         % Nominal region length in m
tau_mn = 0.2e-12;                   % Mean time between collisions
r = 2.82e-15;                       % Electron radius

%Variables
T = 300;                            % Room temperature in K
elecpop = 1000;                     % Number of particles to simulate
dt = 1e-14;                         % Time step in s
x = zeros(elecpop,1);               % Electron y positions
y = zeros(elecpop,1);               % Electron x positions
oldx = zeros(elecpop,1);            % Previous electron x position
oldy = zeros(elecpop,1);            % Previous electron y position
nx = zeros(elecpop,1);              % x velocity distribution
ny = zeros(elecpop,1);              % y velocity distribution
Vth = sqrt(2*k*T/mn);               % Thermal velocity in m/s
n = zeros(elecpop,1);               % Velocity distribution 
mfp = Vth*tau_mn;                   % Mean free path in m
samplepop = 10;                     % # of particles to plot
samp = randi(elecpop,samplepop,1);  % Random particles to observe
iter = 100;                         % # of iterations 
Pscatter = 1 - exp(-dt/tau_mn);     % Scattering probability

% Initializing positions and velocities 
    x(:,1) = rand(elecpop,1)*W;
    y(:,1) = rand(elecpop,1)*L;
    vx = Vth*randn(elecpop,1);
    vy = Vth*randn(elecpop,1);  
    for a = 1:elecpop   % Velocity distribution calculations
        nx(a,1) = sqrt(mn/(2*pi*k*T))*exp(-mn*vx(a,1)^2/(2*k*T));
        ny(a,1) = sqrt(mn/(2*pi*k*T))*exp(-mn*vy(a,1)^2/(2*k*T));
        n(a,1) = sqrt( nx(a,1)^2 + ny(a,1)^2 );
    end
    figure(1)
    histogram(n)
    
    for j = 1:iter 
        oldx = x;
        oldy = y;
        if Pscatter > rand()
            vx = Vth*randn(elecpop,1);
            vy = Vth*randn(elecpop,1); 
        end 
            x = oldx + vx*dt;       % Previous x position + delta L
            y = oldy + vy*dt;       % Previous y position + delta L
            oldx(x<0) = W;          % Making the particles on the left boundary, appear on the right
            oldx(x>W) = 0;          % Making the particles on the right boundary, appear on the left
            x(x<0) = x(x<0) + W;    % All points passed the left get pushed to the other side
            x(x>W) = x(x>W) - W;    % All points passed the right get pushed to the other side
            vy(y>L) = -vy(y>L);     % Particle Y direction gets flipped if it hits the top
            vy(y<0) = -vy(y<0);     % Particle Y direction gets flipped
        for m = 1:samplepop         % For plotting
            figure(2)
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
    figure(3)   % Plotting T
    plot(1:samplepop,NewT)
    title('Temperature at random points')
    xlabel('random occurence')
    ylabel('Temperature')