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
elecpop = 10000;                     % Number of particles to simulate
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
vx = Vth*randn(elecpop,1)/sqrt(2);
vy = Vth*randn(elecpop,1)/sqrt(2); 
v = sqrt( vx.^2 + vy.^2 );
    figure(1)
    histogram(v,100)
    title('Histogram of the Velocity Distribution')
    xlabel('Baskets')
    ylabel('Velocities (m/s)')
    
    PrevScatter = zeros(elecpop,1);
    TScatter = zeros(elecpop,1);
    AvgScatter = zeros(iter,1);
    MFPCalc = zeros(elecpop,1);
    MFP = zeros(iter,1);
    for j = 1:iter 
        oldx = x;
        oldy = y;
        PartScatter = Pscatter > rand(elecpop,1);
        TScatter(PartScatter) = (j - PrevScatter(PartScatter))*dt;
        PrevScatter(PartScatter) = j;
        AvgScatter(j) = mean(TScatter(PartScatter));
        MFPCalc(PartScatter) = sqrt( vx(PartScatter).^2 + vy(PartScatter).^2 ).*TScatter(PartScatter);
        MFP(j) = mean(MFPCalc(PartScatter));
        vx(PartScatter) = Vth*randn(sum(PartScatter),1)/sqrt(2);
        vy(PartScatter) = Vth*randn(sum(PartScatter),1)/sqrt(2); 
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
            title('Sample Particles Trajectories')
            xlabel('W (m) or x position')
            ylabel('L (m) or y position')
            axis([0 200e-9 0 100e-9])
            hold on
        end
            V = mean( sqrt( vx.^2+vy.^2 ) );
            NewT = V.^2*mn/(2*k);  
            figure(3)% Plotting T
            plot(j,NewT,'o')
            hold on
            title('Temperature at random points')
            xlabel('random occurence')
            ylabel('Temperature (K)')
    end 
    mfpCalc = mean(MFP);
    tau = mean(AvgScatter);