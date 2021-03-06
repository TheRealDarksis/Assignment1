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
Vth = sqrt(2*k*T/mn);               % Thermal velocity in m/s
mfp = Vth*tau_mn;                   % Mean free path in m
samplepop = 10;                     % # of particles to plot
samp = randi(elecpop,samplepop,1);  % Random particles to observe
iter = 100;                         % # of iterations 
Pscatter = 1 - exp(-dt/tau_mn);     % Scattering probability

% Initializing positions and velocities 
    y(:,1) = rand(elecpop,1)*L;
    vx = Vth*randn(elecpop,1)/sqrt(2);
    vy = Vth*randn(elecpop,1)/sqrt(2); 
    for z = 1:elecpop   % Generating points on the left side
       x(z,1) = 0.2e-7 + (0.8e-7 - 0.2e-7)*rand();
    end
    for z = elecpop/2+1:elecpop     % Generating points on the right
       x(z,1) = 1.2e-7 + (2e-7 - 1.2e-7)*rand();        
    end
    
%Boxes
Boxes = {};
Boxes{1}.X = [0.8e-7, 0.4e-7];
Boxes{1}.Y = [0e-7, 0.4e-7];
Boxes{2}.X = [0.8e-7, 0.4e-7];
Boxes{2}.Y = [0.6e-7, 1e-7];

% circle for enhancement  
angle = 0:0.05:2*pi;
circx = cos(angle)';
circy = sin(angle)';
circle = [circx,circy];

    for j = 1:iter 
        %Creating box reflection
        vx(oldx>=1.2e-7 & x<=1.2e-7 & (y>0.6e-7 | y<0.4e-7)) = -vx(oldx>=1.2e-7 & x<=1.2e-7 & (y>0.6e-7 | y<0.4e-7));
        vx(oldx<=0.8e-7 & x>=0.8e-7 & (y>0.6e-7 | y<0.4e-7)) = -vx(oldx<=0.8e-7 & x>=0.8e-7 & (y>0.6e-7 | y<0.4e-7));        
        vy(x>0.8e-7 & x<1.2e-7 & oldy>=0.4e-7 & y<=0.4e-7) = -vy(x>0.8e-7 & x<1.2e-7 & oldy>=0.4e-7 & y<=0.4e-7);
        vy(x>0.8e-7 & x<1.2e-7 & oldy<=0.6e-7 & y>=0.6e-7) = -vy(x>0.8e-7 & x<1.2e-7 & oldy<=0.6e-7 & y>=0.6e-7);
        %Code for trajectory and collision
        oldx = x;
        oldy = y;
        PartScatter = Pscatter > rand(elecpop,1);
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
%         figure(2)        
%         circle2(0,0.5e-7,0.2e-7);   % Plotting circle
%         for m = 1:samplepop
%             for a = 1:2
%                 rectangle('Position',[Boxes{a}.X(1), Boxes{a}.Y(1), Boxes{a}.X(2), Boxes{a}.Y(2)])
%                 axis([0 W 0 L])
%             end            
%             plot([oldx(samp(m)), x(samp(m))],[oldy(samp(m)),y(samp(m))],'SeriesIndex',m)
%             title('Sample Particles Trajectories')
%             xlabel('W (m) or x position')
%             ylabel('L (m) or y position')
%             axis([0 200e-9 0 100e-9])
%             hold on 
%         end
            V = sqrt( vx.^2+vy.^2 ); 
        pause(0.01)
    end

% % Plotting electron density maps
% figure
% hist3([x(:), y(:)],'CDataMode','auto')
% xlabel('W (m) or x position')
% ylabel('L (m) or y position')
% title('Electron density map')

%temp map!!!!
nx = 20; ny = 30;
BinX = ceil( x*nx/W );
BinY = ceil( y*ny/L );
for i = 1:nx
    for j = 1:ny
        match = BinX == i & BinY == j;
        Tbox(i,j) = sum( V(match).^2/(2*k) );
    end
end
figure
surf(Tbox)
title('Temperature Density')
xlabel('Width (m)')
ylabel('Length (m)')
zlabel('Temperature (K)')