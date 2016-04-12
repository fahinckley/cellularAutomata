%--------------------------------
% Evaluates the flow around an airfoil
%   Modified from Iain Haslam's code
%--------------------------------
% Franklin Hinckley
% 11 April 2016
%--------------------------------
% Comments from Ian Haslam's code for the possible states
%  c4  c3   c2  D2Q9 model. At each timestep, particle densities propagate
%    \  |  /    outwards in the directions indicated in the figure. An
%  c5 -c9 - c1  equivalent 'equilibrium' density is found, and the densities
%    /  |  \    relax towards that state, in a proportion governed by omega.
%  c6  c7   c8      Iain Haslam, March 2006.
%--------------------------------

%% Clean up
clearvars
close all
clc

%% Define time scale and initial density
omega = 1.0;    % = 1/tau where tau=relaxation time scale
density = 1.0;  % initial density

%% Define density weights
% These are the weights for equilibrium distribution
t1 = 4/9;   % Rest particle
t2 = 1/9;   % "Straight" directions
t3 = 1/36;  % Diagonal directions

% Set up c-squared scale factor
c_squ = 1/3;  

%% Boundaries
% Generate domain
x = linspace(0,2,60);
y = linspace(-0.5,0.5,60);
nx = length(x);
ny = length(y);

% Get segment end-points for NACA airfoil
numP = 100;
[tmpx,tmpy] = NACA(4,4,08,0,1,numP);
tmpx = tmpx + 0.5;

% Mark wing shape in boundaries matrix
BOUND = zeros(length(y),length(x));
% Top half
for ii = 1:length(x)  
    % Find index in wing-x array
    xw = find(tmpx(1:numP) < x(ii),1,'last');
    if ~isempty(xw)
        for jj = 1:length(y)   
            % Check if y coordinate is less than wing surface
            if (y(jj) > 0) && (y(jj) < tmpy(xw))
                BOUND(jj,ii) = 1;
            end        
        end    
    end
end
% Bottom half
for ii = 1:length(x)  
    % Find index in wing-x array
    xw = find(tmpx(numP+1:end) < x(ii),1,'first');
    if ~isempty(xw)
        for jj = 1:length(y)   
            % Check if y coordinate is less than wing surface
            if (y(jj) < 0) && (y(jj) > tmpy(xw+numP))
                BOUND(jj,ii) = 1;
            end        
        end    
    end
end

% Find matrix index of each occupied node
BOUND = fliplr(BOUND');
ON = find(BOUND); 

% There are nine possible states at each node, each "layer" in F tracks
%   one of the possible states
F = repmat(density/9,[nx ny 9]); 

% Mark equilibrium state
FEQ = F; 

% Get layer size
msize = nx*ny; 
CI = [0:msize:msize*7];

% Identify points at which the flow will reflect
TO_REFLECT = [ON+CI(1) ON+CI(2) ON+CI(3) ON+CI(4) ...
              ON+CI(5) ON+CI(6) ON+CI(7) ON+CI(8)];
REFLECTED =  [ON+CI(5) ON+CI(6) ON+CI(7) ON+CI(8) ...
              ON+CI(1) ON+CI(2) ON+CI(3) ON+CI(4)];
        
avu = 1;         % Average x-directed velocity (used to assess equilibrium)
prevavu = 1;     % Previous average velocity (used to assess equilibrium)
ts = 0;          % Time-step counter
deltaU = 1e-7;   % Velocity added to first row (left side) to drive simulation
numactivenodes = sum(sum(1-BOUND)); % double sum covers both dimensions

%% Main loop
% Runs at least 100 time steps,
%   terminates after either fixed iterations or at equilibrium
while ((ts < 5000) && (1e-10 < abs((prevavu-avu)/avu))) || (ts < 100)
    
    % Propagate positions based on state at each node
    F(:,:,4) = F([2:nx 1],[ny 1:ny-1],4);     % up and right
    F(:,:,3) = F(:,[ny 1:ny-1],3);            % right
    F(:,:,2) = F([nx 1:nx-1],[ny 1:ny-1],2);  % down and right
    F(:,:,5) = F([2:nx 1],:,5);               % up
    F(:,:,1) = F([nx 1:nx-1],:,1);            % propagate down
    F(:,:,6) = F([2:nx 1],[2:ny 1],6);        % up and left
    F(:,:,7) = F(:,[2:ny 1],7);               % left
    F(:,:,8) = F([nx 1:nx-1],[2:ny 1],8);     % down and left
    
    % Identify nodes to reflect (boundary points)
    BOUNCEDBACK = F(TO_REFLECT);
    
    % Sum across layers (3rd dimension of array) to get density at each node
    DENSITY = sum(F,3);
    
    % Calculate X velocity ("left" - "right")
    UX = (sum(F(:,:,[1 2 8]),3) - sum(F(:,:,[4 5 6]),3))./DENSITY;
    
    % Calculate Y velocity ("up" - "down")
    UY = (sum(F(:,:,[2 3 4]),3) - sum(F(:,:,[6 7 8]),3))./DENSITY;
    
    % Apply boundary conditions:
    % Add velocity to left side (row 1).
    % Set velocity and density at boundary points to zero.
    % Pre-compute squares, etc.
    UX(1,1:ny) = UX(1,1:ny) + deltaU; %Increase inlet pressure
    UX(ON) = 0; 
    UY(ON) = 0; 
    DENSITY(ON) = 0;
    U_SQU = UX.^2 + UY.^2; 
    U_C2  = UX + UY; 
    U_C4  = -UX + UY; 
    U_C6  = -U_C2; 
    U_C8  = -U_C4;
    
    % Calculate equilibrium distribution: stationary
    FEQ(:,:,9) = t1*DENSITY.*(1-U_SQU/(2*c_squ));
    
    % nearest-neighbours
    FEQ(:,:,1) = t2*DENSITY.*(1+UX/c_squ+0.5*(UX/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,3) = t2*DENSITY.*(1+UY/c_squ+0.5*(UY/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,5) = t2*DENSITY.*(1-UX/c_squ+0.5*(UX/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,7) = t2*DENSITY.*(1-UY/c_squ+0.5*(UY/c_squ).^2-U_SQU/(2*c_squ));
    % next-nearest neighbours
    FEQ(:,:,2) = t3*DENSITY.*(1+U_C2/c_squ+0.5*(U_C2/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,4) = t3*DENSITY.*(1+U_C4/c_squ+0.5*(U_C4/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,6) = t3*DENSITY.*(1+U_C6/c_squ+0.5*(U_C6/c_squ).^2-U_SQU/(2*c_squ));
    FEQ(:,:,8) = t3*DENSITY.*(1+U_C8/c_squ+0.5*(U_C8/c_squ).^2-U_SQU/(2*c_squ));
    
    % Compute new particle velocities
    F = omega*FEQ + (1-omega)*F;
    
    % Bounce-back from boundaries
    F(REFLECTED) = BOUNCEDBACK;
    
    % Update time step and check for equilibrium
    prevavu = avu;
    avu = sum(sum(UX))/numactivenodes; 
    ts = ts+1;
end

%% Plot final state
figure
colormap(gray(2))
image(2-BOUND')
hold on
quiver(2:nx,1:ny,UX(2:nx,:)',UY(2:nx,:)')
h = streamline(2:nx,1:ny,UX(2:nx,:)',UY(2:nx,:)',2*ones(length(1:3:ny),1),1:3:ny);
for ii = 1:length(h)
    plot(h(ii).XData,h(ii).YData,'r');
end
title(['Flow field after ',num2str(ts),'\deltat'])
xlabel('x')
ylabel('y')
