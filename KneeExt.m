function [ t,F,X,V,f,v,G,A ] = KneeExt
%% KNEEEXT A simple function to simulate a knee extension on a leg-press
% 1 Make sure the function is in the Matlab search path or current folder.
% 2 Call the function using 'KneeExt;' from the Command Window.
% 3 Use '[ t,F,X,V,f,v,G,A ] = KneeExt;' to access the output time series
%   in the main Matlab workspace, e.g., run: 
%    [ t,F,X,V,f,v,G,A ] = KneeExt;
%    plot(t,F,'b-',t,f,'r-');   grid on;
%    xlabel('time [s]');        ylabel('force [N]');
%    legend({'External Force','Muscle Force'})
% Inputs:
%   Input parameters are described in section 'INPUT PARAMETERS' below.
%   The user might use the provided but commented 'input' command to 
%   manually change/type (some) parameter values via the Command Window.
% Outputs:
%   t  ...  Time [s]
%   F  ...  External force [N]
%   X  ...  Position: distance from the prox. end of the model-thigh to
%           the distal end of the model-shank [m]
%   V  ...  Velocity of the accelerated point-mass [m/s]
%   f  ...  Force of the model mucscle [N]
%   v  ...  Contraction-velocity of the model mucscle [m/s]
%   G  ...  Values of the function of geometry [-]
%   A  ...  Values of the function of activation dynamics [-]
% $Date: Februay 6, 2018 by H. Penasso and S. Thaller                v.1.02
% ________________________________________
    disp('Knee extension');disp('KneeExt.m; H.Penasso, S.Thaller')
  %% INPUT PARAMETERS
  % Leg
    lt = 0.52; % input ('Length of thigh  (0.52) [m]: ');
    ls = 0.45; % input ('Length of shank  (0.45) [m]: ');
    kt = 0.52; % input ('Length of muscle (0.52) [m]: ');
    ks = 0.06; % input ('Patella center to tibial tub. (0.06) [m]: ');
    r  = 0.06; % input ('Joint radius     (0.06) [m]:');
  % Muscle
    a  = 5000; % input ('Hill's parameter a (5000)    [N]  : ');
    b  = 0.1;  % input ('Hill's parameter b (0.1)     [m/s]: ');
    c  = 4186; % input ('Hill's parameter c (4186)    [W]  : ');
    S  = 7;    % input ('Activation rate constant (7) [1/s]: ');
  % Environment
    m     = 95;   % input ('Moved mass (95) [kg]   : ');
    alpha = 45;   % input ('Leg-press inclination angle (45) [deg]: ');
    g     = 9.81; % input ('Gravity (9.81)  [m/s^2]: ');
  % Initial Conditions
    X0 = 0.65; % input ('Initial position (0.65) [m]: ');
    V0 = 0;    % Set zero [m/s]  
  % Time range for integration
    t0 = 0;    % input ('Start integration at (0) [s]  : ');
    te = 3;    % input (' End  integration at (3) [s]  : ');
  % ODE45 solver settings (see the ode45 Matlab documentaton)
    reltol  = 1e-8;   % ode45 option 
    maxstep = 0.005;  % ode45 option
    refine  = 4;      % ode45 option
    opts    = odeset('MaxStep',maxstep,'RelTol',reltol,'Events',@events,...
        'Refine',refine); % To show the ODE plot add: ,'OutputFcn',@odeplot
  %% Calculations
    vmax = c/a-b; % Calculate vmax of muscle
    gi = sin(alpha*pi/180)*g; % fraction of gravity relative to inclination
  % Calculate initial magnitude of activation dynamics
    gm = GX(lt,ls,kt,ks,r); % Geometrical ratio: 0 to 180 deg extension
    GX0 = spline(gm(:,1),gm(:,2),X0); % Geometrical ratio at X0
  % Check if a concentric movement is possible
    if GX0 * V0 > vmax  % Otherwise the movement is eccentic
       error('Muscle contraction velocity is greater than vmax');
    end
    if m*gi >  GX0*(c/(GX0*V0 + b) -a) % Otherwise the movement is excentic
       error('Eccentric movement');
    end    
    A0 = m*gi/(GX0*(c/(GX0*V0 + b) -a)); % The initial activation
    if A0 < 1            % Shift muscle activation to match A0 at the start
       tau = 1/S*log(1-A0);
    else                 % Assume full activiation otherwise
       tau = -1000;
    end
%% ODE solving
    [t,x] = ode45(@ODE,[t0,te],[X0,V0],opts); % Call solver
    X = x(:,1); % Position data [m]
    V = x(:,2); % Velocity data [m]
%% Output
    A = 1-exp(-S*(t-tau));        % Function of muscle activation [-,ratio]
    G = spline(gm(:,1),gm(:,2),X); % Function of geometical relations [-]
    v = G .* V;                   % Muscle contraction velocity    [m/s]
    f = (c./(v + b) - a).*A;      % Muscle force                    [N]
    F  = G .* f;                  % Exteral force                   [N]
  % Print to Command Window
    fprintf(['\nXend    [m]    =  ',num2str(X(end)),'\n']) % End position
    disp(['Vmax    [m/s]  =  ',num2str(max(V))]) % External max. velocity  
    disp(['Vend    [m/s]  =  ',num2str(V(end))]) % Ext. velocity at end       
    disp(['tend    [s]    =  ',num2str(t(end))]) % Pushing duration           
    disp(['tau     [s]    =  ',num2str(tau)])    % Shift muscle activation
    disp(['int Fdt [Ns]   =  ',num2str(trapz(t,F))]) % Impulse                
    disp(['Fmax    [N]    =  ',num2str(max(F))])     % External max. force
    disp(['int Pdt [J]    =  ',num2str(trapz(t,F.*V))]) % Kinetic energy      
    disp(['Pmax    [W]    =  ',num2str(max(F.*V))])  % External max. power
%% Nested Functions  
    function x = ODE(t,X)
    %% ODE: the right side of the ordinary differential equation to solve
    % Inputs:                   Outputs:
    %  t ... Time [s]            t      ... Time [s]
    %  X ... Position [m]        x[:,1] ... Position [m]
    %                            x[:,2] ... Velocity [m/s]
        gx1 = spline(gm(:,1),gm(:,2),X(1)); % Function of geometry at X(1)
        x(1) = X(2);                    % Substituting 2nd derivative of X
       % The right side of the equation:
        x(2) = -gi + (gx1/m)*(c/(gx1*X(2) + b) -a)*(1-exp(-S*(t -tau)));
        x = x';
    end
    function [value,isterminal,direction] = events(~,x)
    %% EVENTS: locates events during the integration of the ode and aborts
    % Integration if the condition(s) are met. An "odeset(...)" function.
    % Condition(s) for abort:
      %   01  Contraction velocity is greater than vmax
      %   02  Position data X is greater than lo+lu
      % Detect these conditions at:
        value = [vmax - x(2)*spline(gm(:,1),gm(:,2),x(1)),lt+ls-x(1)]; 
        isterminal = [1,1]; % Stop the integration
        direction  = [0,0]; % Detect always (increasing & decreasing)
    end
end
%% Subfunctions
function xgx = GX(lt,ls,kt,ks,r)
%% GX: calculates the geometrical relations
% Inputs:
%  lt ... Length of the thigh [m]
%  ls ... Length of the shank [m]
%  kt ... Length of the muscle [m]
%  ks ... Length from patella center to tibial tuberosity [m]
%  r  ... Joint radius [m]
% Outputs:
%  g[:,1] ... Vector of center of mass positions [m]
%  g[:,2] ... Vector of the gemetrical relation at each position [-]
    bet   = 0:0.01:pi/2; % Allocate beta for 0-180 deg knee extension angle
    % Calculate the corresponding knee extension angles
    sigma = 2*bet +asin((r/kt)*sin(bet)) +asin((r/ks)*sin(bet));
    sigma = sigma(sigma<pi); % Remove values greater 180 deg
    beta  = bet(1:length(sigma)); % Keep the useful part of beta only
    % Calculate the corresponding position
    X = sqrt(lt*lt + ls*ls - 2*lt*ls*cos(sigma));
    % Calculate the geometrical relation at each position and angle
    G    = ((r*sin(beta))./(lt*ls*sin(sigma))).*X;
	G(1) = 0; % Avoid division by 0 at sigma = 0
    xgx   = [X',G']; % Pack position and geometrical relation together
end