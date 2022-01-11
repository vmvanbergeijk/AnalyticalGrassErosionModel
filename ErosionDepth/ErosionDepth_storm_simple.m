% Date: 10/9/2019
% Author: Vera van Bergeijk
% 
% -------------------------------------------------------------------------
% THE FLOW VELOCITY AND THE EROSION DEPTH ALONG THE DIKE PROFILE
% -------------------------------------------------------------------------
%
% The flow velocity and the dike cover erosion along the dike profile are
% calculated. The required input parameters are:
% - the threshold flow velocity
% - number of overtopping waves
% - average overtoppin discharge
% - the storm duration
% - the dike geometry
%
% First, the volume ditribution of the overtopping waves is build based on 
% the overtopping discharge and the number of overtopping waves. Next, the
% erosion depth for every individual wave is calculated with the following
% steps:
% 1 The boundary conditions for the flow velocity model are determined
%   using the overtopping volume and the ermpirical formulas by Van der
%   Meer et al. (2010).
% 2 The flow velocity along the crest and landward slope is calculated
%   using the analytical formulas derived by Van Bergeijk et al.(2019).
% 3 The erosion depth along the profile is calculated using the analytical
%   Grass-Erosion model (Van Bergeijk et al., 2021).
%
% In the end, the total erosion depth along the profile during the storm is
% calculated by summing over all the individual erosion depths.
%
% -------------------------------------------------------------------------
% Flow velocity
% -------------------------------------------------------------------------
% The formula for the flow velocity on the horizontal parts of the profile
% Uh = 1/(f x/2 Q +1/u0) 
% with  f the bottom friction parameter
%       x the cross dike coordinate (m)
%       Q the momentary discharge, assumed to be constant (m^3/s)
%       u0 the flow velocity at the start of the horizontal part (m/s)
%
% The formula for the slope is 
% Uslope = alpha/beta + mu*exp(-3alpha*beta^2*x/cos(phi))
% with alpha = (g sin(phi))^(1/3)
%      beta  = (f/2Q)^(1/3)
%      mu    = us0-alpha/beta
% with phi the slope angle (radians) and us0 the flow velocity at the start 
% of the slope (m/s)
%
% -------------------------------------------------------------------------
% Erosion Model
% -------------------------------------------------------------------------
% The erosion depth zm is calculated using the analytical grass-erosion 
% model for U > Ut
% zm = [omega^2 U^2(x) - Ut^2] T0 Ce 
% with  omega = turbulence parameter
%       U(x)  = modelled flow velocity (m/s)
%       Ut    = threshold flow velocity (m/s)
%       T0    = overtopping period (s)
%       Ce    = inverse strength parameter (s/m)
%
% The turbulence parameter is written as
% omega = 1.5 + 5 r0 
% with the turbulence intensity r0. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Default options
set(0,'defaultaxesfontname','TimesNewRoman')
set(0,'defaulttextfontname','TimesNewRoman')
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',14)
set(0,'DefaultTextFontWeight','normal')
set(0,'DefaultAxesFontWeight','normal')

% Settings
Ut      = 10;            % threshold flow velocity (m/s)
Now     = 2500;          % number of overtopping waves   
qm      = 10/1000;       % average overtopping discharge (m3/s/m)
t_storm = 6*3600;        % storm duration (s)

% Dike geometry
f     = 0.01;         % bottom friction coefficient of grass
Bc    = 5;            % crest width (m)
phi   = atan(1/3);    % landward slope angle in radians
H     = 5;            % crest height (m)
L     = H/sin(phi);   % length of lower slope (m)
dx    = 0.01;         % spatial step (m)
x_c   = [0:dx:Bc];    % cross-dike coordinate array on crest (m)
x_s   = [0:dx:L];     % cross-dike coordinate array on slope (m) 

% Parameters flow-erosion model
g       = 9.81;                     % gravitational acceleration (m/s2)
r0      = 0.1;                      % turbulence intensity
omega   = (1.5+5.*r0);              % turbulence parameter
alpha   = (g*sin(phi)).^(1/3);      % constant in the flow velocity formula
Ce      = 1e-6;                     % inverse strength parameters (s/m)

%% Calculation of the erosion depth after the storm

% Volume distribution
a     = (0.84.*qm*t_storm)/Now; % shape parameter volume distribution
b     = 0.75;                   % shape parameter volume distribution
V_max = a*log(Now)^(1/b);       % maximum volume (m3/m)
dy    = 0.001;                  % interval wave volume vector (m3/m)
V     = [0:dy:V_max];           % volume array
P     = exp(-(V./a).^b);        % probability function of volumes

% From exceedance probability to probability that V is in delta_V
for n = 2:length(P)
    K(n) = P(n-1)-P(n); 
end
K(1) = 1-sum(K);

% Calculate flow velocity and erosion depth per wave
for j = 1:Now
    
    % Boundary condition
    Wv = randsample(V,1,true,K);   % wave volume (m3/m) random sampled from the distribution
    T0 = 3.9.*Wv.^(0.46);          % wave overtopping period [s]
    u0 = 4.5.*Wv.^(0.3);           % flow velocity at x = 0 m [m/s]
    h0 = 0.133.*Wv.^(0.5);         % layer thickness at x = 0 m [m]
    Q  =  u0.*h0;                  % overtopping discharge at x = 0 [m3/s/m] 

    % Flow velocity 
    u_c  = (f*x_c/(2*Q)+1/u0).^(-1);                                        % flow velocity on the crest
    beta  = (f/(2*Q)).^(1/3);                                               % parameter for the  flow velocity
    mu    = u_c(end)-alpha./beta;                                           % parameter for the  flow velocity
    u_s   = alpha./beta+mu.*exp(-3*alpha*beta^2.*x_s./cos(phi));            % flow velocity on the slope
    u_tot = [u_c u_s];                                                      % flow velocity along the dike profile

    % Calculate the erosion depth
    for n=1:length(u_tot)
        % only erosion when the flow velocity Utot is larger than the
        % threshold velocity Ut
        if u_tot(n)>Ut
            zm(j,n)=-Ce*T0*((omega.^2*u_tot(n)).^2-(Ut^2));       % erosion depth of one wave
        else
            zm(j,n) = 0;
        end
    end

end

% Total erosion depth along the dike profile after one storm
Etot = max(sum(zm));


