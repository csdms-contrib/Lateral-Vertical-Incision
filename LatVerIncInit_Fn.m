function [A, Ai, Aeq, AWnorm, Fx, Hh, S, Si, W, WFP, cx, cz, dx, IXS, hFP, k, kA, lambda, maxi, phi, u, v, w, x, xC, y, yC] = LatVerIncInit_Fn

% Initial parameters for the model of channel incision with bank feedback.
% corresponds to development v. 3_0
% Fn is designed for interdependent cross section models that co-evolve
% along the same reach. It uses arrays for the parameters individual to
% each model, the shared parameters (e.g. physical laws) are unique values
%
% Copyright (C) 2017 Luca C. Malatesta
% Developer can be contacted at lcmalate@ucsc.edu and
% lucamalatesta.weebly.com
%  
% /!\ initial geometries are identical for every cross-section at the
% moment, i.e. vectors where only first value is called (e.g. W(1))
%
% A:        Cross section area removed by erosion
% Ai:       Erosion capacity at channel geometry of initial timestep
% Aeq:      Part of erosion capacity used for sediment bypass
% AWnorm:   Normalization factor for AW (max=1) 
% Fx:       Position of fan toe
% Hh:       Ratio of between cliff and talus heights (limits talus size)   
% S:        Slope
% Si:       Initial slope
% W:        Channel width
% WFP:      Floodlpain width
% cx:       x-position of cross section
% cz:       z-position of cross section
% dx:       cell width along x-axis
% IXS:      index of floodplain vs. rest
% hFP:      relief of floodplain (below is floodplain)
% k:        number of time steps
% kA:       factor of erosivity
% lambda:   distriubtion factor of bank vs. floodplain erosion
% maxi:     maximum incision (numerical stability)
% phi:      angle of repose
% u:        exponent for width function
% v:        exponent for width function
% w:        offset factor for width function
% x:        x-topography
% xC:       x-position of channel
% y:        y-topography
% yC:       y-position of channel




% geometry ________________________________________________________________

% fan - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% the longitudinal profile of the fan is explicit to compute evolution of
% the channel slope during incision


%     |
%     |'-._
%  cz_| _  _'-._
%     |         |'-._
%     |               '-._ _____ S
%     |         |          '-._ | 
%     |                         '-._
%     |         |                    '-._
%     |_______________________________________
%               |                           |
%               cx                          Fx
 

NX = 20 ;           % number of cross-sections
Fx = 30000 ;        % position of the fan toe
kA = 0.1 ;          % factor of erosivity
hFP= 1 ;            % relief of floodplain

cx    = ones(NX,1) ;
cx(:) = (1:NX) * (Fx/(NX+1)) ;       % positions of the cross-sections

Si = 0.02  ;        % initial slope
cz = ones(NX,1) .* (ones(NX,1).*Fx-cx).*Si; % initial elevation of the profile


% CROSS SECTION - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

width = 15000 ;     % width of cross-section
dx = 1;             % cell width
Nx = width/dx ;     % number of cells in x-direction

W   = ones(NX,1) .* 600 ;       % channel width
Wo  = 100 ;                     % optimal width
S   = zeros(NX,1) + (-Si) ;     % initial slope
Ai  = ones(NX,1) .* 4000 ;      % initial excavated surface area for each step, L^2
Aeq = 1000 ;                    % value of excavation potential used to implicitly transport upstream sediment flux
maxi = ones(NX,1) .* 401 ;      % maximum incision   <= needs to be modified to adapt maximal depth geometry

Hh = 4 ;        % ratio between cliff height and maximum talus height to determine ucut value
lambda = 0.7 ;  % factor of incision fraction when river erodes bank compared to pure floodplain.


% TIME  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

k = 3000 ;      % number of iterations


% D E F I N E   A R R A Y S  ______________________________________________

x = zeros(NX,Nx);
y = zeros(NX,Nx);

for i=1:NX
    x(i,:) = 1:1:Nx ;           % X-topography
    y(i,:) = cz(i) ;            % Y-topography
end


% initial model ___________________________________________________________

xC  = ones(NX,1) .* round(width/2) ;  % X-position of channel center
yC  = cz-5 ;                          % Y-position of channel center

Rx = zeros(NX,round(W(1)/dx)+1) ;
Ry = zeros(NX,round(W(1)/dx)+1) ;
for i = 1:NX
    Rx(i,:) = (xC(i) - W(i)/2) : 1 : (xC(i) + W(i)/2) ;     % X-river
    Ry(i,:) = zeros(1,round(W(1)/dx)+1) + yC(i) ;           % Y-river
end

% width function - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% scales transport capacity A with width W: 0 at origin, peaking (=1) at optimal
% width W=Wo, goest to 0 at infinite

u = 0.9 ;                         % exponent for falling limb
v = 0.6 ;                         % exponent for rising limb
w = Wo * ((u-v)/v) ;              % offset coefficient to place AW peak at Wo
AWnorm = ((Wo^v)/((w+Wo)^u))^-1 ; % normalization of AW so that max=1

AW = AWnorm * (W(1)^v)/(w+W(1))^u ;   % transport capacity as a function of channel width

A = zeros(NX,1) + ((-S./Si) .* AW .* Ai - Aeq) ; % transport capacity available for erosion
A = 0.7.*A ;
A(A<0) = 0 ;
% (linearly depending on slope and according to W with fn. AW)
% after bypass load (AW) is transported

% initial topography - - - - - - - - - - - - - - - - - - - - - - - - - - -

IXS   = zeros(NX,Nx) ;
indLB = zeros(NX,1) ; indRB = zeros(NX,1) ;
WFP   = zeros(NX,1) ;

phi = 30*pi/180; % angle of repose [rad]

for i = 1:NX
    y(i,xC(i)-3*W(i) : 1 : xC(i)+3*W(i)) = yC(i) + 0.1 ;  % floodplain width
    y(i,xC(i)-W(i)/2 : 1 : xC(i)+W(i)/2) = yC(i) ;
    IXS ( i ,  y(i,:)>=yC(i)  &  y(i,:)<=yC(i)+1.5/W(i) ) = 1 ;      % floodplain = 1, rest = 0
    indLB(i) = xC(i)/dx - find( IXS(i, xC(i)/dx : -1 : 1   )==1, 1, 'last') + 1 ;      % index of left bank (position on top of bank), starting from channel
    indRB(i) = xC(i)/dx + find( IXS(i, xC(i)/dx :  1 : end )==1, 1, 'last') - 1 ;      % index of right bank (position on top of bank), starting from channel
    WFP(i) = indRB(i) -1 - indLB(i) + 1 ;
end


end