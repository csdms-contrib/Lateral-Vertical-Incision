function [A, AL, AR, iE, I1, I2, Itot, IXS, Qs, S, W, WFP, dCout, dir, dLb, dRb, indC, indB, partition, talusout, ucut, xC, y, yC, way] = ...
    LatVerInc_Fn (A, Ai, Aeq, AWnorm, IXS, Fx, Hh, S, Si, W, WFP,cx, dx, hFP, kA, lambda, phi, u, v, w, xC, y, yC)


% this public code corresponds to development version 3_4
%
% Copyright (C) 2017 Luca C. Malatesta
% Developer can be contacted at lcmalate@ucsc.edu and
% lucamalatesta.weebly.com
%
% RANDOM WALK CHANNEL INCISION WITH BANK FEEDBACKS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% go left or right 50/50 chance,
% BIG CHANGE EXPLAIN NEW EROSION MECHANISM HERE
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% A      = effective transport capacity (area)
% Ai     = initial transport capacity
% Aeq    = transport capacity at steady-state (= sediment transport from upstream)
% AWnorm = factor normalizing AW to max=1, for A(width)
% iE     = entrenchment index
% IXS    = cross section index
% S      = slope
% Si     = initial slope
% W      = channel width
% WFP    = floodplain width
% cx     = along-stream position of cross-sections
% dt     = timestep
% dx     = cell size
% phi    = angle of repose
% u      = exponent for falling limb
% v      = exponent for rising limb
% w      = offset coefficient for AW to peak at Wo (optimal width)
% xC     = x-coordinates of channel centre
% y      = cross-section topography
% yC     = channel elevation
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% REACH CONFIGURATION
%    ____         ____
%        |   .---'
%         '--        ___            ____                         __ yC(1)
%                       '--      __|
%                          '--__|      _____             _____   __ yC(i)
%                                           |_     _...''
%                                             |___|              __ yC(NX)
% 
%         cx(1)     ...      cx(i)     ...   cx(NX)  
%
%  S(1)=(yC(1)-yC(2))/DX           S(NX)=(yC(NX-1)-yC(NX))/DX
%                  S(i)=(yC(i-1)-yC(i+1))/2DX
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% CHANNEL GEOMETRY
%             
%            _dLb__ _____W_____ __dRb__
%           '      '           '       '
%   ________                            ________ 
%           |                          |
%           '------_____________-------'        _ _ _ yC
%           |      |     |     |       |
%         indLB  indLC   xC  indRC   indRB
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% TALUS GEOMETRY 
%      ___                   ____
%         |                 |
%         |                 |           0 = substrate/bedrock
%          `-._             |           1 = floodplain
%              ``-.___......|           2 = channel (not yet implemented)
%         '-------'                     3 = talus
%            txL
% IXS: 0000333333332221111111000
%
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% timestep incision or all at once?
% stepwise = 0 ;  % 1: timesep, 0: all at once

NX = length(A) ;                        % number of cross section

%  A R R A Y   F A C T O R Y  _____________________________________________
      
indB  = zeros(NX,2) ;   indC  = zeros(NX,2) ;             % index of banks and cliff's edges (left and right)
txL   = zeros(NX,1) ;   txR   = zeros(NX,1) ;             % horizontal extent of left and right taluses
partition = zeros(NX,1) ;                                       % fraction of the 'way' covered until the bank
AL = zeros(NX,1); AR = zeros(NX,1);       % array for respective lateral and downward excavation
a  = zeros(NX,1) ;                                              % erosion product
I2    = zeros(NX,1) ;                                           % incision during stage 2
At    = zeros(NX,1) ;                                           % area of talus
Estep = zeros(NX,1) ;                                           % number of steps until lateral erosion stops
Afp = zeros(NX,1) ;        % array for area removed by "free" lateral migration in floodplain, to be moved to the other side of channel 
RedTwoStandingBy  = 0 ; RedFiveStandingBy = 0 ;                 % flag for left and right erosion fully completed
ucut = zeros(NX,1) ;
% T = 1 ;     % total time (non dimensional), made of a # of dt steps that are defined in initialisation code


%  R A N D O M   M A C H I N E  ___________________________________________

dir = randi(2,NX,1) ;                           % randomly issue dir = 1 or 2 to move left or right
way = round(random('Normal',4.*W,3.*W,NX,1)) ;   % random pathlength of the lateral migration as a normal distribution around a channel width

way(way<1) = 1 ;   % get rid of negative and zero pathlength
% way=way.*0;

% I N D E X  ______________________________________________________________

for i = 1:NX        % loop through the profiles
    indB(i,1) = xC(i)/dx - find( IXS(i, xC(i)/dx : -1 : 1   )~=1, 1, 'first') + 1 ;      % index of left bank (position on top of bank), starting from channel
    indB(i,2) = xC(i)/dx + find( IXS(i, xC(i)/dx :  1 : end )~=1, 1, 'first') - 1 ;      % index of right bank (position on top of bank), starting from channel
    indC(i,1) = xC(i)/dx - find( IXS(i, xC(i)/dx : -1 : 1   )==0, 1, 'first') + 1 ;      % index of left cliff's edge (past the talus), starting from channel
    indC(i,2) = xC(i)/dx + find( IXS(i, xC(i)/dx :  1 : end )==0, 1, 'first') - 1 ;      % index of right cliff's edge (past the talus), starting from channel
end

indLC = round(xC./dx - ceil(W./2)./dx) ;        % index of left channel extremity, in water
indRC = round(xC./dx + ceil(W./2)./dx) ;        % index of right channel extremity, in water

dLb = round((indLC - indB(:,1)) .* dx) ; % distance between left  bank (i.e. above floodplain) and left  edge of channel
dRb = round((indB(:,2) - indRC) .* dx) ; % distance between right bank (i.e. above floodplain) and right edge of channel

iLBe = find(dir==1 & way<=dLb) ;        % index for leftward erosion NOT reaching the bank
iRBe = find(dir==2 & way<=dRb) ;        % index for rightward erosion NOT reaching the bank
iLBa = find(dir==1 & way>dLb) ;         % index for leftward erosion reaching the bank
iRBa = find(dir==2 & way>dRb) ;         % index for rightward erosion reaching the bank

% the case of Zero erosion   .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
way(dir==1 & way>dLb & A==0) = dLb(dir==1 & way>dLb & A==0) ;    %  leftward path reaching the bank without erosion capacity is limited to floodplain (way<=dLb)
way(dir==2 & way>dRb & A==0) = dRb(dir==2 & way>dRb & A==0) ;    % rightward path reaching the bank without erosion capacity is limited to floodplain (way<=dRb)

% xCold = xC ;     % save position of previous channel center
yCold = yC ;     % save values of previous channel elevation
Wold  = W ;      % save values of width to correct channel position below


% U N D E R C U T  ________________________________________________________

ucut(dir==1) = round( (y( sub2ind(size(y), find(dir==1), indC(dir==1,1)) ) - yC(dir==1)) ./ (2*Hh^2*tan(phi)) ) ;  % set number of steps of undercut before ...
ucut(dir==2) = round( (y( sub2ind(size(y), find(dir==2), indC(dir==2,2)) ) - yC(dir==2)) ./ (2*Hh^2*tan(phi)) ) ;  % failure as a function of cliff height

ucutmax = 20 ;
ucut(ucut>ucutmax) = ucutmax ;    % limit the maximum size of the undercut
ucut(ucut<2 ) = 2 ;


% F L O O D P L A I N   B U S I N E S S  __________________________________

% compute the area of floodplain material moved by left-/rightward channel migration, to be distributed on the other side of the channel
for i=1:length(iLBe)
    Afp(iLBe(i)) = sum( y(iLBe(i), indLC(iLBe(i))-1-floor(way(iLBe(i))/dx) : indLC(iLBe(i))-1) - yCold(iLBe(i)) * dx )  ; 
end          % using 'floor' here to avoid sampling the cliff by wrong rounding
for i=1:length(iRBe)
    Afp(iRBe(i)) = sum( y(iRBe(i), indRC(iRBe(i))+1 : indRC(iRBe(i))+1+floor(way(iRBe(i))/dx)) - yCold(iRBe(i)) * dx ) ;
end
for i = 1:length(iLBa)
    Afp(iLBa(i)) = sum( ( y( iLBa(i), indB(iLBa(i),1) + 1 : indLC(iLBa(i)) - 1 )  - yCold(iLBa(i))) * dx ) ;
end
for i = 1:length(iRBa)
    Afp(iRBa(i)) = sum( ( y( iRBa(i), indRC(iRBa(i)) + 1  : indB(iRBa(i),2) - 1 ) - yCold(iRBa(i))) * dx ) ;
end


% L A T E R A L   E R O S I O N  __________________________________________

% areas to be eroded (ideally) off the cliff
AL(iLBa) = A(iLBa) .* ( (lambda-1) .* ((dLb(iLBa)./way(iLBa))-1) ) ;
AR(iRBa) = A(iRBa) .* ( (lambda-1) .* ((dRb(iRBa)./way(iRBa))-1) ) ;

for j = 1:1000    % start eroding laterally until the ceiling value (AL or AR) is reached 
    
    if ~any(a(iLBa)<AL(iLBa))            % stop left lateral erosion when all cross sections have used their full erosion capacity
        RedTwoStandingBy = true ;       % and raise a flag
    end
    if ~any(a(iRBa)<AR(iRBa))           % stop right lateral erosion when all cross sections have used their full erosion capacity
        RedFiveStandingBy = true ;      % and raise a flag
    end
    
    % if Red Five is standing by, and Red Two is standing by, then lock S-foils in attack position
    if RedTwoStandingBy == true && RedFiveStandingBy == true  ...
            || isempty(iRBa) && RedTwoStandingBy == true ...
            || isempty(iLBa) && RedFiveStandingBy == true 
        break                           % break the loop
    end
    
    iLBaE  = find(dir==1 & way>dLb & a<AL) ;    % index for leftward erosion reaching the bank and not yet reaching full erosion capacity
    iRBaE = find(dir==2 & way>dRb & a<AR) ;    % index for rightward erosion reaching the bank and not yet reaching full erosion capacity
    
    kE = zeros(NX,1);   % reset switch for lateral erosion
    kE (cat(1,iLBaE, iRBaE)) = 1 ;   % turn on switch for active erosion arrays

    % compute erosion product of new increment and add it to what's been eroded so far:
    if ~isempty(iLBaE)
        for i = 1 : length(iLBa)
            a(iLBa(i))  = a(iLBa(i))  + kE(iLBa(i)) * dx * ...
                ( y( iLBa(i),   indB(iLBa(i),1)  - Estep(iLBa(i))  ) - yCold(iLBa(i))  )  ;
        end
    end
    
    % compute erosion product of new increment and add it to what's been eroded so far:
    if ~isempty(iRBaE)
        for i = 1 : length(iRBa)
            a(iRBa(i)) = a(iRBa(i)) + kE(iRBa(i)) * dx * ...
                ( y( iRBa(i),  indB(iRBa(i),2) + Estep(iRBa(i)) )  - yCold(iRBa(i)) ) ;
        end
    end
    
    Estep(iLBaE) = Estep(iLBaE) + 1 ;    % steps of left lateral erosion
    Estep(iRBaE) = Estep(iRBaE) + 1 ;    % steps of right lateral erosion
   
end


dCi   = Estep ;                     % ideal lateral erosion
dCeff = ceil(dCi./ucut).*ucut ;     % effective lateral undercut

dCout = dCi; dCout(5:end,:)=0;

% calculate Area eroded between dCi and dCeff

for i=1:length(iLBa)          % loop through leftward cross sections
    if dCeff(iLBa(i))>dCi(iLBa(i))
        At(iLBa(i)) = sum( y(iLBa(i), indB(iLBa(i),1)-dCeff(iLBa(i)) : indB(iLBa(i),1)-dCi(iLBa(i))-1) - yCold(iLBa(i)) ) ;
    end
end

for i=1:length(iRBa)          % loop through rightward cross sections
    if dCeff(iRBa(i))>dCi(iRBa(i))
        At(iRBa(i)) = sum( y(iRBa(i), indB(iRBa(i),2)+dCi(iRBa(i))+1 : indB(iRBa(i),2)+dCeff(iRBa(i))) - yCold(iRBa(i)) ) ;
    end
end


% I N C I S I O N  ________________________________________________________ 

I1 = A ./(way+W) ;        % stage 1 vertical incision distributed over channel migration path

% Stage 2 vertical incision when in contact with left bank
I2(dir==1 & way>=dLb) = lambda .* (A(dir==1 & way>=dLb)./(W(dir==1 & way>=dLb)+dCeff(dir==1 & way>=dLb))) .* ( 1 - (dLb(dir==1 & way>=dLb)./way(dir==1 & way>=dLb)) ) ; 
% Stage 2 vertical incision when in contact with right bank
I2(dir==2 & way>=dRb) = lambda .* (A(dir==2 & way>=dRb)./(W(dir==2 & way>=dRb)+dCeff(dir==2 & way>=dRb))) .* ( 1 - (dRb(dir==2 & way>=dRb)./way(dir==2 & way>=dRb)) ) ; 

Itot = I1 + I2 ;    % total incision
yC = yC-Itot;       % subtract directly the value of total incision without timestep subdivision


% T A L U S  ______________________________________________________________

txL(iLBa) = round( sqrt((2.* At(iLBa)) ./ tan(phi)) ) ;     % compute horizontal extent of left talus
txR(iRBa) = round( sqrt((2.* At(iRBa)) ./ tan(phi)) ) ;     % compute horizontal extent of right talus

talusout = txL + txR ;  % all taluses for export

% C H A N N E L ___________________________________________________________

xC(iLBe) = xC(iLBe) - way(iLBe)  ;                      % update channel center for leftward  movement in floodplain
xC(iRBe) = xC(iRBe) + way(iRBe) ;                      % update channel center for rightward movement in floodplain
space=zeros(NX,1); space(txL>0)=1;  % offset the channel by one unit to leave room for the talus (if present)
xC(iLBa) = xC(iLBa) - dLb(iLBa) - dCeff(iLBa) + txL(iLBa) + space(iLBa)  ;    % update the position of the channel center
space=zeros(NX,1); space(txR>0)=1;  % offset the channel by one unit to leave room for the talus (if present)
xC(iRBa) = xC(iRBa) + dRb(iRBa) + dCeff(iRBa) - txR(iRBa) - space(iRBa) ;    % update the position of the channel center

xC = round(xC) ;


% C H A N N E L   N A R R O W I N G  ______________________________________

for i = 1:length(iLBa)
    if WFP(iLBa(i))-txL(iLBa(i))-1 < W(iLBa(i))    % narrow channel width if constrained by talus: WFP - txL < W
        W(iLBa(i))  = round(WFP(iLBa(i)) - txL(iLBa(i))) ;              % define channel width as the floodplain width
        W(W<50)=50 ;                                                                   % minimal size (shouldn't be directly constrained in future)
        xC(iLBa(i)) = xC(iLBa(i)) - (Wold(iLBa(i))-W(iLBa(i)))/2 ; % replace channel center with respect to narrowed channel
        xC(iLBa(i)) = round(xC(iLBa(i))) ;                                   
    end  
end

for i = 1:length(iRBa)
    if WFP(iRBa(i))-txR(iRBa(i))-1 < W(iRBa(i))   % narrow channel width if constrained by talus: WFP - txR < W
        W(iRBa(i))  = round(WFP(iRBa(i)) - txR(iRBa(i))) ;               % define channel width as the floodplain width
        W(W<50)=50 ;                                                                       % minimal size (shouldn't be directly constrained in future)
        xC(iRBa(i)) = xC(iRBa(i)) - (Wold(iRBa(i))-W(iRBa(i)))/2 ; % replace chanel center with respect to narrowed channel
        xC(iRBa(i)) = round(xC(iRBa(i))) ;
    end   
end


% N E W  T O P O G R A P H Y   ____________________________________________

% b e h i n d   c h a n n e l  ............................................

for i=1:length(iLBe)
    y(iLBe(i),  indRC(iLBe(i))-way(iLBe(i)) : indRC(iLBe(i))+1)  = ...   % new floodplain geometry after distribution of erosion  across 'way'
        yCold(iLBe(i)) + Afp(iLBe(i))/way(iLBe(i)) - I1(iLBe(i)) ;      % (I1) and redistribution of Afp on other side of channel (Afp/way)
end
for i=1:length(iRBe)
    y(iRBe(i), indLC(iRBe(i))-1 : indLC(iRBe(i))+way(iRBe(i))) = ... % new floodplain geometry after distribution of erosion  across 'way'
        yCold(iRBe(i)) + Afp(iRBe(i))/way(iRBe(i)) - I1(iRBe(i)) ;  % (I1) and redistribution of Afp on other side of channel (Afp/way)
end
for i=1:length(iLBa)       
    y(iLBa(i),  indRC(iLBa(i))-dLb(iLBa(i))-dCeff(iLBa(i)) : indRC(iLBa(i))+1)  = ... % new floodplain geometry after distribution of erosion  across 'way'
        yCold(iLBa(i)) + Afp(iLBa(i))/way(iLBa(i)) - I1(iLBa(i)) ;                       % (I1) and redistribution of Afp on other side of channel (Afp/way)
end
for i=1:length(iRBa)      
    y(iRBa(i), indLC(iRBa(i))-1 : indLC(iRBa(i))+dRb(iRBa(i))+dCeff(iRBa(i))) = ... % new floodplain geometry after distribution of erosion  across 'way' 
        yCold(iRBa(i)) + Afp(iRBa(i))/way(iRBa(i)) - I1(iRBa(i)) ;                      % (I1) and redistribution of Afp on other side of channel (Afp/way)
end

% c h a n n e l  ..........................................................

for i=1:NX    % update the topography with the river [ xC - W/2 ] : [ xC + W/2 ]
    y(i,xC(i)-ceil(W(i)/2) : xC(i)+ceil(W(i)/2)) = yC(i) ;
end

% t a l u s  ..............................................................

for i = 1:length(iLBa)
    if round(txL(i))>0      
        y (iLBa(i),...     % update the topography to include left talus
                xC(iLBa(i)) - ceil(W(iLBa(i))/2) - txL(iLBa(i)) ...
              : xC(iLBa(i)) - ceil(W(iLBa(i))/2) - 1 ) ...
                                = -tan(phi) * ( 1 : txL(iLBa(i))/dx ) + yCold(iLBa(i)) - I1(iLBa(i)) + tan(phi) * txL(iLBa(i)) ;                 
        IXS(iLBa(i),...    % update the index and add the new talus
                xC(iLBa(i)) - ceil(W(iLBa(i))/2) - txL(iLBa(i)) ...
              : xC(iLBa(i)) - ceil(W(iLBa(i))/2) - 1 ) = 3 ;
    end
    y (iLBa(i), y(iLBa(i),:)<yC(iLBa(i))) = yC(iLBa(i)) ;  % correct talus toe to channel height if too low
end

for i = 1:length(iRBa)
    if round(txR(iRBa(i)))>0
        y (iRBa(i),...    % update the topography to include right talus
                xC(iRBa(i)) + ceil(W(iRBa(i))/2) + 1 ...
              : xC(iRBa(i)) + ceil(W(iRBa(i))/2) + txR(iRBa(i))  ) ...
                                = tan(phi) * ( 1 : txR(iRBa(i))/dx ) + yCold(iRBa(i)) - I1(iRBa(i)) ;
        IXS(iRBa(i),...   % update the index and add the new talus
                xC(iRBa(i)) + ceil(W(iRBa(i))/2) + 1 ...
              : xC(iRBa(i)) + ceil(W(iRBa(i))/2) + txR(iRBa(i))  ) = 3 ;
    end
    y (iRBa(i), y(iRBa(i),:)<yC(iRBa(i))) = yC(iRBa(i));  % correct talus toe to channel height if too lo
end


% C E L L   T Y P E   I N D E X  __________________________________________
% 0) bedrock, 1) floodplain, 2) channel, 3) talus

IXS(IXS==1) = 0 ; % reset all floodplain in cross section index to rock before identifying floodplain again

for i = 1:NX
    IXS (i ,   y(i,:)<=yC(i)+hFP) = 1 ;                        % floodplain = 1, rest = 0, with a unique FP definition
    indB(i,1)  = xC(i)/dx - find( IXS(i, xC(i)/dx : -1 : 1   )~=1, 1, 'first') + 1 ;      % index of left bank (position on top of bank), starting from channel
    indB(i,2)  = xC(i)/dx + find( IXS(i, xC(i)/dx :  1 : end )~=1, 1, 'first') - 1 ;      % index of right bank (position on top of bank), starting from channel
    WFP(i) = (indB(i,2)-1)*dx - (indB(i,1)+1)*dx ;                                         % width of floodplain from bank to bank
end


% F I N A L   V A L U E S   F O R   E X P O R T  __________________________

% s l o p e  ..............................................................  

Dcx = diff(cx) ;                   % distance between each X-section
dev = diff(xC) ;                   % lateral offset of channel between each X-section
DxC = sqrt(Dcx.^2+dev.^2) ;        % horizontal distance between each channel centers

S(1:end-1) =  - ( diff(yC(1:end) ) ./ (DxC(1:end)  .*dx) ) ;        % slope of the cells, upstream scheme
S(end)     =   ( yC(end-1) - 0 ) / ((DxC(end)+Fx-cx(end))*dx) ;     % slope of the last (downstream) cell, downstream scheme


% T R A N S P O R T    C A P A C I T Y ____________________________________

AW  = AWnorm .* (W.^v)./(w+W).^u ; % modulates transport capacity with respect to width, according to some gamma distribution type rule
Qs = AW .* ((S.^1.5)./(Si^1.5)) .* Ai ;     % transport capacity linearly depending on slope and according to W with fn. AW
% Qs = AW .* ((S)./(Si)) .* Ai ;

%   /!\
Qs(Qs<Aeq)=Aeq; %crude equilibrium control

DQs   = diff([Aeq; Qs]) ;
A = DQs ;
A = A .* kA ;     % TEST TO SEE EFFECTS OF REDUCING EROSION
A=round(A);
A(A<0) = 0 ;            % FIND A SOLUTION FOR A=0


% A S P E C T   R A T I O  ________________________________________________

    
ixs = 1:4 ;   % pick cross section

x = 1 : length(y(1,:)) ;    % vector of the positions on the x axis
iE = zeros(1,length(ixs)) ; % array for entrenchment index

for i = 1:length(ixs)

    % sliding threshold from base to top of cross section avoiding the extreme points
    yref = (ceil(yC(ixs(i))) + 1 : floor(y(ixs(i),1))-4.5)' ;

    Y    = repmat(y(ixs(i),:), length(yref), 1) ;  % create repeat matrix (lines)  of elevation at ixs 
    Yref = repmat(yref, 1, length(x)) ;            % create repeat matrix (column) of yref of same size as Y

    H = zeros(size(Y)) ;
    H(Yref-Y>0) = 1 ;        % when Yref>Y, H=1, when Yref<Y, H=0
    Iw = sum(H,2).*dx ;      % vector of width under the sliding threshold
    Iw = cumsum(Iw) ;        % cumulative width
    Iw = Iw ./ max(Iw) ;     % normalize Iw by the total area

    iE(i) = 2 * trapz(Iw ) / length(yref) ; % entrenchment index (0 - 1)
    
    clear yref Y Yref H Iw

end



end