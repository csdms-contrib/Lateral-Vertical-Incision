function  LatVerInc_Plot_Geom(cx,dx,IXS,W,xC,y,yC)

% Plot the final model geometry in 3D
%
% Copyright (C) 2017 Luca C. Malatesta
% Developer can be contacted at lcmalate@ucsc.edu and
% lucamalatesta.weebly.com

z     = y ;                 % for sake of clarity, y-axis becomes z-axis in 3D
NX    = length(cx);
cx    = cx.*dx ;      % cx in km
xplot = find(z(1,:));       % x-coordinates
xplot = xplot.*dx ;   % x-coordinates in km
yplot = ones(NX,length(z)); % array for y-coordinates
indLB = zeros(1,NX); indRB = zeros(1,NX);   % arrays for left and right banks

for i = 1:NX
    indLB(i) = xC(i)/dx - find( IXS(i, xC(i)/dx : -1 : 1   )~=1, 1, 'first') + 1 ;        % index of left bank (position on top of bank)
    indRB(i) = xC(i)/dx + find( IXS(i, xC(i)/dx :  1 : end )~=1, 1, 'first') - 1 ;        % index of right bank (position on top of bank)
    
    yplot(i,:) = yplot(i,:) .* cx(i) ;                 % y-coordinates
    yplot(i,:) = yplot(i,:) ;             % y-coordinates in km
    
    fill3( cat(2,xplot,[xplot(end) xplot(1)]) , cat(2,yplot(i,:),[cx(i) cx(i)]), cat(2,z(i,:),[0 0]) ,[0.9, 0.9, 0.9] )

    hold on
end

indLB = indLB  .*dx  ;  % left banks in km
indRB = indRB  .*dx ;  % right banks in km

hFP = fill3(cat(2,indLB,fliplr(indRB)) ,  cat(2,cx',fliplr(cx')) ,  cat(2,yC',fliplr(yC')),'c') ;    % floodplain plane
alpha(hFP,0.3)          % transparency

hCN = fill3(cat(1,xC-(W/2),flipud(xC+(W/2)))' ,  cat(2,cx',fliplr(cx')) ,  cat(2,yC'+0.2,fliplr(yC')+0.2),'b') ;    % channel plane
alpha(hCN,0.6)          % transparency

set(gcf,'Color', [1,1,1])    
% axis equal

% xlim([1500 2500])

% view(150+0.0*plotcount,60)        % viewpoint (AZ, EL) in degrees
view(150,78) 

xlabel('distance across valley','Fontsize',13)
ylabel('distance along fan','Fontsize',13)
zlabel('elevation','Fontsize',13)
title('Entrenchment geometry','fontsize',15,'FontWeight','bold')

legend([hFP,hCN],'active floodplain','channel','Location','NorthEast')

set(gca,'XDir','Reverse')

hold off