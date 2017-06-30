function LatVerInc_plot_entrenchmentGeom(iEC,talusoutC,WC,WFPC)

% PLOT CHARACTERISTICS DESCRIBING ENTRENCHMENT OF RIVER
% SUBPLOT ONE compares channel width, floodplain width, and talus extent in
% the 4 upstream cross-sections
% SUBPLOT TWO compares the entrenchment indices of the 4 upstream
% cross-sections
%
% Copyright (C) 2017 Luca C. Malatesta
% Developer can be contacted at lcmalate@ucsc.edu and
% lucamalatesta.weebly.com

colours = ([ 0.1 0.75 0.2 ; 0.1 0.1 0.95 ; 0.8 0.3 0 ; 0.5 0 0.5 ]) ;   % vector of 4 colours

subplot(3,1,(1:2))
for i=1:size(iEC,1)
   plot(1:20:size(WFPC,2),WFPC(i,1:20:end),'Color',colours(i,:))
   hold on
   plot(1:20:size(WC,2),WC(i,1:20:end),'Color',colours(i,:))
   plot(1:20:size(talusoutC,2),talusoutC(i,1:20:end),'Color',colours(i,:))
end
xlabel('time','Fontsize',13)
ylabel('width','Fontsize',13)    
title('channel, floodplain, and talus widths','fontsize',15,'FontWeight','bold')

subplot(3,1,3)
for i=1:size(iEC,1)
    plot(1:20:size(iEC,2),iEC(i,1:20:end),'Color',colours(i,:))
    hold on
end
legend ('1','2','3','4','Location','southeast')
ylabel('Entrenchment index')
xlabel('time')
title('Funnel index','fontsize',15,'FontWeight','bold')