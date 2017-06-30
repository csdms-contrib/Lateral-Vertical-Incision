function LatVerInc_plot_xs_snapshot(time1,time2,time3,time4,ysnap0,ysnap1,ysnap2,ysnap3,ysnap4)

% PLOT SNAPSHOTS OF FOUR PROFILES DURING THEIR INCISION
%
% timex  : time of snapshots
% ysnapx : topography of snapshot
%
% Copyright (C) 2017 Luca C. Malatesta
% Developer can be contacted at lcmalate@ucsc.edu and
% lucamalatesta.weebly.com

if length(ysnap0)>5000  % if model width > 5000, then crop the plots
    center = length(ysnap0)/2 ;
    xL = center - 2500 ;    % left edge of plot
    xR = center + 2500 ;    % right edge of plot
else        % no crop
    xL = 1 ;               
    xR = length(ysnap0) ;
end

% SUBPLOT ONE
subplot(2,2,1)
plot(ysnap0(1,xL:xR))
hold on
plot(xL:xR,ysnap1(1,xL:xR))
plot(xL:xR,ysnap2(1,xL:xR))
plot(xL:xR,ysnap3(1,xL:xR))
plot(xL:xR,ysnap4(1,xL:xR))
legend('0',num2str(time1),num2str(time2),num2str(time3),num2str(time4),'Location','southeast')
xlim([xL xR])
title('X-section 1','fontsize',15,'FontWeight','bold')
ylabel('elevation')


% SUBPLOT TWO
subplot(2,2,2)
plot(xL:xR,ysnap0(2,xL:xR))
hold on
plot(xL:xR,ysnap1(2,xL:xR))
plot(xL:xR,ysnap2(2,xL:xR))
plot(xL:xR,ysnap3(2,xL:xR))
plot(xL:xR,ysnap4(2,xL:xR))
legend('0',num2str(time1),num2str(time2),num2str(time3),num2str(time4),'Location','southeast')
xlim([xL xR])
title('X-section 2','fontsize',15,'FontWeight','bold')


% SUBPLOT THREE
subplot(2,2,3)
plot(xL:xR,ysnap0(3,xL:xR))
hold on
plot(xL:xR,ysnap1(3,xL:xR))
plot(xL:xR,ysnap2(3,xL:xR))
plot(xL:xR,ysnap3(3,xL:xR))
plot(xL:xR,ysnap4(3,xL:xR))
legend('0',num2str(time1),num2str(time2),num2str(time3),num2str(time4),'Location','southeast')
xlim([xL xR])
title('X-section 3','fontsize',15,'FontWeight','bold')
ylabel('elevation')
xlabel('distance')


% SUBPLOT FOUR
subplot(2,2,4)
plot(xL:xR,ysnap0(4,xL:xR))
hold on
plot(xL:xR,ysnap1(4,xL:xR))
plot(xL:xR,ysnap2(4,xL:xR))
plot(xL:xR,ysnap3(4,xL:xR))
plot(xL:xR,ysnap4(4,xL:xR))
legend('0',num2str(time1),num2str(time2),num2str(time3),num2str(time4),'Location','southeast')
xlim([xL xR])
title('X-section 4','fontsize',15,'FontWeight','bold')
xlabel('distance')