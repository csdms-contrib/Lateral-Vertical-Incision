%
%
% ------___________                                 _________________------
% L A T E R A L -  \     <--o-->       ____________/ - - - - - - - - - - -
%  V S . - - - - -  \__     |        _/ - - - - - - - - - - - - - - - - - -
% - - V E R T I C A L  \    |    ___/  - -   Malatesta, Prancevic, Avouac
%  - - I N C I S I O N  \   V   / - - - -    JGR Earth Surface, 2017
% - - - - - - - - - - -  \_____/ - - - - -   doi:10.1002/2015JF003797
%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%                   this public code corresponds to development version 2_0
%
% Copyright (C) 2017 Luca C. Malatesta
% Developer can be contacted at lcmalate@ucsc.edu and
% Luca C. Malatesta
% Earth & Planetary Sciences
% University of California Santa Cruz
% 1156 High St., Santa Cruz, CA 95064
% lucamalatesta.weebly.com
%
%
% A model to explore how increasingly tall valley walls constrain the
% river lateral erosion and promote vertical incision. Each run is unique
% as a random walk controls the lateral migration of the channel. To store
% and compare repeated runs with identical parameters, there is a built in
% system to save the results of each run.
% The code is organized around the present master file and the necessary
% functions are:
% - LatVerIncInit_Fn.m    initializes all the variables, change parameters
%                         directly in that function
% - LatVerInc_Fn.m   is the looped function containing the incision step
%                    and the random walk
% - LatVerInc_Plot_Geom.m    plots the final geometry in 3D
% - LatVerInc_plot_xs_snapshot.m   plots one cross section at four diferent
%                                  times during the run
% - LatVerInc_plot_entrenchmentGeom.m   plots the characteristics of
%                                       entrenchment of the run
%
% _________________________________________________________________________
% 
% C O P Y R I G H T
%
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 2 of the License, or (at your 
% option) any later version.
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
% You should have received a copy of the GNU General Public License along 
% with this program; if not, write to the Free Software Foundation, Inc., 
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
% _________________________________________________________________________

clear ;

runno = 1 ;     % run number for saving file
nruns = 1 ;     % number of individual runs with same parameters:
                % choose 1 for plots at the end of the run
                % choose 9 for 3x3 comparison plots

for l=1:nruns

% I N I T I A L I S A T I O N _____________________________________________

% call initialization function, check the function for details on variables
[A, Ai, Aeq, AWnorm, Fx, Hh, S, Si, W, WFP, cx, ~, dx, IXS, hFP, k, ...
    kA, lambda, maxi, phi, u, v, w, ~, xC, y, yC] = LatVerIncInit_Fn ;

ysnap0 = y(1:4,:) ;     % record the initial geometry for 1st snapshot

% Create arrays for plotting at the end of single run:
elevC = zeros(length(yC),k) ;  % array for channel elevations through time
iEC   = zeros(4,k) ;           % array for entrenchment index
IeffC = zeros(length(yC),k) ;  % array for total net incision through time
I1C   = zeros(length(yC),k) ;  % array for incision without wall contact
I2C   = zeros(length(yC),k) ;  % array for incision while against wall
dCC   = zeros(length(yC),k) ;  % array for lateral erosion if no wall fdbck
WFPC  = zeros(length(yC),k) ;  % array for width of the floodplain
WC    = zeros(length(yC),k) ;  % array for channel width
SC    = zeros(length(yC),k) ;  % array for channel slope
talusoutC =zeros(length(yC),k);% array for horizontal extent of L&R taluses

snapshot = 1 ;  % counter for the geometry snapshots

% M O D E L _______________________________________________________________

for i = 1:k     % main loop
   
   % run the main function at each time step, check function for variables
   [A, AL, AR, iE, I1, I2, Itot, IXS, Qs, S, W, WFP, dC, dir, dLb, dRb,...
       indC, indB, partition, talusout, ucut, xC, y, yC, way] = ...
    LatVerInc_Fn (A, Ai, Aeq, AWnorm, IXS, Fx, Hh, S, Si, W, WFP,...
       cx, dx, hFP, kA, lambda, phi, u, v, w, xC, y, yC) ;
   
    % store values for unique plot at the end:
    elevC(:,i) = yC ;
    iEC(:,i)   = iE ;
    IeffC(:,i) = Itot;
    I2C(:,i)   = I2 ;
    I1C(:,i)   = I1;
    dCC(:,i)   = dC;
    WFPC(:,i)  = WFP ;
    WC(:,i)    = W ;
    SC(:,i)    = S ;
    talusoutC(:,i) = talusout ;


    % SAVE SNAPSHOTS OF THE RUN - - - - - - - - - - - - - - - - - - - - - -

    if mod(i, round(k/4)) == 0      % opens 4 times over entire run
        if snapshot == 1
            ysnap1 = y(1:4,:) ;   time1  = i ;
        elseif snapshot == 2
            ysnap2 = y(1:4,:) ;   time2  = i ;
        elseif snapshot == 3
            ysnap3 = y(1:4,:) ;   time3  = i ;
        elseif snapshot == 4
            ysnap4 = y(1:4,:) ;   time4  = i ;
        end    
        snapshot = snapshot+1 ;
    end

end

% F I N A L    P L O T S   ________________________________________________

if nruns==1
%---
    figure(1); clf; % plot the end geometry of the run
    LatVerInc_Plot_Geom(cx,dx,IXS,W,xC,y,yC)
%---
    figure(2); clf; % plot 4 snapshots of the 2nd x-section from top
    LatVerInc_plot_xs_snapshot(time1,time2,time3,time4,ysnap0,ysnap1,...
                               ysnap3,ysnap2,ysnap4)
%---
    figure(3); clf;     % plot the elevation of the channel through time
    for i=1:length(yC)
       plot(1:20:k,elevC(i,1:20:k),'k')
       hold on
    end
    xlabel('time','Fontsize',13)
    ylabel('elevation','Fontsize',13)    
    title('Elevation of channels through time','fontsize',15,...
          'FontWeight','bold')
%---
    figure(4); clf;     % plot the entrenchment index though time
    LatVerInc_plot_entrenchmentGeom(iEC,talusoutC,WC,WFPC)

end

% save output of the model run for reference
savefile = ['incisionData_run',num2str(runno),'_',num2str(l),'.mat'] ;
save(savefile,'elevC','iEC','IeffC','I1C','I2C','dCC','SC','talusoutC',...
     'WC','WFPC') ;

if nruns>1
    % clear variables at the end of the run
    clear A Ai Aeq AWnorm Fx Hh Nplot S Si W WFP cx cz dt dx IXS lambda...
          max phi u v w x xC y yC
end

end


% F I N A L   M U L T I P L E    P L O T S   F O R   3 x 3 
% S I M U L A T I O N S (if there are 9 runs) _____________________________
if nruns == 9

    figure(2); clf;
    for i=1:9
        subplot(3,3,i)  % plot 
        data=load(['./incisionData_run',num2str(runno),'_',...
                   num2str(i),'.mat']);
            for j=1:5
                scatter(1:k,data.IeffC(j,1:k)./data.dCC(j,1:k),...
                        'filled','k')
                hold on
            end
        clear data
        ylim([0 1])
        xlim([0 k])
        if i == 4
            ylabel('vertical incision / lateral erosion','Fontsize',13)
        elseif i == 8
            xlabel('time','Fontsize',13)
        end    
    end

    set(gcf,'NextPlot','add');
    axes;
    h = title('Ratio of vertical to lateral erosion','fontsize',15,...
              'FontWeight','bold');
    set(gca,'Visible','off');
    set(h,'Visible','on');


    figure(3); clf;
    for i = 1:9 
        subplot(3,3,i)
        data=load(['./incisionData_run',num2str(runno),'_',...
                   num2str(i),'.mat']);
            for j=1:length(Qs)
               plot(1:k,data.elevC(j,1:k),'k')
               hold on
            end
        clear data
        if i == 4
            ylabel('elevation','Fontsize',13)
        elseif i == 8
            xlabel('time','Fontsize',13)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
        end
    end

    set(gcf,'NextPlot','add');
    axes;
    h = title('Elevation of channels through time','fontsize',15,...
              'FontWeight','bold');
    set(gca,'Visible','off');
    set(h,'Visible','on');

    figure(4); clf;
    for i = 1:9 
        subplot(3,3,i)
        data=load(['./incisionData_run',num2str(runno),'_',...
                   num2str(i),'.mat']);
            for j=1:size(data.iEC,1)
               plot(data.iEC(j,:))
               hold on
            end
            legend ('1','2','3','4','Location','southeast')
        clear data
        ylim([0.3 1])
        if i == 4
            ylabel('funnel index','Fontsize',13)
        elseif i == 8
            xlabel('time','Fontsize',13)
        end
    end

    set(gcf,'NextPlot','add');
    axes;
    h = title('Funnelity of the channel','fontsize',15,...
              'FontWeight','bold');
    set(gca,'Visible','off');
    set(h,'Visible','on');

end
