clear all;
clc;
close all;
dp = 10; % Azimuthal resolution.
figcount = 1;
fil_index = [1 10 19 28];
% Setting up the plot options.
%             helMR_FWTR     MR_TR_FW        nn         nb        
plotoptions = ['y';             'n';        'n'   ;    'n'];
save = 'y';

%% Plotting the IWGEOM.data from standalone TR FW & HeliUM runs.
speed = [50:5:85];
rad = [8.18 1.676];
radr = rad(2)/rad(1);
offset = [0 0 0
    10.15138 -0.36401 0.25194];
theta = linspace(0, 2*pi, 100); % Parameter to define the circle

for rr = 1:length(rad)
    circle_x(rr,:) = offset(rr,1) + rad(rr) * cos(theta); % x-coordinates of the circle
    circle_y(rr,:) = offset(rr,2) + rad(rr) * sin(theta); % y-coordinates of the circle
    circle_z(rr,:) = offset(rr,3) + zeros(size(theta)); % z-coordinates of the circle (0 for xy-plane)
end

%% Plotting the FWGEOM.dat for each speed for the 2 rotor cases
sol_cases = {'nb'};
turns = 8; % Number of Free Wake Turns
nzt = turns*2*pi/(deg2rad(dp))+1; % Number of Collocation Points on a Filament.

% This loop plots the wake geometries of both MR & TR when they are ran in
% the 2 rotor dissimilar geometry case.
% 1. Iterates of the 'solution methodology' used in the freewake algorithm.
%   This concernes the 'trm' and the 'initial' options in the user.input
%   file.
% 2. Iterates over the range of speeds for which the wake geometries are to
%   be plotted.
% 3. Extracts the wake geometry data using the wakefilamentplot.m function.
% 4. Plots the wake geometries with a number of filaments defined by the
% 'fil_index' array defined above.

for s = 1:length(sol_cases)
    % Setting the directory in which the plots need to be saved.
    saveDir = ['C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\' num2str(turns) 'turns\figures\' sol_cases{s} '\'];
    
    for i = 1:length(speed)
        figure(figcount)
        for j = 1:length(fil_index)
            [pcx, pcy, pcz] = wakefilamentplot(2, ['C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\' num2str(turns) 'turns\d' num2str(speed(i)) '\' sol_cases{s} '\FWGEOM.DAT'],turns);
            plot3(pcx(fil_index(j),1:nzt,1),pcy(fil_index(j),1:nzt,1),pcz(fil_index(j),1:nzt,1),'Color',[0 0.4470 0.7410],'LineWidth',0.6)
            hold on;
            scatter3(pcx(fil_index(j),1:nzt,1),pcy(fil_index(j),1:nzt,1),pcz(fil_index(j),1:nzt,1),1,[0 0.4470 0.7410],'*')
            plot3(pcx(fil_index(j),1:nzt,2)*radr,pcy(fil_index(j),1:nzt,2)*radr,pcz(fil_index(j),1:nzt,2)*radr,'Color',[0.8500 0.3250 0.0980],'LineWidth',0.6)
            scatter3(pcx(fil_index(j),1:nzt,2)*radr,pcy(fil_index(j),1:nzt,2)*radr,pcz(fil_index(j),1:nzt,2)*radr,1,[0.8500 0.3250 0.0980],'*')
        end
        blade_endpoints_MR = [pcx(fil_index,1,1)*rad(1),pcy(fil_index,1,1)*rad(1),pcz(fil_index,1,1)*rad(1)]/rad(1);
        blade_endpoints_TR = [pcx(fil_index,1,2)*rad(2),pcy(fil_index,1,2)*rad(2),pcz(fil_index,1,2)*rad(2)]/rad(1);
        
        for h = 1:2
            plot3(blade_endpoints_MR([h h+2],1),blade_endpoints_MR([h h+2],2),blade_endpoints_MR([h h+2],3),'k-','LineWidth',2)
            plot3(blade_endpoints_TR([h h+2],1),blade_endpoints_TR([h h+2],2),blade_endpoints_TR([h h+2],3),'k-','LineWidth',2)
        end
        axis equal
        % plot3(circle_x(1,:), circle_y(1,:), circle_z(1,:), 'k-', 'LineWidth', 2);
        hold on;
        set(gca,'YMinorTick','on')
        set(gca,'XMinorTick','on')
        ax                      = gca;     % gca = get current axes, and store this information in ax
        ax.FontSize             = 14;      % set the property 'FontSize' to 14
        ax.LineWidth            = 1.05;    % set the box around the figure to line width 1.05
        ax.XAxis.Exponent       = 0;       % force the x-axis exponent
        ax.YAxis.Exponent       = 0;       % force the y-axis exponent
        ax.TickLabelInterpreter = 'latex'; %
        % scatter3(10.15138,-0.36401,0.25194,36,'filled','r')
        % plot3(circle_x(2,:), circle_y(2,:), circle_z(2,:), 'r-', 'LineWidth', 2);
        xlabel('$\frac{x}{R}$','Interpreter','latex')
        xlimit = ceil(max(max(pcx(:,:,2)))*radr);
        xlim([-1.2 xlimit])
        ylabel('$\frac{y}{R}$','Interpreter','latex')
        zlabel('$\frac{z}{R}$','Interpreter','latex')
        title(['Two Rotor Case Wake Geometries - ' num2str(speed(i)) ' Solution: ' sol_cases(s)],'Interpreter','latex')
        legend('Main Rotor','', 'Tail Rotor', 'Interpreter', 'latex')
        set(gca, 'XTick', 0:100:1000)
        
        figcount = figcount + 1;

        % Saving the figure to the desired directory
        if save == 'y'
            fileName = [saveDir 'GEOM2R' num2str(speed(i)) sol_cases{s}];
            print(gcf,fileName,'-dpdf')
        end

    end
end

if save == 'y'
    disp('The figures plotted here have been saved to the desired directory');
end

if plotoptions(2) == 'y'
    for i = 1:length(speed)
        figure(figcount)
        for j = 1:length(fil_index)
            [pcx, pcy, pcz] = wakefilamentplot(2, ['C:\Users\hauti\heli-fwmrtr2\MyCases2rot\d' num2str(speed(i)) '\nn-nzt577\FWGEOM.DAT'],16);
            plot3(pcx(fil_index(j),1:577,1)*rad(1),pcy(fil_index(j),1:577,1)*rad(1),pcz(fil_index(j),1:577,1)*rad(1),'Color',[0 0.4470 0.7410],'LineWidth',1.5)
            hold on
            plot3(pcx(fil_index(j),1:577,2)*rad(2),pcy(fil_index(j),1:577,2)*rad(2),pcz(fil_index(j),1:577,2)*rad(2),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5)
        end
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title(['Freewake - Dissimilar Rotor Case ' num2str(speed(i)) ' knots Starting Wake Solution - 16 Wake Turns - nb'])
        legend('Main Rotor', 'Tail Rotor')
        
        figcount = figcount + 1;
    end
end


% trial_speed = [10 20 40 50 60 80 100];
% fil_index = [1 10 19 28];
% for i = 1:length(trial_speed)
%         figure(figcount)
%         [pcx, pcy, pcz] = wakefilamentplot(2,['C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\8turns\d' num2str(trial_speed(i)) '\fb\FWGEOM.DAT'],8);
%     for j = 1:length(fil_index)
%         plot3(pcx(fil_index(j),1:289,1)*rad(1),pcy(fil_index(j),1:289,1)*rad(1),pcz(fil_index(j),1:289,1)*rad(1),'Color',[0 0.4470 0.7410],'LineWidth',0.6)
%         hold on
%         plot3(pcx(fil_index(j),1:289,2)*rad(2),pcy(fil_index(j),1:289,2)*rad(2),pcz(fil_index(j),1:289,2)*rad(2),'Color',[0.8500 0.3250 0.0980],'LineWidth',0.6)
%     end 
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     % xlim([-rad(1) 2.5*rad(1)])
%     % zlim([-rad(1)/2 rad(1)/2])
%     %         view([0 0 1])
%     % % set(gca, 'XDir', 'reverse'); % Reverse the X-axis direction
%     legend('Main Rotor','Tail Rotor')
%     title(['Freewake - Dissimilar Rotor Case ' num2str(range_of_speed(i)) ' knots Starting Wake Solution - ' num2str(turns) ' Wake Turns - fb'])
%     figcount = figcount+1;
% end


% trial_speed = [10 20 50 100];
% % fil_index = [1 10]
% for i = 1:length(trial_speed)
%         figure(figcount)
%         [pcx, pcy, pcz] = wakefilamentplot(2,['C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\8turns\d' num2str(trial_speed(i)) '\fb\FWGEOM.DAT'],8);
%     for j = 1:length(fil_index)
%         plot3(pcx(fil_index(j),1:289,1)*rad(1),pcy(fil_index(j),1:289,1)*rad(1),pcz(fil_index(j),1:289,1)*rad(1),'b','LineWidth',0.6)
%         hold on
%         plot3(pcx(fil_index(j),1:289,2)*rad(2),pcy(fil_index(j),1:289,2)*rad(2),pcz(fil_index(j),1:289,2)*rad(2),'r','LineWidth',0.6)
%     end 
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     % xlim([-rad(1) 2.5*rad(1)])
%     % zlim([-rad(1)/2 rad(1)/2])
%     %         view([0 0 1])
%     % % set(gca, 'XDir', 'reverse'); % Reverse the X-axis direction
%     legend('Main Rotor','Tail Rotor')
%     title([num2str(trial_speed(i)) 'kts'])
%     figcount = figcount+1;
% end