clear; 
clc;
close all;
save = 'y'; % Figure Saving Option
dp = 10;
fil_index = [1 10 19 28];
turns = 8;
nzt = turns*2*pi/(deg2rad(dp))+1; % Number of Collocation Points on a Filament.
rad = [8.18 1.676];
radr = rad(2)/rad(1);
radscale = [radr 1/radr];
offset = [0 0 0
    10.15138 -0.36401 0.25194];
sol_cases = {'nb'};
figcount = 1;
speed = [10:5:85];
% speed = [25 30 35 45 55 65 75 85];

coord_offset = {[0,0,0],[-10.15138,0.36401,-0.25194]};
rotors = {'Main Rotor', 'Tail Rotor'};
radius = [8.18, 1.676];
om = [27, 124.6];

cc = 1;
for s = 1:length(sol_cases)
    % Setting the directory in which the plots need to be saved.
    saveDir = ['C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\' num2str(turns) 'turns\figures\' sol_cases{s} '\'];
    for i = 1:length(speed)
        
        for j = 1:length(fil_index)
            figure(figcount)
            [pcx, pcy, pcz] = wakefilamentplot(2, ['C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\' num2str(turns) 'turns\d' num2str(speed(i)) '\' sol_cases{s} '\FWGEOM.DAT'],turns);
            plot3(pcx(fil_index(j),1:nzt,1),pcy(fil_index(j),1:nzt,1),pcz(fil_index(j),1:nzt,1),'Color',[0 0.4470 0.7410],'LineWidth',0.6)
            hold on;
            scatter3(pcx(fil_index(j),1:nzt,1),pcy(fil_index(j),1:nzt,1),pcz(fil_index(j),1:nzt,1),1,[0 0.4470 0.7410],'*')
            plot3(pcx(fil_index(j),1:nzt,2)*radr,pcy(fil_index(j),1:nzt,2)*radr,pcz(fil_index(j),1:nzt,2)*radr,'Color',[0.8500 0.3250 0.0980],'LineWidth',0.6)
            scatter3(pcx(fil_index(j),1:nzt,2)*radr,pcy(fil_index(j),1:nzt,2)*radr,pcz(fil_index(j),1:nzt,2)*radr,1,[0.8500 0.3250 0.0980],'*')
        end
        
        hold on;
        % set(gca,'YMinorTick','on')
        % set(gca,'XMinorTick','on')
        % ax                      = gca;     % gca = get current axes, and store this information in ax
        % ax.FontSize             = 14;      % set the property 'FontSize' to 14
        % ax.LineWidth            = 1.05;    % set the box around the figure to line width 1.05
        % ax.XAxis.Exponent       = 0;       % force the x-axis exponent
        % ax.YAxis.Exponent       = 0;       % force the y-axis exponent
        % ax.TickLabelInterpreter = 'latex'; %
        % % scatter3(10.15138,-0.36401,0.25194,36,'filled','r')
        % % plot3(circle_x(2,:), circle_y(2,:), circle_z(2,:), 'r-', 'LineWidth', 2);
        % xlabel('$\frac{x}{R}$','Interpreter','latex')
        % xlimit = ceil(max(max(pcx(:,:,2)))*radr);
        % xlim([-1.2 xlimit])
        % ylabel('$\frac{y}{R}$','Interpreter','latex')
        % zlabel('$\frac{z}{R}$','Interpreter','latex')
        % title(['Two Rotor Case Wake Geometries - ' num2str(speed(i)) ' Solution: ' sol_cases(s)],'Interpreter','latex')
        % legend('Main Rotor','', 'Tail Rotor', 'Interpreter', 'latex')
        % set(gca, 'XTick', 0:100:1000)


        %% Introducing the Disk Inflow Plots

        [~, VBZPSImat, inflownorm] = importVBZPSI(1,speed(i),[sol_cases{s} '\']);
        inflows_2rot_VBZPSI(:,i,s) = inflownorm; % SA stands for 'standalone'
        VBZPSImatout(:,:,1) = VBZPSImat(1:length(VBZPSImat(:,1))/2,:);
        VBZPSImatout(:,:,2) = VBZPSImat(length(VBZPSImat(:,1))/2+1:length(VBZPSImat(:,1)),:);
        inflowcoords = VBZPSImatout(:,1:3,:);
        vbzi = VBZPSImatout(:,6,:);
        scale(1,1,1:2) = [om(1)*rad(1), om(2)/rad(2)];
        scale_matrix(1,1:length(vbzi(1,:,1)),1) = ones*om(1)*rad(1);
        scale_matrix(1,1:length(vbzi(1,:,1)),2) = ones*om(2)*rad(2);
        vscaled = -vbzi./scale_matrix;
        medscale(1,1:2,i) = median(vscaled);
        avgV1(1,1:2,i) = mean(vscaled);
        medscale(1,3,i) = medscale(1,2,i)/medscale(1,1,i);
        avgV1(1,3,i) = avgV1(1,2,i)/avgV1(1,1,i);
        
        % Necessary initialisation for the position of the zlabel on the
        % plots.
        zZposold = -10;
        xZposold = 0;

        for j = 1:length(VBZPSImatout(1,1,:))
    
            x = inflowcoords(:,1,j)/rad(1);
            X = reshape(x,[32,36]);
            y = inflowcoords(:,2,j)/rad(1);
            Y = reshape(y,[32,36]);
            z = inflowcoords(:,3,j)/rad(1);
            Z = reshape(z,[32,36]);
            Z_surface = max(max(Z))*ones([32,37]);
            
            if j == 1
                % vbzi(:,:,j) = vbzi(:,:,j) * avgV1(1,3,i);
                vbzi(:,:,j) = vbzi(:,:,j) * medscale(1,3,i);
            end
            V = -reshape(vbzi(:,:,j),[32,36])./(om(j)*radius(j));
            
            % if vscaledr == V
            %     cc
            %     disp('yes')
            %     cc = cc + 1;
            % end

            avgV2(1,j,i) = mean(mean(V));
            X(:,37) = X(:,1);
            Y(:,37) = Y(:,1);
            Z(:,37) = Z(:,1);
            V(:,37) = V(:,1);
            

            % Defining the levels of dicretization needed for the inflow
            % plots.
            levels = linspace(min(min(V)),max(max(V)),60);
            V_discrete = interp1(levels, levels, V, 'nearest', 'extrap');
            surf(X,Y,Z,V,'EdgeColor','none','FaceAlpha',0.85)
            shading interp;
            % Set the colormap
            colormap(jet(length(levels)-1)); % Using jet colormap, change if desired

            % Adjust the color axis to match the levels
            caxis([min(levels) max(levels)]);

            % Add colorbar with the specified levels
            cb = colorbar;
            cb.Ticks = levels(1:10:length(levels));
            % cb.TickLabels = arrayfun(@num2str, levels, 'UniformOutput', false);
            contour3(X,Y,Z,V,levels(1:3:end),'k'); % Contour lines in black
            % axis equal;
            ax                      = gca;     % gca = get current axes, and store this information in ax
            ax.FontSize             = 14;      % set the property 'FontSize' to 14
            ax.LineWidth            = 1.1;    % set the box around the figure to line width 1.05
            ax.XAxis.Exponent       = 0;       % force the x-axis exponent
            ax.YAxis.Exponent       = 0;       % force the y-axis exponent
            ax.TickLabelInterpreter = 'latex'; %
            xl = xlabel('$\frac{x}{R_{MR}}$','Interpreter','latex','FontSize',25);
            xlimit = ceil(max(max(pcx(:,:,2)))*radr);
            xlimvec = [-1 xlimit];
            xlim(xlimvec);
            ylimvec = [-1 1.5];
            ylim(ylimvec)
            
            % Setting the zlims of the plot depending on max and min values
            % of pcz scaled wrt each rotor.
            for d = 1:length(pcz(1,1,:))
                zveclim(d,:) = [min(min(pcz(:,:,d))), max(max(pcz(:,:,d)))];

                if d == 2
                    zveclim(2,:) = zveclim(2,:).*radr;
                end
            end
            zveclim(1,1) = min(zveclim(:,1));
            zveclim(1,2) = max(zveclim(:,2));
            zveclim(2,:) = [];
            zlimround = [floor(zveclim(1)), ceil(zveclim(2))];
            for zi = 1:length(zlimround)
                if zi == 1 && zlimround(zi)+0.5 <= zveclim(zi)
                    zlimround(zi) = zlimround(zi)+0.5;
                elseif zi == 2 && zlimround(zi)-0.5 >= zveclim(zi)
                    zlimround(2) = zlimround(zi)-0.5;
                end
            end
            zlim(zlimround)


            yl = ylabel('$\frac{y}{R_{MR}}$','Interpreter','latex','FontSize',25);
            zl = zlabel('$\frac{z}{R_{MR}}$','Interpreter','latex','FontSize',25);
            
            % zZpos(j) = median(Z(:));
            % xZpos = min(min(min(pcz(:,:,:))));
            % yXpos = min(xlimvec)-0.3;
            % yYpos = min(min(min(pcy(:,:,:))));
            % 
            % if j == 1
            %     zYpos = max(Y(:));
            % end
            % 
            % 
            % if j == 2
            %     xXpos = median([-1 xlimit]);
            %     xYpos = 1.4;
            %     set(xl, 'Position', [xXpos, xYpos, xZpos]); % Adjust position
            %     set(yl, 'Position', [yXpos, 0, xZpos]); % Adjust position
            %     zXpos = max(xlimvec)+0.4;
            %     set(zl, 'Position', [zXpos, zYpos, min(zZpos)-1.3]); % Adjust position
            %     % set(zl, 'HorizontalAlignment', 'left'); % Align the text to the right
                set(zl, 'Rotation', 0); % Set rotation to 0 degrees (horizontal)
            % end
            % axis equal    
            title(['Two Rotor Case - ' num2str(speed(i)) ' Solution: ' sol_cases(s)],'Interpreter','latex','FontSize',16)
            legend('Main Rotor','', 'Tail Rotor', 'Interpreter', 'latex','Fontsize',13)
            set(gca, 'XTick', -10:0.5:10)
            set(gca, 'YTick', -10:0.5:10)
            set(gca, 'ZTick', -10:0.5:10)
            view(-170, 20); % Isometric view (default)
            
            % Formatting & Figure Saving
            % Set axes background color to white
            set(gca, 'Color', 'w');
            % Set axes background color to white
            set(gca, 'Color', 'w');
            % Remove white space around the figure
            % set(gca, 'Units', 'normalized', 'Position', [0 0 1 1]);
            % Define the desired figure size in pixels
            width = 840;  % Width of the figure
            height = 550; % Height of the figure
            
            % Adjust the figure position and size
            set(gcf, 'Position', [100, 20, width, height]); % [left, bottom, width, height]


            if save == 'y'
                
            end
        end
        figcount = figcount + 1;    
    end
end
    
        % % Saving the figure to the desired directory
        % if save == 'y'
        %     % fileName = fullfile(saveDir,['2R' num2str(speed(i)) sol_cases{s}]);
        %     fileName = [saveDir 'GEOM2R' num2str(speed(i)) sol_cases{s}];
        %     print(gcf,fileName,'-dpdf')
        % end
   