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
speed = [15:5:100];
speed = 100;
% speed = [10 15 25 35 40 55 65 70 75 85 90 95 100]
speed = [10 40 70 100];
coord_offset = {[0,0,0],[-10.15138,0.36401,-0.25194]};
rotors = {'Main Rotor', 'Tail Rotor'};
radius = [8.18, 1.676];
om = [27, 124.6];
y_disp = -[0 5 8 11 14];
z_disp = [0 0 0 0];
colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560],[0.6350 0.0780 0.1840]};
colors = {[0 0.4470 0.7410],[0 0.4470 0.7410],[0 0.4470 0.7410],[0 0.4470 0.7410]};

cc = 1;
for s = 1:length(sol_cases)
    % Setting the directory in which the plots need to be saved.
    saveDir = ['C:\Users\hauti\OneDrive - Imperial College London\Documents\London 23-24\Research Project\Result Figures\WakeGEOandInflow\SA\'];
    for i = 1:length(speed)

        [pcx, pcy, pcz] = wakefilamentplot(2, ['C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\' num2str(turns) 'turns\d' num2str(speed(i)) '\IWGEOM.DATA'],turns);
        % figure(figcount)
        for j = 1:length(fil_index)
            
            % plot3(pcx(fil_index(j),1:nzt,tile_index),pcy(fil_index(j),1:nzt,tile_index),pcz(fil_index(j),1:nzt,tile_index),'Color',colors{tile_index},'LineWidth',0.6)
            % hold on;
            % scatter3(pcx(fil_index(j),1:nzt,tile_index),pcy(fil_index(j),1:nzt,tile_index),pcz(fil_index(j),1:nzt,tile_index),1,colors{tile_index},'*')
            
            [pcx, pcy, pcz] = wakefilamentplot(2, ['C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\' num2str(turns) 'turns\d' num2str(speed(i)) '\IWGEOM.DATA'],turns);
            % plot3(pcx(fil_index(j),1:nzt,1),pcy(fil_index(j),1:nzt,1)+y_disp(i),pcz(fil_index(j),1:nzt,1),'Color',colors{i},'LineWidth',0.6)
            % hold on;
            % scatter3(pcx(fil_index(j),1:nzt,1),pcy(fil_index(j),1:nzt,1)+y_disp(i),pcz(fil_index(j),1:nzt,1),1,colors{i},'*')
            plot3(pcx(fil_index(j),1:nzt,2)*radr+coord_offset{2}(1)/rad(1),(pcy(fil_index(j),1:nzt,2)+y_disp(i))*radr,(pcz(fil_index(j),1:nzt,2)+z_disp(i))*radr,'Color',[0.8500 0.3250 0.0980],'LineWidth',0.6)
            hold on;
            scatter3(pcx(fil_index(j),1:nzt,2)*radr+coord_offset{2}(1)/rad(1),(pcy(fil_index(j),1:nzt,2)+y_disp(i))*radr,(pcz(fil_index(j),1:nzt,2)+z_disp(i))*radr,1,[0.8500 0.3250 0.0980],'*')
        end
        axis equal
        % hold on;
        
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

        [VBZPSImatout, VBZPSImat, inflownorm] = importVBZPSI(1,speed(i),[sol_cases{s} '\']);
        %inflows_2rot_VBZPSI(:,i,s) = inflownorm; % SA stands for 'standalone'
        % VBZPSImatout(:,:,1) = VBZPSImat(1:length(VBZPSImat(:,1))/2,:);
        % VBZPSImatout(:,:,2) = VBZPSImat(length(VBZPSImat(:,1))/2+1:length(VBZPSImat(:,1)),:);
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

        % for j = 1:length(VBZPSImatout(1,1,:))
    
            % x = inflowcoords(:,1,2)/rad(1);
            % X = reshape(x,[32,36]);
            % y = inflowcoords(:,2,2)/rad(1);
            % Y = reshape(y,[32,36]);
            % z = inflowcoords(:,3,2)/rad(1);
            % Z = reshape(z,[32,36]);
            % Z_surface = max(max(Z))*ones([32,37]);
            x = inflowcoords(:,1,2)/rad(1);
            X = reshape(x,[32,36]);
            y = inflowcoords(:,2,2)/rad(1);
            Y = reshape(y,[32,36]);
            z = inflowcoords(:,3,2)/rad(1);
            Z = reshape(z,[32,36]);
            Z_surface = max(max(Z))*ones([32,37]);
           
            
            % if j == 1
            %     % vbzi(:,:,j) = vbzi(:,:,j) * avgV1(1,3,i);
            %     vbzi(:,:,1) = vbzi(:,:,1) * medscale(1,3,i);
            % end
           

            V = -reshape(vbzi(:,:,2),[32,36])./(om(2)*radius(2));
            V = V/max(max(V))

            % figure(40)
            % X(:,37) = X(:,1);
            % Y(:,37) = Y(:,1);
            % Z(:,37) = Z(:,1);
            % V(:,37) = V(:,1);
            % levels = linspace(min(min(V)),max(max(V)),60);
            % contourf(X,Y,V,levels,'LineStyle','none')
            % % levels = linspace(min(min(V)),max(max(V)),60);
            % colormap(jet(length(levels)-1)); % Using jet colormap, change if desired
            % 
            % % Adjust the color axis to match the levels
            % clim([min(levels) max(levels)]);
            % colorbar
            % view(180,-90)
            % hold off;
            % figure(41)
            % surf(X,Y,Z,V)
            % hold on;
            % plot(X(:,10),Y(:,10),'g','LineWidth',1.6)
            % plot(X(:,19),Y(:,19),'b','LineWidth',1.6)
            % plot(X(:,28),Y(:,28),'k','LineWidth',1.6)
            % levels = linspace(min(min(V)),max(max(V)),60);
            % colormap(jet(length(levels)-1)); % Using jet colormap, change if desired
            % 
            % % Adjust the color axis to match the levels
            % clim([min(levels) max(levels)]);
            % colorbar
            % % cb.Limits = [-0.6 1]
            % shading interp
            % view(180,-90)
            % hold off;

            % vscaled2 = reshape(vscaled(:,:,j),[32 36]);
            % if vscaled2 == V
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
            figure(figcount)
            levels = linspace(-0.6,max(max(V)),60);
            V_discrete = interp1(levels, levels, V, 'nearest', 'extrap');
            % surf(X+coord_offset{2}(1)/rad(1),Y+y_disp(i)*radr,Z+z_disp(i)*radr,V,'EdgeColor','none','FaceAlpha',0.85)
            % surf(X,Y,Z,V,'EdgeColor','none','FaceAlpha',0.85)
            % shading interp;
            % Set the colormap
            % colormap(jet(length(levels)-1)); % Using jet colormap, change if desired

            % Adjust the color axis to match the levels
            % clim([min(levels) max(levels)]);

            % Add colorbar with the specified levels
            % cb = colorbar;
            % colorbar('Location','westoutside')
            % cb.Limits = [-0.6 1]
            % pos = get(cb, 'Position');

            % Modify the position (example: move it to the right)
            % pos(1) = pos(1) +0.2; % Move to the right
            % set(cb, 'Position', pos);
            % cb.Ticks = levels(1:10:length(levels));
            % Move the color bar outside the plot area
            % cb.Position = [0.9, 0.2, 0.03, 0.8]; % Adjust these values as needed
            % cb.Position = [pos(1), pos(2), pos(3), 0.5]; % Adjust the height value (0.5) as needed
            % cb.TickLabels = arrayfun(@num2str, levels, 'UniformOutput', false);
            % contour3(X+coord_offset{2}(1)/rad(1),Y,Z,V,levels(1:3:end),'k'); % Contour lines in black
            % axis equal;
            ax                      = gca;     % gca = get current axes, and store this information in ax
            ax.FontSize             = 26;      % set the property 'FontSize' to 14
            ax.LineWidth            = 1.1;    % set the box around the figure to line width 1.05
            ax.XAxis.Exponent       = 0;       % force the x-axis exponent
            ax.YAxis.Exponent       = 0;       % force the y-axis exponent
            ax.TickLabelInterpreter = 'latex'; %
            xl = xlabel('$\frac{x}{R_{MR}}$','Interpreter','latex','FontSize',46);
            xlimit = ceil(max(max(pcx(:,:,2)))*radr);
            zlimit = find( pcx(:,:,2) == max(max(pcx(:,:,2))) );
            % set(gca, 'YDir', 'reverse');
            % zlim([pcx(zlimit) zlimit]) 
            xlimvec = [-0.5 xlimit];
            % xlim([-0.25 2]);
            ylimvec = [-2.75 0.25];
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
            % zlim(zlimround)
           zlim([-0.4 0.4])

            yl = ylabel('$\frac{y}{R_{MR}}$','Interpreter','latex','FontSize',54);
            % labelPosition = get(yl,'Position');
            % labelPosition(2) = 0;  % Adjust the 0.5 as needed
            zl = zlabel('$\frac{z}{R_{MR}}$','Interpreter','latex','FontSize',54);
            % labelPosition = get(zl,'Position');
            % labelPosition(3) = labelPosition(3)+1  % Adjust the 0.5 as needed

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
            %     set(yl, 'Position', [-1.5, 0.4, -0.6]); % Adjust position
            % %     zXpos = max(xlimvec)+0.4;
            %     set(zl, 'Position', [xlimvec(end)+0.5, ylimvec(end), -0.1]); % Adjust position
            % %     % set(zl, 'HorizontalAlignment', 'left'); % Align the text to the right
                set(zl, 'Rotation', 0); % Set rotation to 0 degrees (horizontal)
            % % end
            % axis equal    
            % title(['Two Rotor Case - ' num2str(speed(i)) ' Solution: ' sol_cases(s)],'Interpreter','latex','FontSize',16)
            % legend('off')
            % legend(num2str(speed(1)),'',num2str(speed(2)), 'Interpreter', 'latex','Fontsize',13,'Location','northwest')
            set(gca, 'XTick', -10:0.5:20)
            set(gca, 'YTick', -10:0.5:10)
            set(gca, 'ZTick', -10:0.25:10)
            % set(gca, 'YDir', 'reverse');
            view(90, 0); % Isometric view (default)
            ytext = [-0.5 -1.3 -1.9 -2.5]
            text(0, ytext(i), 0.32, [num2str(speed(i)) ' kts'], 'FontSize', 36, 'Color', 'k','Interpreter','latex');
            
            % Formatting & Figure Saving
            % Set axes background color to white
            % set(gca, 'Color', 'w');
            % Remove white space around the figure
            % set(gca, 'Units', 'normalized', 'Position', [0 0 1 1]);
            % Define the desired figure size in pixels
            width = 800;  % Width of the figure
            height = 600; % Height of the figure
            
            % Adjust the figure position and size
            set(gcf, 'Position', [0, 0, width, height]); % [left, bottom, width, height]
        
            if save == 'y'
                
            end
        end
        % figcount = figcount + 1;    

        % Saving the figure to the desired directory
        % if save == 'y'
        %     fileName = [saveDir 'GEOMMR70100.pdf'];
        %     exportgraphics(gcf, fileName, 'ContentType', 'vector');
        % end
end


% end


% % Create a new figure with a tiled layout
% figure;
% t = tiledlayout(2, 2);
% 
% % Capture each existing figure and display in the tiled layout
% 
% % Tile 1
% nexttile;
% img1 = getframe(figure(1)); % Capture the content of figure1
% imshow(img1.cdata); % Display the captured image in the tile
% 
% % Tile 2
% nexttile;
% img2 = getframe(figure(2)); % Capture the content of figure2
% imshow(img2.cdata); % Display the captured image in the tile
% 
% % Tile 3
% nexttile;
% img3 = getframe(figure(3)); % Capture the content of figure3
% imshow(img3.cdata); % Display the captured image in the tile
% 
% % Tile 4
% nexttile;
% img4 = getframe(figure(4)); % Capture the content of figure4
% imshow(img4.cdata); % Display the captured image in the tile



    
        
   