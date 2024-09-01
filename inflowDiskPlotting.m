clear;
clc;
close all;

%% Output matrix of the VBZPSI
speed = [10 25 40];
speed = [10 25 40];
figcount = 1;

coord_offset = {[0,0,0],[-10.15138,0.36401,-0.25194]};
rotors = {'Main Rotor', 'Tail Rotor'};
radius = [8.18, 1.676];
om = [27, 124.6];
save = 'y';

plotoptions = {'n' , 'y'};

myMap = containers.Map({'10','25','40', '55', '70', '100'}, {[], [], [
    1.38137, -0.00972141, -0.0610419;
    1.33351, -0.0218192, -0.0230485;
    1.28976, -0.0285641, -0.00122973;
    1.19851, -0.0411523,  0.0318039;
    1.1399,  -0.0495153,  0.0511694;
    1.08828, -0.0611731,  0.0845571
], [
    1.3909, -0.0195566, -0.0337604;
    1.35308, -0.0259656, -0.0130284;
    1.28084, -0.0354985,  0.0185869;
    1.14575, -0.0490369,  0.0490718;
    1.20992, -0.0415083,  0.0329377;
    1.08719, -0.0575809,  0.0685203
], [
    1.36424, -0.035999, 0.0172404;
    1.27856, -0.0256666, 0.0332948;
    1.2207, -0.020672, 0.0436792;
    1.16363, -0.0160247, 0.0547439;
    1.10781, -0.0108176, 0.0660352
],[ 1.39812, -0.00788588, 0.00702376;
    1.36554, -0.0336386, 0.00357248;
    1.24258, -0.0338401, 0.0124773;
    1.31873, -0.00988401, 0.0102793;
    1.16376, -0.0376186, 0.0173825;
    1.07787, -0.00648282, 0.0196033
]});



if plotoptions{1} == 'y'
    %% Plotting the Inflow Distribution over the Standalone Rotor Disks
    tiles1 = tiledlayout(3,1);
    for i = 1:length(speed)
        [VBZPSImatout, VBZPSImat, inflownorm] = importVBZPSI(1,speed(i),'');
        inflows_SA_VBZPSI(:,i) = inflownorm; % SA stands for 'standalone'
        inflowcoords = VBZPSImatout(:,1:3,:);
        vbzi = VBZPSImatout(:,6,:);
        
        for j = 1:length(VBZPSImatout(1,1,:))

            x = inflowcoords(:,1,j)+coord_offset{j}(1);
            y = inflowcoords(:,2,j)+coord_offset{j}(2);
            z = inflowcoords(:,3,j)+coord_offset{j}(3);
    
            r = sqrt(x.^2+y.^2+z.^2);
            rad = reshape(r,[32, 36])/radius(j);
            psi_or = deg2rad(0:10:360);
            psi = zeros(32,36);
            for index = 1:length(psi(1,:))
                psi(:,index) = psi_or(index);
            end
            
            % Convert polar coordinates to Cartesian coordinates
            X = rad .* cos(psi);
            % X = reshape(x,[32 36]);
            X(:,37) = X(:,1);
            Y = rad .* sin(psi);
            % Y = reshape(y,[32 36]);
            Y(:,37) = Y(:,1);
            V = -reshape(vbzi(:,:,j),[32,36])./(om(j)*radius(j));
            V(:,37) = V(:,1);
            V = V/max(max(V));
            % Plot the data
            % figure(1);
            % scatter(X(:), Y(:), 'filled');
            axis equal;  % Ensure the plot is circular
            % xlabel('X');
            % ylabel('Y');
            % title('Disk Plot of Radial and Azimuthal Points');
            % grid on;
            % hold off;
            figure(figcount)
            h=polar([0 2*pi],[0,1]);
            set(gca, 'FontSize', 14);
            view([180 -90])
            delete(h)
            hold on
            levels = linspace(-0.6,1,60);
            contourf(X,Y,V,levels,'LineStyle','none')
            Z = zeros(size(X));
            % surf(X,Y,Z,V,'EdgeColor','none','FaceAlpha',1)
            % shading interp;
            % Set the colormap
            colormap(jet(length(levels)-1)); % Using jet colormap, change if desired

            % Adjust the color axis to match the levels
            clim([min(levels) max(levels)]);
            cb = colorbar
            cb.Limits = [-0.6 1]
            hold on; scatter(0,0)
            figcount = figcount + 1;
            title([rotors{j} ' Disk Inflow Distribution - ' num2str(speed(i)) ' kts']);
            
            
            if j == 1 || j == 3 || j ==5
                figure(7)
                nexttile
               
                
                h=polar([0 2*pi],[0,1]);
                textHandles = findall(gcf,'Type','text'); % Find all text objects
                set(textHandles, 'FontSize', 18);
                axis equal
                view([180 90])
                % Get the current axes
                
                % Set the overall font size for the axes (including tick labels)
                % set(ax, 'FontSize', 10);
                
                % Find all text objects within the axes and adjust their font size
                % texts = findall(ax, 'Type', 'text');
                % set(texts, 'FontSize', 10); % Adjust the font size as needed
                delete(h)
                hold on
                contourf(X,Y,V,levels,'LineStyle','none')
                colormap(jet(length(levels)-1)); % Using jet colormap, change if desired
                % title([num2str(speed(i)) ' kts'],'Interpreter','latex','FontSize',24);
                % Adjust the color axis to match the levels
                clim([min(levels) max(levels)]);
                % cb = colorbar
                % Create a figure
                fig = figure(7);

                % Set figure size (Width x Height in pixels)
                % fig.Position = [0, -200, 300, 950]; % [left, bottom, width, height]
            else
                figure(8)
                nexttile
                
                h=polar([0 2*pi],[0,1]);
                axis equal
                ax = gca;
                % thetaLabels = ax.ThetaTickLabel;
                textHandles = findall(gcf,'Type','text'); % Find all text objects
                set(textHandles, 'FontSize', 18); % Set font size of all text objects
                % Loop through each text object and adjust its position slightly outward
                % for th = 1:length(textHandles)
                % pos = get(textHandles(th), 'Position');
                % % Increase the radius component to move the labels outward
                % pos(1) = pos(1) * 1.07; % Adjust the factor to move the labels outward more or less
                % set(textHandles(th), 'Position', pos);
                % end
                % Define angle positions (in radians) for labels
                view([180 90])
                % Get the current axes
                ax = gca;
                
                % Set the overall font size for the axes (including tick labels)
                % set(ax, 'FontSize', 10);
                
                % Find all text objects within the axes and adjust their font size
                % texts = findall(ax, 'Type', 'text');
                % set(texts, 'FontSize', 10); % Adjust the font size as needed
                delete(h)
                hold on
                contourf(X,Y,V,levels,'LineStyle','none')
                colormap(jet(length(levels)-1)); % Using jet colormap, change if desired
                % title([num2str(speed(i)) ' kts'],'Interpreter','latex','FontSize',24);
                % Adjust the color axis to match the levels
                clim([min(levels) max(levels)]);
                hold on;
                % cb = colorbar
                % Create a figure
                fig = figure(8);

                % Set figure size (Width x Height in pixels)
                % fig.Position = [0, -200, 300, 950]; % [left, bottom, width, height]

            end

        end
    end

    figure(8)
    ztext = [1 2 3];
    % for sp = 1:3
    %     text(0, 0, ztext(sp), [num2str(speed(sp)) ' kts'], 'FontSize', 24, 'Color', 'k','Interpreter','latex');
    %     hold on;
    % end
end

if plotoptions{2} == 'y'
    %% Plotting the Inflow Distribution over the Interacting Rotors Disks
sol_cases = {'nb'};
for s = 1:length(sol_cases)
    tiles = tiledlayout(1,3);
    figure(figcount)
    for i = 1:length(speed)
        [~, VBZPSImat, inflownorm] = importVBZPSI(2,speed(i),[sol_cases{s} '\']);
        inflows_2rot_VBZPSI(:,i,s) = inflownorm; % SA stands for 'standalone'
        VBZPSImatout(:,:,1) = VBZPSImat(1:length(VBZPSImat(:,1))/2,:);
        VBZPSImatout(:,:,2) = VBZPSImat(length(VBZPSImat(:,1))/2+1:length(VBZPSImat(:,1)),:);
        inflowcoords = VBZPSImatout(:,1:3,:);
        vbzi = VBZPSImatout(:,6,:);

        
        for j = 1:length(VBZPSImatout(1,1,:))
    
            x = inflowcoords(:,1,j)+coord_offset{j}(1);
            y = inflowcoords(:,2,j)+coord_offset{j}(2);
            z = inflowcoords(:,3,j)+coord_offset{j}(3);
            % x = inflowcoords(:,1,j);
            % y = inflowcoords(:,2,j);
            % z = inflowcoords(:,3,j);


            % Scattering the points of main rotor interactions with the
            % tail rotor
            interaction_coord_array = myMap(num2str(speed(i)));
    
            r = sqrt(x.^2+y.^2+z.^2);
            rad = reshape(r,[32, 36])/radius(j);
            psi_or = deg2rad(0:10:360);
            psi = zeros(32,36);
            for index = 1:length(psi(1,:))
                psi(:,index) = psi_or(index);
            end
            
            % Convert polar coordinates to Cartesian coordinates
            X = rad .* cos(psi);
            X(:,37) = X(:,1);
            Y = rad .* sin(psi);
            Y(:,37) = Y(:,1);
            V = -reshape(vbzi(:,:,j),[32,36])./(om(j)*radius(j));
            V(:,37) = V(:,1);
            V = V/max(max(V));

            % x = inflowcoords(:,1,j)/rad(1);
            % X = reshape(x,[32,36]);
            % X(:,37) = X(:,1);
            % y = inflowcoords(:,2,j)/rad(1);
            % Y = reshape(y,[32,36]);
            % Y(:,37) = Y(:,1);
            z = inflowcoords(:,3,j)/rad(1);
            Z = reshape(z,[32,36]);
            Z_surface = max(max(Z))*ones([32,37]);
            
            % Plot the data
            % figure(1);
            % scatter(X(:), Y(:), 'filled');
            % axis equal;  % Ensure the plot is circular
            % xlabel('X');
            % ylabel('Y');
            % title('Disk Plot of Radial and Azimuthal Points');
            % grid on;
            % hold off;
            % figure(figcount)
            if j == 2
                nexttile
                
                h=polar([0 2*pi],[0,1]);
                textHandles = findall(gcf,'Type','text'); % Find all text objects
                set(textHandles, 'FontSize', 24);
                view([180 90])
                delete(h)
                hold on
                levels = linspace(-1,1,60);
                contourf(X,Y,V,levels,'LineStyle','none')
                hold on;
                axis equal
                % if isempty(interaction_coord_array) == 0
                %     for sc = 1:length(interaction_coord_array(:,1))
                %         scatter(xscat(sc),yscat(sc),36,[0 0 0],'filled','+');
                %         hold on;
                %     end
                % end
                % Set the colormap to match the number of levels
                colormap(jet(length(levels) - 1))
                
                % Add a colorbar and set its limits
                cb = colorbar;
                set(cb, 'FontName', 'Arial', 'FontSize', 22);
                caxis([-0.6 1]) % This sets the limits for the colormap and colorbar
                cb.Limits = [-0.6 1]; % Explicitly setting the colorbar limits
                hold on; scatter(0,0)
                if isempty(interaction_coord_array) == 0 && j == 2
                    for sc = 1:length(interaction_coord_array(:,1))
                        % xscat(sc) = interaction_coord_array(sc,1)*radius(1)+coord_offset{j}(1);
                        % xscat(sc)
                        % yscat(sc) = interaction_coord_array(sc,2)*radius(1)+coord_offset{j}(2);
                        % zscat(sc) = interaction_coord_array(sc,3)*radius(1)+coord_offset{j}(3);
                        % zscat(sc)
                         xscat = (interaction_coord_array(sc,1)*radius(1)+coord_offset{j}(1))/radius(2);
                        xscat
                        yscat = (interaction_coord_array(sc,2)*radius(1)+coord_offset{j}(2))/radius(2);
                        zscat = (interaction_coord_array(sc,3)*radius(1)+coord_offset{j}(3)/radius(2));
                        zscat
                        scatter(xscat,-zscat,200,[0 0 0],'*','LineWidth',2.5);
                        scatter(0,0,1,[0 0 0],'*')
                        legend('','','','','','Main Rotor Vortex Interaction Points','Interpreter','Latex','FontSize',28)
                        hold on;
                    end
                end

               

                
                % % title([rotors{j} ' Disk Inflow Distribution - ' num2str(speed(i)) ' kts - ' sol_cases{s}]);
                % title([ num2str(speed(i)) ' kts'],'Interpreter','latex');
            
            end

            % saveDir = ['C:\Users\hauti\OneDrive - Imperial College London\Documents\London 23-24\Research Project\Result Figures\WakeGEOandInflow\2rot\'];
            % if save == 'y'
            %     rotor = {'MR','TR'};
            %     fileName = [saveDir rotor{j} 'inflow' num2str(speed(i)) '.pdf'];
            %     exportgraphics(gcf, fileName, 'ContentType', 'vector');
            % end
            % figcount = figcount + 1;
            % title([rotors{j} ' Disk Inflow Distribution - ' num2str(speed(i)) ' kts - ' sol_cases{s}]);
            
            




        end
            % Define the desired figure size in pixels
            width = 2000;  % Width of the figure
            height = 800; % Height of the figure
            
            % Adjust the figure position and size
            set(gcf, 'Position', [0, 0, width, height]); % [left, bottom, width, height]
            text(0.27, -1.3, [num2str(speed(i)) ' kts'], 'FontSize', 40, 'Color', 'k','Interpreter','latex');
            
            % Optionally, adjust the marker size in the legend (not usually necessary if scatter size is set)
            % legendChildren = findobj(leg, 'type', 'patch');
            % set(legendChildren, 'Markersize', 1000); % Adjust if needed
    end
end
end


