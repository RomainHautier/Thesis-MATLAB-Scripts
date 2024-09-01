clear;
clc;
close all;

%% Output matrix of the VBZPSI
speed = [10 70 100];
figcount = 1;

coord_offset = {[0,0,0],[-10.15138,0.36401,-0.25194]};
rotors = {'Main Rotor', 'Tail Rotor'};
radius = [8.18, 1.676];
om = [27, 124.6];

plotoptions = {'y' , 'n'};

if plotoptions{1} == 'y'
    %% Plotting the Inflow Distribution over the Standalone Rotor Disks
    tiles1 = tiledlayout(3,1);
    for i = 1:length(speed)
        [VBZPSImatout, VBZPSImat, inflownorm] = importVBZPSI(1,speed(i),'');
        inflows_SA_VBZPSI(:,i) = inflownorm; % SA stands for 'standalone'
        inflowcoords = VBZPSImatout(:,1:3,:);
        vbzi = VBZPSImatout(:,6,:);
        [CLr1, CLr2] = importCL(2,speed(i));
        CLr2 = transpose(CLr2(:,:,1));
        CLr1 = transpose(CLr1(:,:,1));
        CLr1(:,37) =  CLr1(:,1);
        CLr2(:,37) =  CLr2(:,1);
        CLs(:,:,1) = CLr1;
        CLs(:,:,2) = CLr2;


        for j = 1:length(VBZPSImatout(1,1,:))

            % x = inflowcoords(:,1,j)+coord_offset{j}(1);
            % y = inflowcoords(:,2,j)+coord_offset{j}(2);
            % z = inflowcoords(:,3,j)+coord_offset{j}(3);
            % 
            % r = sqrt(x.^2+y.^2+z.^2);
            % rad = reshape(r,[32, 36])/radius(j);
            % psi_or = deg2rad(0:10:360);
            % psi = zeros(32,36);
            % for index = 1:length(psi(1,:))
            %     psi(:,index) = psi_or(index);
            % end
            % 
            % % Convert polar coordinates to Cartesian coordinates
            % X = rad .* cos(psi);
            % % X = reshape(x,[32 36]);
            % X(:,37) = X(:,1);
            % Y = rad .* sin(psi);
            % % Y = reshape(y,[32 36]);
            % Y(:,37) = Y(:,1);
            % V = -reshape(vbzi(:,:,j),[32,36])./(om(j)*radius(j));
            % V(:,37) = V(:,1);
            % V = V/max(max(V));
            % % Plot the data
            % % figure(1);
            % % scatter(X(:), Y(:), 'filled');
            % axis equal;  % Ensure the plot is circular
            % % xlabel('X');
            % % ylabel('Y');
            % % title('Disk Plot of Radial and Azimuthal Points');
            % % grid on;
            % % hold off;
            % Define the number of segments
            nAzimuth = 36;  % Number of azimuthal divisions
            nRadius = 32;   % Number of radial divisions
            
            % Create a vector for the azimuth angles
            theta = linspace(0, 2*pi, nAzimuth);
            
            % Create a vector for the radial positions
            r = linspace(0.069, 1, nRadius);
            
            % Create a grid using meshgrid
            [Theta, R] = meshgrid(theta, r);
            
            % Convert polar coordinates (R, Theta) to Cartesian coordinates (X, Y)
            X = R .* cos(Theta);
            X(:,37) = X(:,1);
            Y = R .* sin(Theta);
            Y(:,37) = Y(:,1);



            figure(figcount)
            h=polar([0 2*pi],[0,1]);
            view([180 90])
            delete(h)
            hold on
            levels = linspace(-0.6,1,60);
            % contourf(X,Y,V,levels,'LineStyle','none')
            contourf(X,Y,CLs(:,:,j),levels,'LineStyle','none')
            Z = zeros(size(X));
            % surf(X,Y,Z,V,'EdgeColor','none','FaceAlpha',1)
            % shading interp;
            % Set the colormap
            % colormap(jet(length(levels)-1)); % Using jet colormap, change if desired

            % Adjust the color axis to match the levels
            % clim([min(levels) max(levels)]);
            cb = colorbar
            % cb.Limits = [-0.6 1]
            hold on; scatter(0,0)
            figcount = figcount + 1;
            title([rotors{j} ' Disk Inflow Distribution - ' num2str(speed(i)) ' kts']);
            
            
            if j == 1 || j == 3 || j ==5
                figure(7)
                nexttile
                
                h=polar([0 2*pi],[0,1]);
                ax = gca;
                view([180 90])
                % Get the current axes
                ax = gca;
                
                % Set the overall font size for the axes (including tick labels)
                set(ax, 'FontSize', 10);
                
                % Find all text objects within the axes and adjust their font size
                texts = findall(ax, 'Type', 'text');
                set(texts, 'FontSize', 10); % Adjust the font size as needed
                delete(h)
                hold on
                % CLs(:,:,j) = CLs(:,:,j)/max(max(CLs(:,:,j)));
                % levels = linspace(min(min(CLs(:,:,j))),max(max(CLs(:,:,j))),60);
                contourf(X,Y,CLs(:,:,1),'LineStyle','none')
                % colormap(jet(length(levels)-1)); % Using jet colormap, change if desired
                title([num2str(speed(i)) ' kts'],'Interpreter','latex','FontSize',24);
                % Adjust the color axis to match the levels
                % clim([min(levels) max(levels)]);
                % cb = colorbar
                % Create a figure
                fig = figure(7);

                % Set figure size (Width x Height in pixels)
                fig.Position = [0, -200, 300, 950]; % [left, bottom, width, height]
            else
                figure(8)
                nexttile
                
                h=polar([0 2*pi],[0,1]);
                ax = gca;
                view([180 90])
                % Get the current axes
                ax = gca;
                
                % Set the overall font size for the axes (including tick labels)
                set(ax, 'FontSize', 10);
                
                % Find all text objects within the axes and adjust their font size
                texts = findall(ax, 'Type', 'text');
                set(texts, 'FontSize', 10); % Adjust the font size as needed
                delete(h)
                hold on
                CLs(:,:,2) = CLs(:,:,2)/max(max(CLs(:,:,2)));
                levels = linspace(min(min(CLs(:,:,2))),max(max(CLs(:,:,2))),10);
                contourf(X,Y,CLs(:,:,2),levels,'LineStyle','none')
                % contourf(X,Y,V,levels,'LineStyle','none')
                % colormap(jet(length(levels)-1)); % Using jet colormap, change if desired
                title([num2str(speed(i)) ' kts'],'Interpreter','latex','FontSize',24);
                % Adjust the color axis to match the levels
                clim([min(levels) max(levels)]);
                cb = colorbar
                % Create a figure
                fig = figure(8);

                % Set figure size (Width x Height in pixels)
                fig.Position = [0, -200, 300, 950]; % [left, bottom, width, height]

            end

        end
    end
end

if plotoptions{2} == 'y'
    %% Plotting the Inflow Distribution over the Interacting Rotors Disks
sol_cases = {'nb'};
for s = 1:length(sol_cases)
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
            % Plot the data
            % figure(1);
            % scatter(X(:), Y(:), 'filled');
            % axis equal;  % Ensure the plot is circular
            % xlabel('X');
            % ylabel('Y');
            % title('Disk Plot of Radial and Azimuthal Points');
            % grid on;
            % hold off;
            figure(figcount)
            h=polar([0 2*pi],[0,1]);
            delete(h)
            hold on
            % levels = linspace(1,-0.6,20);
            contour(X,Y,V,'LineStyle','none')
            
            % Set the colormap to match the number of levels
            % colormap(jet(length(levels) - 1))

            % Add a colorbar and set its limits
            cb = colorbar;
            caxis([-0.6 1]) % This sets the limits for the colormap and colorbar
            cb.Limits = [-0.6 1]; % Explicitly setting the colorbar limits
            hold on; scatter(0,0)
            figcount = figcount + 1;
            title([rotors{j} ' Disk Inflow Distribution - ' num2str(speed(i)) ' kts - ' sol_cases{s}]);
            
            




        end
    end
end
end

figure(figcount)

tiles = tiledlayout(3,1)
% axis equal

nexttile
copyobj(allchild(get(5,'CurrentAxes')), gca); % Copy the first figure's content
title('$10$ kts','Interpreter','latex');

nexttile
copyobj(allchild(get(3,'CurrentAxes')), gca); % Copy the first figure's content
title('$25$ kts','Interpreter','latex');

nexttile
copyobj(allchild(get(1,'CurrentAxes')), gca); % Copy the first figure's content
title('$40$ kts','Interpreter','latex');

figcount = figcount + 1;


% figure(figcount)
% 
% tiles = tiledlayout(3,1)
% % axis equal
% 
% nexttile
% copyobj(allchild(get(6,'CurrentAxes')), gca); % Copy the first figure's content
% title('$10$ kts','Interpreter','latex');
% 
% nexttile
% copyobj(allchild(get(4,'CurrentAxes')), gca); % Copy the first figure's content
% title('$25$ kts','Interpreter','latex');
% 
% nexttile
% copyobj(allchild(get(2,'CurrentAxes')), gca); % Copy the first figure's content
% title('$40$ kts','Interpreter','latex');
% 
% figcount = figcount + 1;