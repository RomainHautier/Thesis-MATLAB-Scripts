clear; 
clc;
close all;

speed = [10 20 40 50 60 80 100];
solver = 'F';
figcount = 1;
rotors = {'Main Rotor', 'Tail Rotor'};
What_to_Plot = {'MR', 'TR'};
          
%            Standalone Cases    2 ----- rotor ---- Cases
to_plot = {       'n'         ,           'y'            };

%% Extracting the circulation data from the solver for the standalone Rotor Cases.

if to_plot{1} == 'y'
    for s = 1:length(speed)
        filepath = 'C:\Users\hauti\heli-fwmrtr2\MyCasesNew-8turns\';
        [x,y,CIRC] = gammaplotf(speed(s),What_to_Plot,'F', filepath);
    
        for i = 1:length(x(1,1,:))
    
            figure(figcount)
            contourf(x(:,:,i),y(:,:,i),CIRC(:,:,i))
            colorbar;
            title(['Circulation distribution over the Standalone ' rotors{i} ' at ' num2str(speed(s)) ' kts'])
            figcount = figcount + 1;
    
        end
    
    end
end

%% Extracting the circulation data from the solver for the 2 Rotor Cases.

% Defining the 2 rotor cases for which to plot the circulation distribution
% over the rotor disk.
What_to_Plot = {'nb'};

if to_plot{2} == 'y'
    for s = 1:length(speed)
        filepath = 'C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\8turns';
        [x,y,CIRC] = gammaplotf(speed(s),What_to_Plot,'F', filepath);
    
        for i = 1:length(x(1,1,:))
    
            figure(figcount)
            contourf(x(:,:,i),y(:,:,i),CIRC(:,:,i))
            colorbar;
            title(['Circulation over the ' rotors{i} 'at ' num2str(speed(s)) '  kts - ' What_to_Plot], 'Interpreter', 'latex');
            figcount = figcount + 1;
    
        end
    
    end
end


