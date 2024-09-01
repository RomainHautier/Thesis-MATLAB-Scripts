% This function imports the data stored in the VBZ_PSI.DAT file & Caculates
% the inflow norm at both the Main and Tail rotors for the Freewake cases.
function [VBZPSImatout, VBZpsi, inflownorm] = importVBZPSI(nr,speed,rotcase)
% Note that the '2rotcase' input variable must be left empty in case the
% the standalone rotors are being considered.

%% Initialize variables.
if nr == 1
    rotor = {'MR','TR'};
else
    rotor = {'2'};
end

for r = 1:length(rotor)
    if nr == 1
        filename = ['C:\Users\hauti\heli-fwmrtr2\MyCasesNew-8turns\' rotor{r} '\t' num2str(speed) '\VBZ_PSI.dat'];
    elseif nr == 2
        filename = ['C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\8turns\d' num2str(speed) '\' rotcase 'VBZ_PSI.dat'];
    else % let nr = another number than 1 or 2 to plot the MR inflow from HeliUM. The nr is then set back to 1 to allow plotting thereafter.
        nr = 1;
        filename = ['C:\Users\hauti\heliumx2\HeliUM.sample\uh60a_fw_' num2str(speed) 'kt\VBZ_PSI.dat'];
    end
    
    ns = 32;
    np = 36;
    rad = [8.18, 1.676];
    om = [27, 124.6];
    
    %% Format for each line of text:
    %   column1: double (%f)
    %	column2: double (%f)
    %   column3: double (%f)
    %	column4: double (%f)
    %   column5: double (%f)
    %	column6: double (%f)
    %   column7: double (%f)
    % For more information, see the TEXTSCAN documentation.
    formatSpec = '%12f%13f%13f%13f%13f%13f%f%[^\n\r]';
    
    %% Open the text file.
    fileID = fopen(filename,'r');
    
    %% Read columns of data according to the format.
    % This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
    
    %% Close the text file.
    fclose(fileID);
    
    %% Create output variable
    VBZPSImat = [dataArray{1:end-1}];
    VBZPSImatout(:,:,r) = VBZPSImat;
    VBZpsi(1:length(VBZPSImat(:,1)),1:length(VBZPSImat(1,:))) = VBZPSImat;
    VBZ(1,:) = VBZpsi(:,6);

    if nr == 1

        s = ns:ns:length(VBZ(1,:));
        k=1;
    
        for p = 1:np
            inflownw(1:ns,p) = VBZ(1,k:s(p));
            k = k+ns;
        end
    
        norm = sum(sum(inflownw(:,:)));
        normmps=-norm/(ns*np);
        inflownorm(r,1) = normmps/(om(r)*rad(r));

    elseif nr == 2

        VBZ_ph = VBZ; % Creating a placeholder for manipulation.
        VBZ = [];
        VBZ(1,1:length(VBZ_ph)/2) = VBZ_ph(1,1:length(VBZ_ph)/2);
        VBZ(2,1:length(VBZ_ph)/2) = VBZ_ph(1,length(VBZ_ph)/2+1:end);

        s = ns:ns:length(VBZ(1,:));
        
        % Iterating over the two row of VBZ, representing the data from the
        % MR and TR (1st and 2nd rotor respectively)
        for z = 1:2
            k = 1;
            for p = 1:np
                inflownw(1:ns,p) = VBZ(z,k:s(p));
                k = k+ns;
            end
        
            norm = sum(sum(inflownw(:,:)));
            normmps=-norm/(ns*np);
            inflownorm(z,1) = normmps/(om(z)*rad(z));
        end
    end

end