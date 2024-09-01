function [inflowmat, inflownorm] = importINFLOW(speed, filepath,r,filend)

    inflowmat = zeros(1800,21,r);
    np = 36;
    dnp = 9;

    for i = 1:length(speed)

        %% Initialize variables.
        filename = [filepath num2str(speed(i)) '\timectcq.dat'];

        if r == 2
           filename = [filepath num2str(speed(i)) '\' filend '\timectcq.dat'];
        end
        % filename = ['C:\Users\hauti\heli-fw\MyCases\t' num2str(speed(i)) '\timectcq.dat'];

        %% Format for each line of text:
        %   column5: double (%f)
        % For more information, see the TEXTSCAN documentation.
        formatSpec = '%*44s%14f%[^\n\r]';

        %% Open the text file.
        fileID = fopen(filename,'r');

        %% Read columns of data according to the format.
        % This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
        dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);

        %% Close the text file.
        fclose(fileID);

        %% Post processing for unimportable data.
        % No unimportable data rules were applied during the import, so no post processing code is included. To generate code which works for unimportable data, select unimportable cells in a file and regenerate the script.

        %% Create output variable
        inflowvec = [dataArray{1:end-1}];

        % Only taking the inflow norm from the last rotor rotation --> here
        % the last 36 values of the inflow values are taken, thus the last
        % values from the last 360°. 
        inflownorm(i) = mean(inflowvec(length(inflowvec)-np+1:end));

        if r == 1
            inflownorm(i) = mean(inflowvec(length(inflowvec)-np+1:end));
        else
            intervals = length(inflowvec)/dnp;
            rotor_indices = [];
            rotor1_indices = [];
            count = 0;
            count1 = 0;
            for z = 1:intervals/2
                rotor_indices = [rotor_indices (1:9)+dnp*2*count1];
                rotor1_indices = [rotor1_indices (10:18)+dnp*2*count1];
                count1 = count1 + 1;
            end
            
            inflownorm(i,1) = mean(inflowvec(rotor_indices));
            inflownorm(i,2) = mean(inflowvec(rotor1_indices));
        end
        
        % Note that the total length of the inflow vector is nfw/4*np.
        if length(inflowvec) ~= length(inflowmat(:,1))
            inflowvec(length(inflowmat(:,1))) = 0;
        end

        inflowmat(:,i) = inflowvec;

    % if speed == 80
    %     inflowverif = [dataArray{1:end-1}];
    %     averageintot = mean(inflowverif);
    %     averagein = mean(inflowverif(length(inflowverif)-np+1:end));
    % end

        %% Clear temporary variables
        clearvars filename formatSpec fileID dataArray ans;
end