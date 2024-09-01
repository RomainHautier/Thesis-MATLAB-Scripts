clear;
clc;

% This script is aimed at saving the data from the runs into the desired
% folders.
turns = 8;
speed = 95;
sol_case = 'fb'; % This defines the solution case used,


directory = ['C:\Users\hauti\heli-fwmrtr2\MyCases2rotFW\' num2str(turns) 'turns\d' num2str(speed)]
enddirectory = [directory '\' sol_case]
% Define the specific file names to be copied
% filesToCopy = {'user.INPUT', 'VBX.DAT', 'VBY.DAT', 'VBZ.DAT', 'VBZ_PSI.DAT','FWGEOM.DAT', 'FWG_b.DAT','flap.DAT','timectcq.DAT','clhist.DAT','lift.DAT', 'timeflap.DAT', 'timeflight.DAT','vrz.DAT'};
filesToCopy = {'user.INPUT', 'VBX.DAT', 'VBY.DAT', 'VBZ.DAT', 'VBZ_PSI.DAT','FWGEOM.DAT', 'FWG_b.DAT','flap.DAT','timectcq.DAT','clhist.DAT', 'timeflap.DAT', 'timeflight.DAT','vrz.DAT'};
if sol_case(1) == 'y'
    filesToCopy{length(filesToCopy)+1} = 'timeflightmodified.DAT';
end

% Loop through each file and copy it to the output directory
for i = 1:length(filesToCopy)
    % Construct full file paths for input and output
    infile = fullfile(directory, filesToCopy{i});
    outfile = fullfile(enddirectory, filesToCopy{i});
    
    % Copy file and overwrite if it already exists
    copyfile(infile, outfile, 'f');
end