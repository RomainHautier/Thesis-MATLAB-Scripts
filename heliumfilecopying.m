clear;
clc;
close all;

% MATLAB script to copy and rename directories
% Define the base directory where the folders are located
baseDir = 'C:\Users\hauti\heliumx2\HeliUM.sample'; % or specify a directory, e.g., 'C:/myfolders'

speed = 50;
filesToCopy = {'helium.exe', 'inflowcorrectiondata.DATA'};
originalPath = {'C:\Users\hauti\heliumx2', 'C:\Users\hauti\MATLAB Drive\Research Project'};

for i = 1:3
    % Get a list of all folders matching the pattern 'uh60a_fw_80kt'
    folderName = fullfile(baseDir, ['uh60a_fw_' num2str(speed) 'ktc']);

    newFolderName = [folderName num2str(i)];
    % originalPath = 'C:\Users\hauti\heliumx2'
    newFolderPath = fullfile(baseDir, newFolderName);

    % Copy the folder to the new location
    copyfile(folderName, newFolderName);

    % Copying the new helium.exe instance and the inflow correction data
    % file to the directories.

    % for j = 1:length(filesToCopy)
    %     srcFile = fullfile(originalPath, filesToCopy{j});
    %     destFile = fullfile(newFolderName, filesToCopy{j});
    %     copyfile(srcFile, destFile);
    %     % fprintf('Copied "%s" to "%s"\n', srcFile, destFile);
    % end



end


% speed_hel_copy = [1 5:5:100];
% 
% for j = 1:length(speed_hel_copy)
%     newFolderName = ['C:\Users\hauti\heliumx2\HeliUM.sample\uh60a_fw_' num2str(speed_hel_copy(j)) 'ktc']
%     srcFile = fullfile(originalPath, 'helium.exe')
%     destFile = fullfile(newFolderName)
%     copyfile(srcFile, destFile);
%     % fprintf('Copied "%s" to "%s"\n', srcFile, destFile);
% end




% % Loop through each folder found
% for k = 1:length(folderList)
%     % Get the folder name
%     oldFolderName = folderList(k).name;
% 
%     % Construct the new folder name by appending 'c'
%     newFolderName = [oldFolderName '_c'];
% 
%     % Full paths for source and destination folders
%     oldFolderPath = fullfile(baseDir, oldFolderName);
%     newFolderPath = fullfile(baseDir, newFolderName);
% 
%     % Copy the folder to the new location
%     copyfile(oldFolderPath, newFolderPath);
% 
%     % Display message
%     fprintf('Copied and renamed folder "%s" to "%s"\n', oldFolderName, newFolderName);
% end
% 
% disp('All folders have been copied and renamed.');