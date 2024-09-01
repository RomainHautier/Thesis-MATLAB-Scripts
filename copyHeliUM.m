clear;
clc;
% close all;

copy_speed = 25;
originalPath = {'C:\Users\hauti\heliumx2','C:\Users\hauti\MATLAB Drive\Research Project'};
filesToCopy = {'helium.exe', 'inflowcorrectiondata.data'};


%% Creating new instances of the HeliUM runs to run the new inflow corrected cases into.

copy = 'n'

if copy == 'y'
    copy_speed = 25;
    
    for i = 1:4
        oldFolderName = ['C:\Users\hauti\heliumx2\HeliUM.sample\uh60a_fw_' num2str(copy_speed) 'kt']
        newFolderName = ['C:\Users\hauti\heliumx2\HeliUM.sample\uh60a_fw_' num2str(copy_speed) 'ktc' num2str(i)]
        
        if ~isfolder(newFolderName)
            mkdir(newFolderName)
        end
        newFolderName
        for x = 1:length(filesToCopy)
            srcFile = fullfile(originalPath{x}, filesToCopy{x})
            destFile = fullfile(newFolderName)
            copyfile(srcFile, destFile);
            fprintf('Copied "%s" to "%s"\n', srcFile, destFile);
        end

        copyfile(oldFolderName, newFolderName,'f');
        fprintf('Copied "%s" to "%s"\n', srcFile, destFile);
    end

    
end

%% Copying the helium.exe instance & the inflowcorrectiondata.data file into the new folders.

for j = 1:length(copy_speed)
    for i = 1:4
        for x = 1:length(filesToCopy)
            newFolderName = ['C:\Users\hauti\heliumx2\HeliUM.sample\uh60a_fw_' num2str(copy_speed(j)) 'ktc' num2str(i)]
            srcFile = fullfile(originalPath{x}, filesToCopy{x})
            destFile = fullfile(newFolderName)
            copyfile(srcFile, destFile);
            fprintf('Copied "%s" to "%s"\n', srcFile, destFile);
        end
    end
end

