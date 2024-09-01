% This function copies the content of the FWGEOM.DAT,.... files from a
% previous 2 rotor solution to transform them as input for another 2 rotor
% case.

% Note that this function only supports the files to be copied and saved in
% another filepath to be in the same directory of 2 rotor cases, i.e. 4, 8
% or 16 turns here (following the way I have arranged my files).

function state = DATtoDATAconverter(filepath, indirectory, outdirectory, folder)
    % Defining the filenames of the files to be copied & transformed into
    % .DATA files to be saved in the desired trajectory.
    infiles = {'FWGEOM.DAT','FWG_b.DAT','flap.DAT'};
    outfiles = {'IWGEOM.DATA', 'IWG_b.DATA','flap.DATA'};
    
    for i = 1:length(infiles)
        % Define the filepath to the input directory.
        infilepath = [filepath '\d' indirectory '\' folder];
        % Define the filepath to the output directory.
        outfilepath = [filepath '\d' outdirectory];
        
        % Full path to the input file.
        infile = fullfile(infilepath, infiles{i})
        
        % Full path to the input file.
        outfile = fullfile(outfilepath, outfiles{i})

        % Open the input file for reading
        fid_in = fopen(infile, 'r');
        if fid_in == -1
            error('Could not open input file: %s', infile);
        end
        
        % Read the content of the input file
        content = fread(fid_in, '*char')';
        
        % Close the input file
        fclose(fid_in);
    
        % Open the output file for writing
        fid_out = fopen(outfile, 'w');
        if fid_out == -1
            error('Could not open output file: %s', outfile(i));
        end
        
        % Write the content to the output file
        fwrite(fid_out, content, 'char');
        
        % Close the output file
        fclose(fid_out);

    end
    
    fid = fopen([filepath '\d' outdirectory '\converted.txt'], 'wt');
    fprintf(fid, '%s', ['The wake files for ' outdirectory 'kts have been converted on the ' datestr(now) '.' ...
        'This was realised by combining the solution for the 2 rotor case ran at a speed of ' indirectory 'kts in the ' folder ' folder.']);
    fclose(fid);
    
    state = 'Done';
end