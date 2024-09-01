function data = importEXCEL(filepath)
    %    Workbook: C:\Users\hauti\MATLAB Drive\Research Project\FWHELsheet.xlsx
    %    Worksheet: Results
    %
    % Auto-generated by MATLAB on 05-Jul-2024 17:03:00

    %% Set up the Import Options and import the data
    opts = spreadsheetImportOptions("NumVariables", 22);

    % Specify sheet and range
    opts.Sheet = "Results";
    opts.DataRange = "A5:V25";

    % Specify co% Specify column names and types
    opts.VariableNames = ["Velocitykts", "CT", "VarName3", "HelicopterVelocityVectorfts", "VarName5", "VarName6", "MRPowerHP", "MRThrustLBS", "TRInflow", "CollectiveControl", "VarName11", "ReadFromAATAILR", "VelocityAtTRInBodyFramefts", "VarName14", "VarName15", "VelocityAtTRInTRFramefts", "VarName17", "VarName18", "MUTR", "VarName20", "VarName21","TTR"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];


    % Import the data
    FWHELsheet = readtable(filepath, opts, "UseExcel", false);

    %% Convert to output type
    data = table2array(FWHELsheet);

    %% Clear temporary variables
    clear opts
end