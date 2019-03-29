% Script containing code to perform UP/DOWN state detection on local field potential recordings
% and analyze firing properties of simultaneous juxtasomal recordings of cortical interneurons.
% Written by Valentina Pasquale (2017)
% Contact: valentina.pasquale@iit.it
%%
% sourceFolder must contain the folders of all experiments and a file called
% "files.txt" containing the list of all files containing parameters (one for each experiment)
sourceFolder = 'D:\Valentina\current_activities\Fellin\full_dataset\UP&DOWNstateDet\Interneuron Analysis\PV-Positive';
txtOutputAllExpFolder = 'D:\Valentina\current_activities\Fellin\full_dataset\UP&DOWNstateDet\Interneuron Analysis\statisticalAnalysisData_PV_20160302';
figOutputAllExpFolder = 'D:\Valentina\current_activities\Fellin\full_dataset\UP&DOWNstateDet\Interneuron Analysis\figOutputData_PV_20160302';
listFilename = [sourceFolder, filesep, 'files_test.txt'];
fileID = fopen(listFilename,'r');
fileList = textscan(fileID,'%s');
fclose(fileID);
fileList = fileList{1,1};
for ii = 1:size(fileList,1)
    parameterFile = [sourceFolder,fileList{ii,1}];
    [outputFolder] = UP_DOWN_DET_analyzeJuxta(parameterFile);
    if ~isempty(outputFolder)
        [~,~,txtFileNames] = dirr(outputFolder,'\.txt\>','name','isdir','0');
        for jj = 1:length(txtFileNames)
            try
                copyfile(txtFileNames{1,jj},txtOutputAllExpFolder)
            catch ME
                error(ME.identifier,ME.message);
                return
            end
        end
        [~,~,figFileName_lockStrength] = dirr(outputFolder,'\lockStrengthVStau.pdf\>','name','isdir','0');
        [~,~,figFileName_phaseFiringDistrib] = dirr(outputFolder,'\phaseFiringDistrib.pdf\>','name','isdir','0');
        try
            copyfile(figFileName_lockStrength{1},figOutputAllExpFolder)
            copyfile(figFileName_phaseFiringDistrib{1},figOutputAllExpFolder)
        catch ME
            error(ME.identifier,ME.message);
            return
        end
    else
        return
    end
end