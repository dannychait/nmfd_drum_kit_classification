function [audioFiles] = getFileNames(directory, ext)

    directoryFiles = dir(directory);

    j = 1;
    
    ext = ['\.' ext];
    
    audioFiles = { };

    for i=1:length(directoryFiles)
        if(regexpi(directoryFiles(i).name, ext))
            audioFiles{j} = directoryFiles(i).name;
            j = j + 1;
        end
    end
end

