function [files,folders] = getFileInFolder2(Folder, Extention)
% returns a list of the files and the folders in the location specified by the
% input argument 'Folder'
% same as getFileInFoder.m, but handles specific file extentions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
files=[];
folders=[];

% Make sure 'folder' actually exists
if not(exist(Folder)==7)
    error([Folder ' does not exist.']);
end

f=dir(Folder);

% set counters
x=1;
y=1;

while x <= length(f)
    if f(x).isdir
        if ~strcmp(f(x).name,'.') && ~strcmp(f(x).name,'..')
            folders{y}=fullfile(Folder,f(x).name);
            y=y+1;
        end
        f(x)=[];
    elseif strcmp(f(x).name,'.') || strcmp(f(x).name,'..')
        f(x)=[];
    else        
        files{x}=f(x).name;
        x=x+1;
    end
end

%remove mac's idiotic additional files
if length(files)~=0
    fNum = 0;
    for i=1:length(files)
        tName = char(files(i));
        if tName(1) ~= '.'
            fNum = fNum + 1;
            tList(fNum) = files(i);
        end
    end
    if fNum == 0
      files = [];
    else
      files = tList;
    end
end


%reduces list to files present in folder with desired extention
if ~isempty(Extention) && length(files)~=0
    file_count = size(files,2);
    nfiles = 0;
    for k=1:file_count
        %generate list of files to process
        CurFile = char(files(k));
        s=strfind(CurFile,Extention);
        if ~(isempty(s))
            nfiles = nfiles+1;
            fileList(nfiles,:) = cellstr(CurFile);
        end %endif
    end %end for (k)
    if nfiles == 0
        fileList = char([]);
    end
    clear files;
    files(:,:) = fileList(:,:);
end

end %end function
