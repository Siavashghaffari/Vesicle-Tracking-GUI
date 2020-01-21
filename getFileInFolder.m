function [files,folders] = getFileInFolder(Folder, Extention)
% returns a list of the files and the folders in the location specified by the
% input argument 'Folder'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
files=[];
folders=[];

% Make sure 'folder' actually exists
if not(exist(Folder)==7)
    error([Folder ' seems not to exist']);
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

