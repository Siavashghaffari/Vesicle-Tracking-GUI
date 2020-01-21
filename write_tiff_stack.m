function write_tiff_stack (filename, stack)
%
% Function to write a x,y,slice stack into an x,y,slice tiff stack
%
% Output file will be 8 or 16 bit tiff (depending on stack type),
% uncompressed.
%
% Inputs:
%   -filename = name to give output file, name should include the .tif
%    extention
%   -stack = stack to write to file
%
% Output:
%   -A tiff formatted file, with the name [filename]
%
%

if nargin ~= 2
   error 'requires file name and input stack'
end

if ndims(stack)>3
    error 'write_tiff_stack can only handle single images and single-dimention stacks'
end

if ndims(stack) == 2
    imwrite(stack(:,:),filename, 'tif', 'Compression','none','WriteMode','overwrite');
else
 wr_mode = 'overwrite';
 for i=1:max(size(stack(1,1,:)))
     imwrite(stack(:,:,i),filename, 'tif', 'Compression','none','WriteMode',wr_mode);
     wr_mode = 'append';
 end
end
end


