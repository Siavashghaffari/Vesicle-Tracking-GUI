function Stack1 = LoadStack(Im1)
%
% Macro to load stack from file
%
%  Input is an image stack.  Stack is assumed to be a 3D stack structured
%  X,Y,Z/t
%
%  Output is an array of type unit16, with the same X/Y/Z dimentions as the
%  input image.
%
%

%% Error checking
if nargin ~= 1
    error('Incorrect number of input variables.  Input is a single image stack file')
end

%make sure stack exists
if ~(exist(Im1)==2)
    error(['Image file ' Im1 ' does not exist!'])
end

stackInfo=imfinfo(Im1);
num_slices = length(stackInfo);

%% Load Stack
% pre-allocate memory based on the x/y resolution and the number of frames.
Stack1=uint16(zeros(stackInfo(1,1).Height,stackInfo(1,1).Width, num_slices));

%load stack
for i=1:num_slices
    Stack1(:,:,i)=imread(Im1,i);
end

end
