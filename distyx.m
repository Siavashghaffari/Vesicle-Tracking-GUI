function [D]=distyx(y,x)

% function: [D]=euclideandist(i,x)
%
% Calculates the Euclidean distances between the position y and all objects in x	 
%								    
% Input: 
% y - an position [x,y,...] 
% x - data matrix (m,n); m-objects, n-variables	    
%                                                                 
% Output: 
% D - Euclidean distance (m,1)

[m,n]=size(x);
if n ~= length(y)
    D = [];
    error (['Position Y and Table X do not have the same dimentionality.  Position Y is ' int2str(length(y)) 'D, while table X is ' int2str(n) 'D.']);
end

D=sqrt(sum((((ones(m,1)*y)-x).^2)'));

if n==1
   D=abs((ones(m,1)*y-x))';
end

end %end function
