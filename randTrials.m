function out = randTrials(n,varargin)
% This function randomizes presentations of stimuli for the linearTDSR function

%INPUTS:
%n is the number of presentations for each stimulus sequence

%the varargin are the unique stimulus sequences to be included 
    %each stimulus is an [s X D] matrix. (s = # of states within a single
    %stimulus sequence, D = # of features 

%OUTPUT: [n*s*(nargin - 1) X D] of randomized stimulus presentations
    
%Note that zeros(1 X D) are added between each presentation

out = [];
ncue = nargin - 1;
curr = 0;
for i = 1:n

    
    A = randperm(ncue,ncue);
    
    for j = 1:ncue
        
        str = size(varargin{A(j)},1);
        
        for k = 1:str
            curr = curr + 1;
            out(curr, :) = varargin{A(j)}(k,:);
        end
        curr = curr + 1;
        out(curr, :) = zeros(1, size(varargin{A(j)},2));
    end
end    
