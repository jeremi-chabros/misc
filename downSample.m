function dsMatrix = downSample(mat, option, dsFactor, method, vst)
% Description: Downsamples matrix bin-wise with specified method.
%----------
% INPUT
% mat      - [n x t]  matrix to be downsampled
% dsFactor - [scalar] downsampling factor
% method   - [string] method to be used: 'mean' for bin-wise mean,
%                     'sum' for within-bin sum
% vst      - [bool]   optional flag specifying the variance-stabilizing
%                     transformation; 1 to use, 0 to ignore
%----------
% OUTPUT
% dsMatrix - [n x t/dsFactor] downsampled matrix
%----------
% @author JJChabros (jjc80@cam.ac.uk), January 2021

if ~exist('vst', 'var')
    vst = 0;
end

if strcmp(option, 'nbins')
    dsFactor = length(mat)/dsFactor;
end

dsMatrix = zeros(size(mat,1),round(length(mat)/dsFactor));

switch method
    
    case 'mean'
        for i = 1:size(mat, 1)
            dsMatrix(i,:) = nanmean(reshape([mat(i,:); nan(mod(-numel(mat(i,:)),dsFactor),1)],dsFactor,[]));
        end
        
    case 'sum'
        for i = 1:size(mat, 1)
            dsMatrix(i,:) = nansum(reshape([mat(i,:); nan(mod(-numel(mat(i,:)),dsFactor),1)],dsFactor,[]));
        end
        
end

if vst
    dsMatrix = sqrt(dsMatrix);
end

end

