function [ outputMaps ] = filter2DMatrices_021521( inputMaps, filterSize, filterStd)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% 090518 I guess this is the gaussian smoothing mentioned by Dr, Nitz.
% however in that case the center of gaussian should be 0.5-0.8 or so
%2D ratemap smoothing function **************************************

if ~iscell(inputMaps)
    inputMaps={inputMaps};
end

halfNarrow = filterSize/2;% -5:5 is too large for 19*15 or 20:12 bins. 1/2 of the boxs are considered as possible regions that mice can go...- -
narrowStdev = filterStd;

narrowGaussian=fspecial('gaussian',halfNarrow*2,narrowStdev);
narrowGaussianNormed=narrowGaussian;

outputMaps=cell(length(inputMaps),1);
for i = 1:length(inputMaps)
    outputMaps{i} = nanconv(inputMaps{i},narrowGaussianNormed, 'nanout');
end

if length(inputMaps)==1
    outputMaps=[];
    outputMaps=nanconv(inputMaps{i},narrowGaussianNormed, 'nanout');
end
end


