function [vslice, hslice] = find_image_resolution_slices(I)
% this function will find the maximum in the target image and extract the 
% horizontal and vertical slices.
% I - image
% 
% vslice - vertical slice
% hslice - horizontal slice
%
% vslice and hslice will be padded with nan's to be the same size

% find dimensions
[Nz, Nx] = size(I);

% before we calculate anything, figure out if we need to pad.
pad_vert = false;
pad_horz = false;

if(Nx > Nz)
    pad_vert = true;
elseif(Nz > Nx)
    pad_horz = true;
end

% find the maximum of the image
[~,max_ind]     = max(I(:));                
[max_r, max_c]  = ind2sub(size(I),max_ind); %maximum is indexed as I(max_r,max_c)

vslice = I(:,max_c)';
hslice = I(max_r,:);

if(pad_vert)
    vslice = cat(2,vslice,nan.*ones(1,Nx-Nz));
elseif(pad_horz)
    hslice = cat(2,hslice,nan.*ones(1,Nz-Nx));
end
