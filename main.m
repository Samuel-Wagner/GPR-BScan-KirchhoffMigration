% Migration of a GPR B-scan
% Samuel Wagner @ UCDavis ECE / MML , Aug. 9, 2021
% - This is an example driver to using the "Bscan_migration_v3.m" function
% - The calls to migrate are hidden behind "Scan" object methods
% - A number of datasets are included for examples
% - Each dataset is a pulsed GPR dataset with three targets
%   "xxps_y.mat": xx is the sigma value (o) for a differential (1st gauss.
%   derivative) input pulse, in ps. y is the dataset ID
% 
%
% to use, run the script. an example of migration with the Scan object
% and migration using the function manually are included
%
%
% utils\progressbar.m is written by Steve Hoelzer at the Mathworks file
% exchange https://www.mathworks.com/matlabcentral/fileexchange/6922-progressbar
%
%  
%
% housekeeping
addpath ('.\utils\');
fixit;
prettygraphs;

%% read in data
data_source_dir = ".\data\";     % directory of data
data_source_fn  = "50ps_1.mat"; % exact filename of scan
perform_timing_test = 1;         % do you want to do a timing test (time_taken vs. N)? set to 0 if no

% load in data. the data has been saved as a "Scan" object.
% it will be loaded in as a variable "S"
load(data_source_dir+data_source_fn);

% perform a background removal
L = 1;          % how many singular values to remove. 1 for low. BW pulses, 2 for hi. BW pulses
                % if L is too high, the image will form horizontal lobes!!
S.BKGR_PCA(L);  % call the PCA removal


% display the data
figure(); colormap(gray(2^12));
subplot(121);
imagesc(S.x,S.t.*1e9,S.Data); hold on;
colorbar;
xlabel('x (m)');
ylabel('Time (ns)');
title('Raw Data');

subplot(122);
imagesc(S.x,S.t.*1e9,S.DataBKGR); hold on;
colorbar;
xlabel('x (m)');
ylabel('Time (ns)');
title('BKGR Data');

%% set migration variables
% some variables will already be set... but I don't remember if they are
% set, or correct, and worth changing.

h       = 21.59e-2; % height of antenna off of ground
tx_dx   = 3e-2;     % distance of tx from array x "center"
rx_dx   = -3e-2;    % distance of rx from array x "center"
er      = 5;        % ground permittivity

dt      = 1.3e-9;   % where is t0, really?
compensate_dt = 1;  % do you want to compensate for t0? (recommended)

% set the variables
S.h      = h;   
S.tx_pos = tx_dx;
S.rx_pos = rx_dx;
S.er     = er;
S.t      = S.t - compensate_dt*dt; 

% the "calculation" vectors.
Nx = 512;
Nz = 256;
xC = linspace(0,max(S.x), Nx);
zC = linspace(0, 40e-2, Nz);
display_progress_bar = 1;

% perform migration
S.migrate(xC, zC, display_progress_bar);

% normalize the image for display purposes
S.Image = normalize(S.Image);

%% display the migrated image next to the BKGR B-scan

DR = 10; % dynamic range (in dB10) with which to display image

% data bkgr
figure(); 
subplot(121);
imagesc(S.x,S.t.*1e9,S.DataBKGR); hold on;
colorbar;
xlabel('x (m)');
ylabel('Time (ns)');
title('BKGR Data');
set(gca,'colormap', gray(2^12));

%image
subplot(122);
imagesc(S.xC,S.zC.*1e2,10.*log10(S.Image)); hold on;
colorbar;
xlabel('x (m)');
ylabel('z (cm)');
title('Migrated Image');
set(gca,'colormap', jet);
caxis([-DR 0]);
colorbar;

%% perform a matched filter and matched filter migration.
% obtain the input pulse from the filename. if not using given dataset,
% this section will definitely cause many errors

% matched filter migration is preferred by some because it produces only 
% "single" dot targets for point targets. but it has a lower dynamic range
% (effectively adds some inter-target clutter). Typically, view non-MF
% migrations with a DR of ~ -10 dB, and a MF-migration with 6-3 dB DR.

split_fn = strsplit(data_source_fn, 'p');
o   = str2double(split_fn(1)).*1e-12;
Ts  = S.t(2)-S.t(1);
N   = 400;
t   = ((-N/2+1):(N/2)).*Ts;

% simulate the DOG pulses
base_gauss      = normalize(exp(-1/2.*(t./o).^2));      % 0th derivative
differential    = normalize(gradient(base_gauss));      % 1st derivative
ricker          = normalize(gradient(differential));    % 2nd derivative
thirddiff       = normalize(gradient(ricker));          % 3rd derivative

% we take the ricker as the matched wfm.
matched_wfm      = ricker;  % selection
flip_matched_wfm = 1;       % does it need +/- flip? probably, yes.

if(flip_matched_wfm)
   matched_wfm = -matched_wfm;
end

take_positive_xcorr_only = 1;

% perform matched filter
S.matched_filter(matched_wfm,take_positive_xcorr_only);

% migrate the matched-filter data
S.migrateMF(xC,zC,1);

% normalize it for disp.
S.ImageMF = normalize(S.ImageMF);

%% display it all
DR = 10;
DR_matched = 6;

% normal data bkgr
figure()
subplot(221);
imagesc(S.x,S.t.*1e9,S.DataBKGR); hold on;
colorbar;
xlabel('x (m)');
ylabel('Time (ns)');
title('BKGR Data');
set(gca,'colormap', gray(2^12));

% matched-filtered data
subplot(222);
imagesc(S.x,S.t.*1e9,S.DataMF); hold on;
colorbar;
xlabel('x (m)');
ylabel('Time (ns)');
title('BKGR/MF Data');
set(gca,'colormap', gray(2^12));

% normal image
subplot(223);
imagesc(S.xC,S.zC.*1e2,10.*log10(S.Image)); hold on;
colorbar;
xlabel('x (m)');
ylabel('z (cm)');
title('Migrated Image');
set(gca,'colormap', jet);
caxis([-DR 0]);
colorbar;

% matched-filtered image
subplot(224);
imagesc(S.xC,S.zC.*1e2,10.*log10(S.ImageMF)); hold on;
colorbar;
xlabel('x (m)');
ylabel('z (cm)');
title('Migrated Image (MF)');
set(gca,'colormap', jet);
caxis([-DR_matched 0]);
colorbar;

%% plot the resolution slices from the maximum target
[vslice, hslice] = find_image_resolution_slices(S.Image);
[vslice_MF, hslice_MF] = find_image_resolution_slices(S.ImageMF);

% take out the nan's inserted by find_image_resolution_slices()
vslice(isnan(vslice))       =[];
hslice(isnan(hslice))       =[];
vslice_MF(isnan(vslice_MF)) =[];
hslice_MF(isnan(hslice_MF)) =[];

% vertical slices
figure();
subplot(121); hold on; grid on;
plot(S.zC.*1e2, vslice,'k');
plot(S.zC.*1e2, vslice_MF,'b');
xlabel('Depth (cm)');
ylabel('Image mag, linear');
legend('Normal KM','MF/KM');
title('Vertical Slices');

% horizontal slices
subplot(122);hold on; grid on;
plot(S.xC.*1e2, hslice,'k');
plot(S.xC.*1e2, hslice_MF,'b');
xlabel('x (m)');
ylabel('Image mag, linear');
legend('Normal KM','MF/KM');
title('Horizontal Slices');

%% do a time-test using manual migration
% this is an example of how to use the migration function manually
if(perform_timing_test)
    % get data
    D = S.DataBKGR;
    er = S.er;
    tx_pos = S.tx_pos;
    rx_pos = S.rx_pos;
    h = S.h;
    t = S.t;
    x = S.x;

    % vector of "N" to measure time taken
    Nt = 30;
    pixel_widths = logspace(log10(2e-2), log10(0.5e-3), Nt);
    times = zeros(1,Nt);
    progressbar('Timing test using manual migration');
    for ii = 1:Nt
        % function prototype
        % I = Bscan_migration_v3(B, h, er, t, x, tx_pos, rx_pos, xQ, zQ, progress)  
        xC = 0:pixel_widths(ii):max(x);
        zC = 1e-3:pixel_widths(ii):40e-2;
        tic;
        I = Bscan_migration_v3(D, h, er, t, x, tx_pos, rx_pos, xC, zC, 0);
        times(ii)=toc;
        progressbar(ii/Nt);
    end
    
    figure(); 
    semilogx(pixel_widths.*1e3,times,'kx-'); hold on; grid on;
    xlabel('Pixel width (mm)');
    ylabel('Time taken (s)');
    set(gca,'XTickLabel',get(gca,'xtick'));
    aa=gca;
    aa.GridColor=[0.1 0.1 0.1];
    aa.GridAlpha=0.6;
    aa.MinorGridAlpha=0.6;
end