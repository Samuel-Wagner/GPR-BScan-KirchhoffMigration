function I = Bscan_migration_v3(B, h, er, t, x, tx_pos, rx_pos, xQ, zQ, progress)
% INPUTS
% B - pre-processed B-scan. should be Nt * Nx 
% h - antenna height off of ground
% er - relative permittivity of ground
% t - time vector (sampling rate)
% x - relative movement of the array
% tx_pos - the x/y position tx
% rx_pos - the x/y position rx
% xQ - calculation domain for x
% yQ - calculation domain for y
% zQ - calculation domain for z
% progress - display progress bar?
%
% OUTPUTS
% I - migrated data for [xQ yQ zQ] calculation domain
%
% !!!! NOTE!!!! The proper Kirchhoff Migration algorithm
% has several weighting factors, including 
% 1/sqrt(R) and cos(theta). Additionally, it uses the half-derivative of 
% B instead of just B. I have neglected all of these for
% 1) speed (slight increase)
% 2) aesthetic (half-derivative looks worse!)
% 3) they don't hurt the end product much.

% for a given box, calculate all possible time-of-arrival indices
c   = 3e8;
NxQ = numel(xQ);
NzQ = numel(zQ);

x_lims = 0.5;                   % how much +/- in scan x of the tx/rx midpoint to include in calculation.
dx     = x(2)-x(1);             % step size of scan x
xc     = -x_lims:dx:x_lims;     % chopped x for evaluation
Nxc    = numel(xc);             % number of x "positions" for the array



%% PART 1 ----
% calculate the indices of time
TTD_inds = zeros(Nxc,NzQ,'uint32');

% each index corresponds to the time index at which an antenna positioned
% at x=0 will see a target at xQc,yQ,zQ for each receiver
if(progress)
    progressbar('Sliding B-scan migration')
end

% if migrating the same data repeatedly, this section
% could be done once for all data.

for ii = 1:Nxc
    % select xq
    xc_q = xc(ii);  

        for jj = 1:NzQ
            % select zq
            zQ_q = zQ(jj);

            % calculate the GPR angles for the given Tx
            [theta_a_tx, theta_g_tx] = GPR_transmission_angles_v4(...
                er,...              % relative perm. of ground
                h,...               % height off of ground
                tx_pos + xc_q,...   % x-value of antenna
                0,...               % y-value of antenna
                0,...               % query x-point
                0,...               % query y-point
                zQ_q);              % query z-point

            % true time delay - from TX to TARG (xq,yq,zq)
            TTD_tx = (h*sec(theta_a_tx) + zQ_q*sec(theta_g_tx)*sqrt(er))/c; 

            % calculate the GPR angles for the given Rx
            [theta_a_rx, theta_g_rx] = GPR_transmission_angles_v4(...
                er,...              % relative perm. of ground
                h,...               % height off of ground
                rx_pos + xc_q,...   % x-value of antenna
                0,...               % y-value of antenna
                0,...               % query x-point
                0,...               % query y-point
                zQ_q);              % query z-point
            
            % true time delay - from TX to RX (combined with 
            TTD = (h*sec(theta_a_rx) + zQ_q*sec(theta_g_rx)*sqrt(er))/c + TTD_tx;
            
            % find the index of time corresponding to this
            TTD_inds(ii,jj) = min1(t,TTD);                    
        end         %z-iter
end                 %x-iter

%% PART 2 ----
% "slide" the indices previously calculated over the entire scan to
% migrate.sliding is much, much faster. - but more complex
% TTD_inds - true-time-delay inds - will trace a hyperbola
% in (z,x) for xa=0

I = zeros(NzQ,NxQ,'single');
for ii = 1:NxQ
    xq = xQ(ii);
    [~,x0_ind]=min(abs(x-xq));
    
    %determine the limits of xa that we will look at.
    xa_lower = max(min(x),x(x0_ind)-x_lims);
    xa_upper = min(max(x),x(x0_ind)+x_lims);

    % find out what indices in x to subset
    size_left  = floor((x(x0_ind)-xa_lower) ./ dx);
    size_right = floor((xa_upper-x(x0_ind)) ./ dx);
    
    ind_left_TTD = floor(Nxc/2)-size_left;
    ind_right_TTD = floor(Nxc/2)+size_right;
    
    [~,xa_lower_ind]=min(abs(x-xa_lower));
    
    C11 = B(:,xa_lower_ind + (0:(size_left+size_right)));
    I2  = (1:size(C11,2))';
    
    
    for jj = 1:NzQ
        % now looking at xq,zq.
        TTD_mat  = squeeze(TTD_inds(ind_left_TTD:ind_right_TTD,jj));
        linInds  = sub2ind(size(C11),TTD_mat(:),I2);
        I(jj,ii) = sum(C11(linInds));
    end
    if(progress)
        progressbar(ii/NxQ);
    end
end

% take the absolute value of the image. Don't discriminate against negative
% pulses!
I = abs(I);