classdef Scan < handle
    properties 
        Data     {mustBeNumeric}; % Raw Data of time by x (rows - time, columns- x)
        DataBKGR {mustBeNumeric}; % BKGR data - background removed
        DataFK   {mustBeNumeric}; % frequency (time/spatial) version of Data or DataBKGR
        DataMF   {mustBeNumeric}; % matched filter
        DataENV  {mustBeNumeric}; % envelope (in time) version of DataBKGR
        Image    {mustBeNumeric}; % migrated image
        ImageMF  {mustBeNumeric}; % migrated image

        xC {mustBeNumeric}; %calculation points (x) for image
        zC {mustBeNumeric}; %calculation points (z) for image

        x {mustBeNumeric}; %x vector
        t {mustBeNumeric}; %time vector
        f {mustBeNumeric}; %frequency vector
        Nx {mustBeNumeric}; % number of x points
        Nt {mustBeNumeric}; % number of t points
        Ts {mustBeNumeric}; %time sampling period
        Fs {mustBeNumeric}; %frequency sampling rate
        Kx {mustBeNumeric}; %spatial sampling rate
        kx {mustBeNumeric}; %spatial frequency vector
        dx {mustBeNumeric}; %spatial sampling period
        matchedWFM {mustBeNumeric}; % waveform used for matched-filtering
        tx_pos {mustBeNumeric}; % position of tx antenna(s) 1D ONLY
        rx_pos {mustBeNumeric}; % position of rx antenna(s) 1D ONLY (can be worked around)

        h {mustBeNumeric}; %height of scan above ground
        er {mustBeNumeric} %relative permittivity of ground
       
        title string;      %title of scan
    end
   
    methods
        %constructor
        function obj = Scan(ScanData,xVec,tVec) 
            if(nargin>0)
                obj.Data=ScanData;

                obj.x  = xVec;
                obj.dx = xVec(2)-xVec(1);
                obj.Kx = 1/obj.dx;
                obj.Nx = numel(xVec);
                obj.kx = linspace(-obj.Kx/2, obj.Kx/2,obj.Nx);

                obj.t  = tVec;
                obj.Ts = tVec(2)-tVec(1);
                obj.Fs = 1/obj.Ts;
                obj.Nt = numel(tVec);
                obj.f  = linspace(-obj.Fs/2,obj.Fs/2,obj.Nt);
            end
        end
       
        %background removal using principal component analysis
        function obj = BKGR_PCA(obj,L)
            %L - number of singular values to remove, 2 is best, 1 is OK 
            [U, S, V] = svd(obj.Data,'econ');
            o = diag(S);
            obj.DataBKGR = U*diag([zeros(L,1); o(L+1:end)])*V';
        end
       
        % function to get the envelope signal
        function obj=envelope(obj,useBKGR)
            if(useBKGR)
                obj.DataENV = abs(hilbert(obj.DataBKGR));
            else
                obj.DataENV = abs(hilbert(obj.Data));
            end
        end
       
        %background removal using mean removal
        function obj = BKGR_MEAN(obj)
            obj.DataBKGR = obj.Data - mean(obj.Data,2);
        end
       
        %migrate the background-removed B-scan
        function obj = migrate(obj,xC,zC,progress)
            %xC - x points
            %zC - z points
            obj.Image = Bscan_migration_v3(obj.DataBKGR,...
                obj.h, obj.er, obj.t, obj.x,obj.tx_pos, obj.rx_pos, xC, zC, progress);
            obj.xC = xC;
            obj.zC = zC;
        end
       
        % function to "matched filter" the scan
        function obj = matched_filter(obj,wfm,positive)
            % wfm - waveform to use as reference
            % positive - do you want to throw out negative xcorrs?

            obj.DataMF = zeros(size(obj.DataBKGR));
            for k = 1:numel(obj.x)
                obj.DataMF(:,k)=conv(obj.DataBKGR(:,k),flip(wfm),'same'); 
            end
            obj.matchedWFM=wfm;

            if(positive)
                obj.DataMF(obj.DataMF<0)=0;
            end
        end
       
        % migrate matched filter data
        function obj=migrateMF(obj,xC,zC,progress)
            obj.ImageMF = Bscan_migration_v3(obj.DataMF,...
                obj.h, obj.er, obj.t, obj.x,obj.tx_pos, obj.rx_pos, xC, zC, progress);
            obj.xC = xC;
            obj.zC = zC; 
        end

        % generate frequency (f-k) representation
        function obj = generate_FK(obj,useBKGR)
            if(useBKGR)
                temp = fftshift(fft(obj.DataBKGR,[],1),1);
            else
                temp = fftshift(fft(obj.Data,[],1),1);
            end
            
            obj.DataFK = fftshift(fft(temp,[],2),2);
        end
    end
end