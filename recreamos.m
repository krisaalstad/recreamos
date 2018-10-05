function rout=recreamos( sat, indir, aoi, dopar )
%RECREAMOS  REad, CRop, REproject, resAMple and mOSaic satellite imagery
% Inputs:
%   sat=Satellite keyword argument. Sentinel-2=Sentinel-2 L1C (TOA), Landsat_8=Landsat 8 L1 (TOA),
%       Landsat_8_SR=Landsat 8 L2 (SR), Landsat_7=Landsat 7 L1, 
%       RapidEye=RapidEye TOA, Planet=PlanetScope TOA.
%   indir=(Relative) path to the input directory that contains the satellite
%       image(s), either as one scene in the directory or multiple scenes in
%       seperate numbered (i.e. 1,2,3 etc...) subdirectories.
%   aoi=An "area of interest" structure that contains a rectangular
%       bounding box defining the area of interest, defined by UTM coordinate
%       vectors aoi.x and aoi.y, as well as the UTM zone.
%   dopar=Switch (1/0) if the resampling should be performed using a
%   parallel pool (1) or not (0). If the parallel switch is on then half
%   the available clusters on the local parallel pool will be used.
% Outputs:
%   rout=An output reflectance structure. Contains the following fields:
%                   r=An Ny x Nx x Nb array
%{
rout.r=r;
rout.x=x;
rout.y=y;
rout.t=mean(t);
rout.SE=mean(SE);
rout.delta=delta;
rout.utmz=rs(1).utmz;
rout.order=rs(1).order;
rout.source=rs(1).source;
%}

% By Kristoffer Aalstad (last revised: September 2018).

indir=[indir '/']; % Pad the input directory string with an extra / (can never have too many).



%% Part 1. Read, crop and reproject.
if strcmp('Planet',sat)
    ris=r_PlanetScope(indir,aoi);
elseif strcmp('RapidEye',sat)
    ris=r_RapidEye(indir,aoi);
elseif strcmp('Sentinel-2',sat)
    ris=r_Sentinel_2(indir,aoi);
elseif strcmp('Landsat_7',sat)
    ris=r_Landsat_7(indir,aoi);
elseif strcmp('Landsat_8',sat)
    ris=r_Landsat_8(indir,aoi);
elseif strcmp('Landsat_8_SR',sat)
    ris=r_Landsat_8_SR(indir,aoi);
elseif strcmp('Landsat_7_SR',sat)
    ris=r_Landsat_7_SR(indir,aoi);
elseif strcmp('Sentinel-2_SR',sat)
    ris=r_Sentinel_2_SR(indir,aoi);
else
    error('Please specify a valid satellite constellation. Current options: Sentinel-2, Sentinel-2 SR, Landsat_8, Landsat_8_SR, Landsat_7, Landsat_7_SR, RapidEye, Planet')
end

%rout=ris;
%return


%% Part 2. Resample and mosaic.

rs=ris; clear ris;
ns=numel(rs);
delta=rs(1).delta;
nb=size(rs(1).r,3);

% Sort the scenes in order from north to south, and west to east if the
% northings are the same.
ulc=zeros(ns,1,'single');
np=0; % Number of pixels.
for s=1:ns
    ulc(s,1)=rs(s).ulx;
    ulc(s,2)=rs(s).uly;
    np=np+numel(rs(s).X(:));
end
[~,order]=sortrows(ulc,[-2 1]); % Sort north to south first then west to east.

% Define the pixel centroids for all the scenes (in the right order).
Xr=zeros(np,1,'double'); Yr=zeros(np,1,'double');
rr=zeros(np,nb,'uint8');

% Find the maximum scale factor for across all the scenes for each band.
sf=zeros(nb,1);
for s=1:ns
    for b=1:nb
        if rs(s).sf(b)>sf(b)
            sf(b)=rs(s).sf(b);
        end
    end
end

% Rescale and concatenate the scenes into a single array (called rr)
% to allow mosaicing.
here=1;
for s=1:ns
    ind=order(s);
    nps=numel(rs(ind).X(:));
    these=here:(here-1+nps);
    Xr(these)=rs(ind).X(:);
    Yr(these)=rs(ind).Y(:);
    for b=1:nb
        temp=single(rs(ind).r(:,:,b)).*rs(ind).sf(b)./sf(b);
        rr(these,b)=uint8(temp(:));
    end
    here=here+nps;
end

% Remove any missing data from the concatenated array.
% OBS! QA "bands" may have values of 0. So limit missing data to 0s in actual bands (b\in 1:6) 
these=all(rr(:,1:6)>0,2); 
Xr=Xr(these);
Yr=Yr(these);
rr=rr(these,:);


% Construct a regular grid and find the nearest neighbors.
mpx=mean(aoi.x); mpy=mean(aoi.y);
d=sqrt((Xr-mpx).^2+(Yr-mpy).^2);
mphere=find(d==min(d),1,'first');
xn=min(aoi.x)-3.*delta; xm=max(aoi.x)+3.*delta;
yn=min(aoi.y)-3.*delta; ym=max(aoi.y)+3.*delta;
xl=(Xr(mphere):-delta:xn)'; xl=flipud(xl);
xu=(Xr(mphere):delta:xm)';
x=[xl(1:end-1); xu];
yl=(Yr(mphere):-delta:yn)'; yl=flipud(yl);
yu=(Yr(mphere):delta:ym)';
y=[yl(1:end-1); yu]; y=flipud(y);

% Crop the grid to where you actually have data.
% I guess this part is problematic for S2? Fixed.
maxx=max(Xr)+0.*delta;  minx=min(Xr)-0.*delta; 
maxy=max(Yr)+0.*delta;  miny=min(Yr)-0.*delta;
x=x(x>=minx&x<=maxx);
y=y(y>=miny&y<=maxy);


[X,Y]=meshgrid(x,y);
clear xl xu yl yu d

% Find NN
disp('Finding NN');


if dopar
    c=parcluster('local');
    nc=ceil(c.NumWorkers./4); % Maybe add factor 0.5 here to avoid using all!
    p=gcp('nocreate');
    if isempty(p)
        delete(p);
        parpool('local',nc);
    else
        if ~p.Connected||p.NumWorkers~=nc
            delete(p);
            parpool('local',nc);
        end
    end
    
    nc=64; % Change number of chuncks (not workers) to avoid
                   % having too large arrays in memory on cores. Not sure
                   % if the problem is the size of Xr,Yr or Xchunck,Ychunck
    p=gcp;
    np=numel(X);
    chunck=floor(np./nc);
    NN=zeros(np,1); dNN=NN; 
    % NOTE! NN shouldn't be "single" since this can end up rounding up
    % indices when they are passed to this array.
                                                            
    % Parallel loop.
    NN_=cell(nc,1); dNN_=NN_;
    parfor n=1:nc
        if n==nc
            these=((n-1)*chunck+1):np;
        else
            these=((n-1)*chunck+1):(n*chunck);
        end
        Xchunck=X(these); Xchunck=Xchunck(:);
        Ychunck=Y(these); Ychunck=Ychunck(:);
        [NN_{n},dNN_{n}]=knnsearch([Xr Yr],[Xchunck,Ychunck]);
    end
    
    % Allocate the output of the parallel loop to NN and dNN.
    for n=1:nc
        if n==nc
            these=((n-1)*chunck+1):np;
        else
            these=((n-1)*chunck+1):(n*chunck);
        end
        NN(these)=NN_{n};
        dNN(these)=dNN_{n};
    end
    
else
    % No parallel loop, just one call to knnsearch. 
    [NN,dNN]=knnsearch([Xr(:) Yr(:)],[X(:) Y(:)]); 
end


ii=numel(y); jj=numel(x);
r=zeros(ii,jj,nb,'uint8');


for b=1:nb
    temp=rr(:,b);
    temp=temp(NN);
    % Something strange with the next line for Sentinel-2 (works fine for the
    % rest ?).
    temp(dNN>sqrt(2).*delta+1)=0; % Ignore neighbors more than one pixel away.
    r(:,:,b)=reshape(temp,ii,jj);
end



% Save the output.
t=zeros(ns,1); SE=t;
for s=1:ns
    t(s)=rs(s).t;
    SE(s)=rs(s).SE;
end

rout.r=r;
rout.sf=sf;
rout.x=x;
rout.y=y;
rout.t=mean(t);
rout.SE=mean(SE);
rout.delta=delta;
rout.utmz=rs(1).utmz;
rout.order=rs(1).order;
rout.source=rs(1).source;




%% Read functions for each satellite type.
    %% PlanetScope
    % PROBLEM reading .tif files that don't contain the string "_clip" for
    % planet. Temporary fix: ignore these files. E.g. 23-Mar-2017 directory
    % 3 is a problem.
    function ris=r_PlanetScope(indir,aoi)
        
        cont=dir(indir); %  <-- Contents of the read directory
        cont=cont(3:end);
        
        if exist([indir '1'],'dir') % Multiple scenes.
            ns=numel(cont);
        else
            ns=1;
        end
        
        sourcecell=cell(ns,1);
        for s=1:ns
            sourcecell{s}='PlanetScope';
        end
        ris=struct('source',sourcecell); % <-- Allocate the output structure.
        
        
        
        % Loop over scenes.
        for s=1:ns
            %  Define the read directory for the given scene.
            if ns==1
                readdir=indir;
            else
                readdir=[indir sprintf('%d/',s)];
            end
            
            filesare=dir(readdir);
            filesare=filesare(3:end);
            
            
            % Read in metadata for the given scene.
            mis=filesare(3).name;
            out=parseXML([readdir mis]);
            utmz=out.Children(10).Children(2).Children(2).Children(2).Children(6).Children(6).Children(1).Data;
            utmz=strsplit(utmz,'/'); utmz=utmz{2};
            utmz=strsplit(utmz,' '); utmz=utmz{4};
            SE=str2double(out.Children(6).Children(2).Children(8).Children(2).Children(8).Children.Data);
            t=out.Children(6).Children(2).Children(8).Children(2).Children(14).Children.Data;
            t=datenum([t(1:10) ' ' t(12:19)],'yyyy-mm-dd HH:MM:SS');
            
            % Read in the images.
            fis=filesare(2).name;
            
            
            if strcmp(fis(end-7:end),'clip.tif') % <-- Temporary fix, ignore problematic images.
                disp('clipped');
                [temp,sr]=geotiffread([readdir fis]); % B, G, R, NIR
                %size(temp)
                rtemp=zeros(size(temp));
                bands=[3; 2; 1; 4]; % Position of the bands in A in the R-G-B-NIR array R (see note at start of code).
                sf=zeros(4,1);
                for b=1:4
                    indis=12+2*(b-1); % For band=1:4 third child index goes 12:2:18
                    sf(b)=str2double(out.Children(10).Children(2).Children(indis).Children(10).Children(1).Data);
                    
                    rtemp(:,:,bands(b))=double(temp(:,:,b)).*sf(b);
                end
                temp=rtemp; clear rtemp;
                
                % Proc domain
                delta=sr.CellExtentInWorldX;
                x=((sr.XWorldLimits(1)+delta/2):delta:(sr.XWorldLimits(2)-delta/2))';
                y=((sr.YWorldLimits(1)+delta/2):delta:(sr.YWorldLimits(2)-delta/2))';
                y=flipud(y);
                dx=delta; dy=delta;
                ulx=x(1); uly=y(1);
                
                %  Define your aoi.
                xn=min(aoi.x)-3.*dx; xm=max(aoi.x)+3.*dx;
                yn=min(aoi.y)-3.*dy; ym=max(aoi.y)+3.*dy;
                
                % Check utm zone.
                utmz=str2double(utmz(1:2)); % UTM zone of scene
                targetutmz=str2double(aoi.utmz(1:2));
                if utmz~=targetutmz % If wrong utm zone.
                    disp('wrong utm zone');
                    
                    % First define a padded AOI in the "wrong" utm zone to speed up the
                    % interpolation/utm zone conversion.
                    cutmstruct = defaultm('utm'); %<--- The "correct" utm zone structure.
                    cutmstruct.zone = aoi.utmz;
                    cutmstruct.geoid = wgs84Ellipsoid;
                    cutmstruct = defaultm(cutmstruct);
                    [Xd,Yd]=meshgrid(aoi.x,aoi.y);
                    [LATd,LONd]=minvtran(cutmstruct,Xd(:),Yd(:));
                    
                    wutmstruct=defaultm('utm'); %<----The "wrong" utm zone structure.
                    wutmstruct.zone=sprintf('%d%s',utmz,aoi.utmz(end)); % Inherit the letter from the correct utm zone.
                    wutmstruct.geoid=wgs84Ellipsoid;
                    wutmstruct=defaultm(wutmstruct);
                    [Xdw,Ydw]=mfwdtran(wutmstruct,LATd,LONd);
                    
                    padcut=3*dx;
                    xdlim=[min(Xdw)-padcut max(Xdw)+padcut];
                    ydlim=[min(Ydw)-padcut max(Ydw)+padcut];
                    
                    
                    
                    thesex=x>=xdlim(1)&x<=xdlim(2);
                    thesey=y>=ydlim(1)&y<=ydlim(2);
                    
                    % Crop
                    % First to the domain.
                    x=x(thesex); y=y(thesey);
                    [X,Y]=meshgrid(x,y);
                    temp=temp(thesey,thesex,:);
                    notmissing=all(temp>0,3);
                    % Then crop out all the junk (missing pixels) to minimize
                    % the image size.
                    xn=min(X(notmissing)); xm=max(X(notmissing));
                    yn=min(Y(notmissing)); ym=max(Y(notmissing));
                    thesex=x>=xn&x<=xm;
                    thesey=y>=yn&y<=ym;
                    x=x(thesex); y=y(thesey);
                    [X,Y]=meshgrid(x,y);
                    temp=temp(thesey,thesex,:);
                    
                    % Now define the location of the centroids of the cut out pixels in the
                    % correct utm zone.
                    [LAT,LON]=minvtran(wutmstruct,X,Y);
                    [X,Y]=mfwdtran(cutmstruct,LAT,LON);
                    
                else
                    % Crop
                    % First to the domain.
                    thesex=x>=xn&x<=xm;
                    thesey=y>=yn&y<=ym;
                    x=x(thesex); y=y(thesey);
                    [X,Y]=meshgrid(x,y);
                    temp=temp(thesey,thesex,:);
                    
                    % Then crop out all the junk (missing pixels) to minimize
                    % the image size.
                    notmissing=all(temp>0,3);
                    xn=min(X(notmissing)); xm=max(X(notmissing));
                    yn=min(Y(notmissing)); ym=max(Y(notmissing));
                    thesex=x>=xn&x<=xm;
                    thesey=y>=yn&y<=ym;
                    x=x(thesex); y=y(thesey);
                    [X,Y]=meshgrid(x,y);
                    temp=temp(thesey,thesex,:);
                end
                
                
                
                ii=size(X,1);
                jj=size(X,2);
                ris(s).X=X;
                ris(s).Y=Y;
                ris(s).r=zeros(ii,jj,4,'uint8'); % PlanetScope is in uint16 by default.
                ris(s).sf=zeros(4,1);
                
                
                % Save the image for the given band to the output structure
                % along with a band dependent scale factor.
                % to scale from uint8 compression to the measured
                % reflectance.
                for b=1:4
                    btemp=temp(:,:,b);
                    ris(s).sf(b)=max(btemp(:))/255;
                    ris(s).r(:,:,b)=uint8(btemp./ris(s).sf(b));
                    % Might need to consider a saturation value like for S2!
                end
                
            else % Just make images you want to ignore into a single pixel at "origin".
                disp('not clipped, ignoring');
                ris(s).r=zeros(1,1,4,'uint8');
                ris(s).sf=zeros(4,1);
                ris(s).X=0;
                ris(s).Y=0;
                ulx=0; uly=0; dx=3;
            end
            
            
            ris(s).t=t;
            ris(s).order={'R'; 'G'; 'B'; 'NIR'};
            ris(s).SE=SE;
            ris(s).utmz=aoi.utmz;
            ris(s).delta=dx;
            ris(s).ulx=single(round(ulx)); % ULC identifiers for UNCROPPED scene.
            ris(s).uly=single(round(uly));
            
            
            
        end % <- End scene (s) loop.
    end % <- End read RapidEye function.

            
        
        

    %% RapidEye
    function ris=r_RapidEye(indir,aoi)
        
        cont=dir(indir); %  <-- Contents of the read directory
        cont=cont(3:end);
        
        if exist([indir '1'],'dir') % Multiple scenes.
            ns=numel(cont);
        else
            ns=1;
        end
        
        sourcecell=cell(ns,1);
        for s=1:ns
            sourcecell{s}='RapidEye';
        end
        ris=struct('source',sourcecell); % <-- Allocate the output structure.
        
        % Exo-atmospheric irradiance for RapidEye radiance to reflectance
        % conversion (B,G,R,RE,NIR)
        EAI=[1997.8; 1863.5; 1560.4; 1395.0; 1124.4];
        
        
        % Loop over scenes.
        for s=1:ns
            %  Define the read directory for the given scene.
            if ns==1
                readdir=indir;
            else
                readdir=[indir sprintf('%d/',s)];
            end
            
            filesare=dir(readdir);
            filesare=filesare(3:end);
            
            
            % Read in metadata for the given scene.
            mis=filesare(2).name;
            out=parseXML([readdir mis]);
            utmz=out.Children(10).Children(2).Children(2).Children(2).Children(8).Children(6).Children.Data;
            utmz=strsplit(utmz,'/'); utmz=utmz{2};
            utmz=strsplit(utmz,' '); utmz=utmz{4};
            SE=str2double(out.Children(6).Children(2).Children(8).Children(2).Children(8).Children.Data);
            t=out.Children(6).Children(2).Children(8).Children(2).Children(14).Children.Data;
            t=datenum([t(1:10) ' ' t(12:19)],'yyyy-mm-dd HH:MM:SS');
            
            % Read in the images.
            fis=filesare(1).name;
            [temp,sr]=geotiffread([readdir fis]);
            rtemp=zeros(size(temp,1),size(temp,2),4);
            bands=[3; 2; 1; NaN; 4]; % Position of the bands in temp in the R-G-B-NIR array R (see note at start of code).
            sf=zeros(5,1);
            year=datevec(t); year=year(1);
            doy=floor(t)-datenum(sprintf('01-Jan-%d',year));
            SunDist=1-0.01672.*cos(deg2rad(0.9856.*(doy-4))); % Estimate of earth sun distance (AU).
            for b=1:5
                if ~isnan(bands(b))
                    sf(b)=str2double(out.Children(10).Children(2).Children(12).Children(14).Children.Data);
                    rad_to_ref=pi.*SunDist.^2./(EAI(b).*cos(deg2rad(90-SE))); % See Planet  labs product specification document.
                    rtemp(:,:,bands(b))=rad_to_ref.*double(temp(:,:,b)).*sf(b);
                end
            end
            temp=rtemp; 
            clear r;
            
            % Proc domain
            delta=sr.CellExtentInWorldX;
            x=((sr.XWorldLimits(1)+delta/2):delta:(sr.XWorldLimits(2)-delta/2))';
            y=((sr.YWorldLimits(1)+delta/2):delta:(sr.YWorldLimits(2)-delta/2))';
            y=flipud(y);
            dx=delta; dy=delta;
            ulx=x(1); uly=y(1);
            
            %  Define your aoi.
            xn=min(aoi.x)-3.*dx; xm=max(aoi.x)+3.*dx;
            yn=min(aoi.y)-3.*dy; ym=max(aoi.y)+3.*dy;
            
            % Check utm zone.
            utmz=str2double(utmz(1:2)); % UTM zone of scene
            targetutmz=str2double(aoi.utmz(1:2));
            if utmz~=targetutmz % If wrong utm zone.
                disp('wrong utm zone');
                
                % First define a padded AOI in the "wrong" utm zone to speed up the
                % interpolation/utm zone conversion.
                cutmstruct = defaultm('utm'); %<--- The "correct" utm zone structure.
                cutmstruct.zone = aoi.utmz;
                cutmstruct.geoid = wgs84Ellipsoid;
                cutmstruct = defaultm(cutmstruct);
                [Xd,Yd]=meshgrid(aoi.x,aoi.y);
                [LATd,LONd]=minvtran(cutmstruct,Xd(:),Yd(:));
                
                wutmstruct=defaultm('utm'); %<----The "wrong" utm zone structure.
                wutmstruct.zone=sprintf('%d%s',utmz,aoi.utmz(end)); % Inherit the letter from the correct utm zone.
                wutmstruct.geoid=wgs84Ellipsoid;
                wutmstruct=defaultm(wutmstruct);
                [Xdw,Ydw]=mfwdtran(wutmstruct,LATd,LONd);
                
                padcut=3*dx;
                xdlim=[min(Xdw)-padcut max(Xdw)+padcut];
                ydlim=[min(Ydw)-padcut max(Ydw)+padcut];
                thesex=x>=xdlim(1)&x<=xdlim(2);
                thesey=y>=ydlim(1)&y<=ydlim(2);
                
                
                % Crop
                % First to the domain.
                x=x(thesex); y=y(thesey);
                [X,Y]=meshgrid(x,y);
                temp=temp(thesey,thesex,:);
                notmissing=all(temp>0,3);
                % Then crop out all the junk (missing pixels) to minimize
                % the image size.
                xn=min(X(notmissing)); xm=max(X(notmissing));
                yn=min(Y(notmissing)); ym=max(Y(notmissing));
                thesex=x>=xn&x<=xm;
                thesey=y>=yn&y<=ym;
                x=x(thesex); y=y(thesey);
                [X,Y]=meshgrid(x,y);
                temp=temp(thesey,thesex,:);
                
                
                % Now define the location of the centroids of the cut out pixels in the
                % correct utm zone.
                [LAT,LON]=minvtran(wutmstruct,X,Y);
                [X,Y]=mfwdtran(cutmstruct,LAT,LON);
                
            else
                % Crop
                % First to the domain.
                thesex=x>=xn&x<=xm;
                thesey=y>=yn&y<=ym;
                x=x(thesex); y=y(thesey);
                [X,Y]=meshgrid(x,y);
                temp=temp(thesey,thesex,:);
                
                % Then crop out all the junk (missing pixels) to minimize
                % the image size.
                notmissing=all(temp>0,3);
                xn=min(X(notmissing)); xm=max(X(notmissing));
                yn=min(Y(notmissing)); ym=max(Y(notmissing));
                thesex=x>=xn&x<=xm;
                thesey=y>=yn&y<=ym;
                x=x(thesex); y=y(thesey);
                [X,Y]=meshgrid(x,y);
                temp=temp(thesey,thesex,:);
                
                
            end
            ii=size(X,1);
            jj=size(X,2);
            ris(s).X=X;
            ris(s).Y=Y;
            ris(s).r=zeros(ii,jj,4,'uint8'); % Landsat 8 is in uint16 by default.
            ris(s).sf=zeros(4,1);
            
            
            % Save the image for the given band to the output structure
            % along with a band dependent scale factor.
            % to scale from uint8 compression to the measured
            % reflectance.
            for b=1:4
                btemp=temp(:,:,b);
                ris(s).sf(b)=max(btemp(:))/255;
                ris(s).r(:,:,b)=uint8(btemp./ris(s).sf(b));
                % Might need to consider a saturation value like for S2!
            end
            
            ris(s).t=t;
            ris(s).order={'R'; 'G'; 'B'; 'NIR'};
            ris(s).SE=SE;
            ris(s).utmz=aoi.utmz;
            ris(s).delta=dx;
            ris(s).ulx=single(round(ulx)); % ULC identifiers for UNCROPPED scene.
            ris(s).uly=single(round(uly));
            
            
            end % <- End scene (s) loop.  
     end % <- End read RapidEye function.
           

    %% Sentinel-2
    function ris=r_Sentinel_2(indir,aoi)
        
        cont=dir(indir); %  <-- Contents of the read directory
        cont=cont(3:end);
        
        sfs2=1e4; % Sentinel-2 scale factor.
        sats2=65535; % Saturated pixel value.
        
        if exist([indir '1'],'dir') % Multiple scenes.
            ns=numel(cont);
        else
            ns=1;
        end
        
        sourcecell=cell(ns,1);
        for s=1:ns
            sourcecell{s}='Sentinel-2';
        end
        ris=struct('source',sourcecell); % <-- Allocate the output structure.
        
        
        % Loop over scenes.
        for s=1:ns
            %  Define the read directory for the given scene.
            if ns==1
                readdir=indir;
            else
                readdir=[indir sprintf('%d/',s)];
            end
            
            filesare=dir(readdir);
            filesare=filesare(3:end);
            
            
            % Read in metadata for the given scene.
            mis=filesare(end).name;
            out=parseXML([readdir mis]);
            these=out.Children(4).Children(2);
            utmz=these.Children(2).Children.Data;
            utmz=strsplit(utmz,' ');
            utmz=utmz{end};
            nr=str2double(these.Children(6).Children(2).Children.Data);
            nc=str2double(these.Children(6).Children(4).Children.Data);
            ulx=str2double(these.Children(12).Children(2).Children.Data);
            uly=str2double(these.Children(12).Children(4).Children.Data);
            dx=str2double(these.Children(12).Children(6).Children.Data);
            dy=str2double(these.Children(12).Children(8).Children.Data);
            x=(ulx+(dx.*(0:(nc-1))))';
            y=(uly+(dy.*(0:(nr-1))))';
            t=out.Children(2).Children(8).Children.Data;
            t=datenum([t(1:10) ' ' t(12:19)]);
            these=out.Children(4).Children(4).Children(4);
            SE=90-str2double(these.Children(2).Children.Data);
            
            % Loop over bands (files in the subdirectory)
            for f=1:6
                disp(sprintf('Sentinel-2 Band %d',f));
                temp=single(imread([readdir filesare(f).name]));
                
                if f==1
                %  Define your aoi. 
                    xn=min(aoi.x)-3.*dx; xm=max(aoi.x)+3.*dx;
                    yn=min(aoi.y)-3.*dy; ym=max(aoi.y)+3.*dy;
                    
                    % Check utm zone.
                    utmz=str2double(utmz(1:2)); % UTM zone of scene
                    targetutmz=str2double(aoi.utmz(1:2));
                    if utmz~=targetutmz % If wrong utm zone.
                        disp('wrong utm zone');
                        
                        % First define a padded AOI in the "wrong" utm zone to speed up the
                        % interpolation/utm zone conversion.
                        cutmstruct = defaultm('utm'); %<--- The "correct" utm zone structure.
                        cutmstruct.zone = aoi.utmz;
                        cutmstruct.geoid = wgs84Ellipsoid;
                        cutmstruct = defaultm(cutmstruct);
                        [Xd,Yd]=meshgrid(aoi.x,aoi.y);
                        [LATd,LONd]=minvtran(cutmstruct,Xd(:),Yd(:));
                        
                        wutmstruct=defaultm('utm'); %<----The "wrong" utm zone structure.
                        wutmstruct.zone=sprintf('%d%s',utmz,aoi.utmz(end)); % Inherit the letter from the correct utm zone.
                        wutmstruct.geoid=wgs84Ellipsoid;
                        wutmstruct=defaultm(wutmstruct);
                        [Xdw,Ydw]=mfwdtran(wutmstruct,LATd,LONd);
                        
                        padcut=3*dx;
                        xdlim=[min(Xdw)-padcut max(Xdw)+padcut];
                        ydlim=[min(Ydw)-padcut max(Ydw)+padcut];
                        thesex=x>=xdlim(1)&x<=xdlim(2);
                        thesey=y>=ydlim(1)&y<=ydlim(2);
                        
                        % Now define the location of the centroids of the cut out pixels in the
                        % correct utm zone.
                        x=x(thesex); y=y(thesey);
                        [X,Y]=meshgrid(x,y);
                        [LAT,LON]=minvtran(wutmstruct,X,Y);
                        [X,Y]=mfwdtran(cutmstruct,LAT,LON);
                        
                    else
                        thesex=x>=xn&x<=xm;
                        thesey=y>=yn&y<=ym;
                        x=x(thesex); y=y(thesey);
                        [X,Y]=meshgrid(x,y);
                    end
                    ii=size(X,1);
                    jj=size(X,2);
                    ris(s).X=X;
                    ris(s).Y=Y;
                    ris(s).r=zeros(ii,jj,6,'uint8'); % Landsat 8 is in uint16 by default.
                    ris(s).sf=zeros(6,1);
            end
                
                
                % Save the image for the given band to the output structure
                % along with a band dependent scale factor.
                % to scale from uint8 compression to the measured
                % reflectance.
                if f<4
                    bis=4-f; % RGB
                else
                    bis=f; % NIR, SWIR1, SWIR2
                end
                if f>4 % % Interpolate SWIR (20 m) to 10 m grid using NN.
                    temp=imresize(temp,[nr,nc],'nearest');
                end
                temp=temp(thesey,thesex); %  Crop
                temp(temp==sats2)=max(temp(temp~=sats2)); % Set saturated values to maximum reflectance.
                scaled=double(temp)./sfs2;
                scaled(scaled<0)=0;
                ris(s).sf(bis)=max(scaled(:))./255;
                ris(s).r(:,:,bis)=uint8(scaled./ris(s).sf(bis));  
                
            end % <- End band (f) loop.
            
                ris(s).t=t;
                ris(s).order={'R'; 'G'; 'B'; 'NIR'; 'SWIR1'; 'SWIR2'};
                ris(s).SE=SE;
                ris(s).utmz=aoi.utmz;
                ris(s).delta=dx;
                ris(s).ulx=single(round(ulx)); % ULC identifiers for UNCROPPED scene.
                ris(s).uly=single(round(uly));    
                
        end % <- End scene (s) loop.  
     end % <- End read Sentinel-2 function.

    %% Landsat 7
    function ris=r_Landsat_7(indir,aoi)
        
        cont=dir(indir); %  <-- Contents of the read directory
        cont=cont(3:end);
        
        
        % Search strings  in the Landsat 7 metadata.
        searchstr=cell(16,1); searchval=zeros(size(searchstr));
        for str=1:numel(searchstr)
            if str==1
                searchstr{str}='DATE_ACQUIRED';
            elseif str==2
                searchstr{str}='SCENE_CENTER_TIME';
            elseif str==3
                searchstr{str}='SUN_ELEVATION';
            elseif str<=9
                b=str-3; % Bands 1-5&7
                b=b.*(b~=6)+7.*(b==6);
                searchstr{str}=sprintf('REFLECTANCE_MULT_BAND_%d',b);
            elseif str<=15
                b=str-9; % Bands 1-5&7
                b=b.*(b~=6)+7.*(b==6);
                searchstr{str}=sprintf('REFLECTANCE_ADD_BAND_%d',b);
            else
                searchstr{str}='UTM_ZONE';
            end
        end
       
        
         
        if exist([indir '1'],'dir') % Multiple scenes.
            ns=numel(cont);
        else
            ns=1;
        end
        
        sourcecell=cell(ns,1);
        for s=1:ns
            sourcecell{s}='Landsat 7';
        end
        ris=struct('source',sourcecell); % <-- Allocate the output structure.
        
        
        % Loop over scenes.
        for s=1:ns
            %  Define the read directory for the given scene.
            if ns==1
                readdir=indir;
            else
                readdir=[indir sprintf('%d/',s)];
            end
            
            filesare=dir(readdir);
            filesare=filesare(3:end);
            
            
            % Read in metadata for the given scene.
            mis=[readdir filesare(end).name];
            fid=fopen(mis);
            t=[];
            
            str=1;
            while str<=numel(searchstr)
                line=fgetl(fid);
                lines=strsplit(line,'=');
                if any(strcmp(searchstr{str},strtrim(lines{1})))
                    if str<=2
                        t=[t lines{end}];
                        if str==2
                            t=datenum([t(1:11) ' ' t(14:21)],'yyyy-mm-dd HH:MM:SS');
                        end
                    else
                        searchval(str)=str2double(lines{end});
                    end
                    str=str+1;
                end
            end
            
          
            % Loop over bands
            for f=1:6
                disp(sprintf('Landsat 7 Band %d',f));
                [temp,sr]=geotiffread([readdir filesare(f).name]);
                % Retrieve georeferencing information for first band (same for
                % all).
                if f==1
                    dx=sr.SampleSpacingInWorldX;
                    dy=sr.SampleSpacingInWorldY;
                    nr=sr.RasterSize(1);
                    nc=sr.RasterSize(2);
                    ulx=sr.XWorldLimits(1);
                    uly=sr.YWorldLimits(1);
                    x=single((ulx+(dx.*(0:(nc-1))))');
                    y=single(flipud((uly+(dy.*(0:(nr-1))))'));
                    
                    
                    %  Define your aoi.
                    xn=min(aoi.x)-3.*dx; xm=max(aoi.x)+3.*dx;
                    yn=min(aoi.y)-3.*dy; ym=max(aoi.y)+3.*dy;
                    
                    utmz=searchval(end); % UTM zone of scene
                    targetutmz=str2double(aoi.utmz(1:2));
                    if utmz~=targetutmz % If wrong utm zone.
                        disp('wrong utm zone');
                        
                        % First define a padded AOI in the "wrong" utm zone to speed up the
                        % interpolation/utm zone conversion.
                        cutmstruct = defaultm('utm'); %<--- The "correct" utm zone structure.
                        cutmstruct.zone = aoi.utmz;
                        cutmstruct.geoid = wgs84Ellipsoid;
                        cutmstruct = defaultm(cutmstruct);
                        [Xd,Yd]=meshgrid(aoi.x,aoi.y);
                        [LATd,LONd]=minvtran(cutmstruct,Xd(:),Yd(:));
                        
                        wutmstruct=defaultm('utm'); %<----The "wrong" utm zone structure.
                        wutmstruct.zone=sprintf('%d%s',utmz,aoi.utmz(end)); % Inherit the letter from the correct utm zone.
                        wutmstruct.geoid=wgs84Ellipsoid;
                        wutmstruct=defaultm(wutmstruct);
                        [Xdw,Ydw]=mfwdtran(wutmstruct,LATd,LONd);
                        
                        padcut=3*dx;
                        xdlim=[min(Xdw)-padcut max(Xdw)+padcut];
                        ydlim=[min(Ydw)-padcut max(Ydw)+padcut];
                        thesex=x>=xdlim(1)&x<=xdlim(2);
                        thesey=y>=ydlim(1)&y<=ydlim(2);
                        
                        % Now define the location of the centroids of the cut out pixels in the
                        % correct utm zone.
                        x=x(thesex); y=y(thesey);
                        [X,Y]=meshgrid(x,y);
                        [LAT,LON]=minvtran(wutmstruct,X,Y);
                        [X,Y]=mfwdtran(cutmstruct,LAT,LON);
                    else
                        thesex=x>=xn&x<=xm;
                        thesey=y>=yn&y<=ym;
                        x=x(thesex); y=y(thesey);
                        [X,Y]=meshgrid(x,y);
                    end
                    nr=size(X,1); 
                    nc=size(X,2);
                    ris(s).X=X; 
                    ris(s).Y=Y; 
                    ris(s).r=zeros(nr,nc,6,'uint8'); % Landsat 8 is in uint16 by default.
                    ris(s).sf=zeros(6,1);
                    
                end
                
                % Save the image for the given band to the output structure
                % along with a band dependent scale factor.
                % to scale from uint16 compression to the measured
                % reflectance.
                if f<4
                    bis=4-f; % RGB
                else
                    bis=f; % NIR, SWIR1, SWIR2
                end
                scaled=(searchval(3+f).*double(temp(thesey,thesex))+searchval(9+f))./sin(deg2rad(searchval(3)));
                scaled(scaled<0)=0;    
                ris(s).sf(bis)=max(scaled(:))./255;
                ris(s).r(:,:,bis)=uint8(scaled./ris(s).sf(bis));
              
            end % <- End band (f) loop.
            
                ris(s).t=t;
                ris(s).order={'R'; 'G'; 'B'; 'NIR'; 'SWIR1'; 'SWIR2'};
                ris(s).SE=searchval(3);
                ris(s).utmz=aoi.utmz;
                ris(s).delta=dx;
                ris(s).ulx=single(round(ulx)); % ULC identifiers for UNCROPPED scene.
                ris(s).uly=single(round(uly));
            
            
        end % <- End scene (s) loop.
        
    end % <- End read Landsat_7 function.

        
        

    %% Landsat 8
    function ris=r_Landsat_8(indir,aoi)
        
        cont=dir(indir); %  <-- Contents of the read directory
        cont=cont(3:end);
        
        % Search strings  in the Landsat 8 metadata.
        searchstr=cell(16,1); searchval=zeros(size(searchstr));
        for str=1:numel(searchstr)
            if str==1
                searchstr{str}='DATE_ACQUIRED';
            elseif str==2
                searchstr{str}='SCENE_CENTER_TIME';
            elseif str==3
                searchstr{str}='SUN_ELEVATION';
            elseif str<=9 % Bands 2,3,4,5,6,7
                searchstr{str}=sprintf('REFLECTANCE_MULT_BAND_%d',str-2);
            elseif str<=15 %  Bands 2,3,4,5,6,7
                searchstr{str}=sprintf('REFLECTANCE_ADD_BAND_%d',str-8);
            else
                searchstr{str}='UTM_ZONE';
            end
        end
        
        
        if exist([indir '1'],'dir') % Multiple scenes.
            ns=numel(cont);
        else
            ns=1;
        end
        
        sourcecell=cell(ns,1);
        for s=1:ns
            sourcecell{s}='Landsat 8';
        end
        ris=struct('source',sourcecell); % <-- Allocate the output structure.
        
        
        % Loop over scenes.
        for s=1:ns
            %  Define the read directory for the given scene.
            if ns==1
                readdir=indir;
            else
                readdir=[indir sprintf('%d/',s)];
            end
            
            filesare=dir(readdir);
            filesare=filesare(3:end);
            
            % Read in metadata for the given scene.
            mis=[readdir filesare(end).name];
            fid=fopen(mis);
            t=[];
            str=1;
            while str<=numel(searchstr)
                line=fgetl(fid);
                lines=strsplit(line,'=');
                if any(strcmp(searchstr{str},strtrim(lines{1})))
                    if str<=2
                        t=[t lines{end}];
                        if str==2
                            t=datenum([t(1:11) ' ' t(14:21)],'yyyy-mm-dd HH:MM:SS');
                        end
                    else
                        searchval(str)=str2double(lines{end});
                    end
                    str=str+1;
                end
            end
            
            % Loop over bands
            for f=1:6
                disp(sprintf('Landsat 8 Band %d',f));
                [temp,sr]=geotiffread([readdir filesare(f).name]);
                % Retrieve georeferencing information for first band (same for
                % all).
                if f==1
                    dx=sr.SampleSpacingInWorldX;
                    dy=sr.SampleSpacingInWorldY;
                    nr=sr.RasterSize(1);
                    nc=sr.RasterSize(2);
                    ulx=sr.XWorldLimits(1);
                    uly=sr.YWorldLimits(1);
                    x=double((ulx+(dx.*(0:(nc-1))))');
                    y=double(flipud((uly+(dy.*(0:(nr-1))))'));
                    
                    
                    %  Define your aoi.
                    xn=min(aoi.x)-3.*dx; xm=max(aoi.x)+3.*dx;
                    yn=min(aoi.y)-3.*dy; ym=max(aoi.y)+3.*dy;
                    
                    utmz=searchval(end); % UTM zone of scene
                    targetutmz=str2double(aoi.utmz(1:2));
                    if utmz~=targetutmz % If wrong utm zone.
                        disp('wrong utm zone');
                        
                        % First define a padded AOI in the "wrong" utm zone to speed up the
                        % interpolation/utm zone conversion.
                        cutmstruct = defaultm('utm'); %<--- The "correct" utm zone structure.
                        cutmstruct.zone = aoi.utmz;
                        cutmstruct.geoid = wgs84Ellipsoid;
                        cutmstruct = defaultm(cutmstruct);
                        [Xd,Yd]=meshgrid(aoi.x,aoi.y);
                        [LATd,LONd]=minvtran(cutmstruct,Xd(:),Yd(:));
                        
                        wutmstruct=defaultm('utm'); %<----The "wrong" utm zone structure.
                        wutmstruct.zone=sprintf('%d%s',utmz,aoi.utmz(end)); % Inherit the letter from the correct utm zone.
                        wutmstruct.geoid=wgs84Ellipsoid;
                        wutmstruct=defaultm(wutmstruct);
                        [Xdw,Ydw]=mfwdtran(wutmstruct,LATd,LONd);
                        
                        padcut=3*dx;
                        xdlim=[min(Xdw)-padcut max(Xdw)+padcut];
                        ydlim=[min(Ydw)-padcut max(Ydw)+padcut];
                        thesex=x>=xdlim(1)&x<=xdlim(2);
                        thesey=y>=ydlim(1)&y<=ydlim(2);
                        
                        % Now define the location of the centroids of the cut out pixels in the
                        % correct utm zone.
                        x=x(thesex); y=y(thesey);
                        [X,Y]=meshgrid(x,y);
                        [LAT,LON]=minvtran(wutmstruct,X,Y);
                        [X,Y]=mfwdtran(cutmstruct,LAT,LON);
                    else
                        thesex=x>=xn&x<=xm;
                        thesey=y>=yn&y<=ym;
                        x=x(thesex); y=y(thesey);
                        [X,Y]=meshgrid(x,y);
                    end
                    nr=size(X,1); 
                    nc=size(X,2);
                    ris(s).X=X; 
                    ris(s).Y=Y; 
                    ris(s).r=zeros(nr,nc,6,'uint8'); % Landsat 8 is in uint16 by default.
                    ris(s).sf=zeros(6,1);
                    
                end
                
                % Save the image for the given band to the output structure
                % along with a band dependent scale factor.
                % to scale from uint16 compression to the measured
                % reflectance.
                if f<4
                    bis=4-f; % RGB
                else
                    bis=f; % NIR, SWIR1, SWIR2
                end
                scaled=(searchval(3+f).*double(temp(thesey,thesex))+searchval(9+f))./sin(deg2rad(searchval(3)));
                scaled(scaled<0)=0;    
                ris(s).sf(bis)=max(scaled(:))./255;
                ris(s).r(:,:,bis)=uint8(scaled./ris(s).sf(bis));
              
            end % <- End band (f) loop.
            
                ris(s).t=t;
                ris(s).order={'R'; 'G'; 'B'; 'NIR'; 'SWIR1'; 'SWIR2'};
                ris(s).SE=searchval(3);
                ris(s).utmz=aoi.utmz;
                ris(s).delta=dx;
                ris(s).ulx=single(round(ulx)); % ULC identifiers for UNCROPPED scene.
                ris(s).uly=single(round(uly));
            
            
        end % <- End scene (s) loop.
        
    end % <- End read Landsat_8 function.


    %% Landsat 8 surface reflectance
    function ris=r_Landsat_8_SR(indir,aoi)
        
        cont=dir(indir); %  <-- Contents of the read directory
        cont=cont(3:end);
        Nb=8;
        
        % Search strings  in the Landsat 8 metadata.
        searchstr=cell(5,1); searchval=zeros(size(searchstr));
        searchstr{1}='LANDSAT_SCENE_ID';
        searchstr{2}='DATE_ACQUIRED';
        searchstr{3}='SCENE_CENTER_TIME';
        searchstr{4}='SUN_ELEVATION';
        searchstr{5}='UTM_ZONE';
        
        % 'Band' file names [the last 2 'bands' are not bands but QA and
        % saturation layers].
        bandstr=cell(Nb,1);
        bandstr{1}='band2.tif';
        bandstr{2}='band3.tif';
        bandstr{3}='band4.tif';
        bandstr{4}='band5.tif';
        bandstr{5}='band6.tif';
        bandstr{6}='band7.tif';
        bandstr{7}='pixel';
        bandstr{8}='radsat';
        
        if exist([indir '1'],'dir') % Multiple scenes.
            ns=numel(cont);
        else
            ns=1;
        end
        
        sourcecell=cell(ns,1);
        for s=1:ns
            sourcecell{s}='Landsat-8 SR';
        end
        ris=struct('source',sourcecell); % <-- Allocate the output structure.
        
        
        % Loop over scenes.
        for s=1:ns
            %  Define the read directory for the given scene.
            if ns==1
                readdir=indir;
            else
                readdir=[indir sprintf('%d/',s)];
            end
            
            filesare=dir(readdir);
            filesare=filesare(3:end);
            
            % Read in metadata for the given scene.
            nf=numel(filesare); ok=0; f=1;
            while ~ok
                fis=filesare(f).name;
                snis=strsplit(fis,'_');
                if strcmp(snis(end),'MTL.txt')
                    ok=1;
                elseif f==nf
                    error('Couldn''nt find MTL file');
                else
                    f=f+1;
                end
            end
            mis=[readdir fis];
            fid=fopen(mis);
            t=[];
            str=1;
            while str<=numel(searchstr)
                
                line=fgetl(fid);
                lines=strsplit(line,'=');
                if any(strcmp(searchstr{str},strtrim(lines{1})))
                    if str==1% Scene ID
                        sceneid=strtrim(lines{1});
                    elseif str<=3 % Time stamp
                        t=[t lines{end}];
                        if str==3
                            t=datenum([t(1:11) ' ' t(14:21)],'yyyy-mm-dd HH:MM:SS');
                        end
                    else
                        searchval(str)=str2double(lines{end});
                    end
                    str=str+1;
                end
            end
            
            % Loop over bands
            for f=1:8
                fprintf('\n Landsat 8 SR Band %d \n',f);
                % Find the relevant file.
                nf=numel(filesare); ok=0; k=1;
                while ~ok
                    fis=filesare(k).name;
                    snis=strsplit(fis,'_');
                    if strcmp(snis(end),bandstr{f})||strcmp(snis(end-1),bandstr{f})
                        ok=1;
                    elseif k==nf
                        error('Couldn''t find band file %d',f);
                    else
                        k=k+1;
                    end
                end
                
                
                % Read in the file.
                [temp,sr]=geotiffread([readdir filesare(k).name]); % k is the actual file index.
                
                % Retrieve georeferencing information for first band (same for
                % all).
                if f==1
                    dx=sr.CellExtentInWorldX;
                    dy=sr.CellExtentInWorldY;
                    nr=sr.RasterSize(1);
                    nc=sr.RasterSize(2);
                    ulx=sr.XWorldLimits(1);
                    uly=sr.YWorldLimits(1);
                    x=double((ulx+(dx.*(0:(nc-1))))');
                    y=double(flipud((uly+(dy.*(0:(nr-1))))'));
                    
                    
                    %  Define your aoi.
                    xn=min(aoi.x)-3.*dx; xm=max(aoi.x)+3.*dx;
                    yn=min(aoi.y)-3.*dy; ym=max(aoi.y)+3.*dy;
                    
                    utmz=searchval(end); % UTM zone of scene
                    targetutmz=str2double(aoi.utmz(1:2));
                    if utmz~=targetutmz % If wrong utm zone.
                        disp('wrong utm zone');
                        
                        % First define a padded AOI in the "wrong" utm zone to speed up the
                        % interpolation/utm zone conversion.
                        cutmstruct = defaultm('utm'); %<--- The "correct" utm zone structure.
                        cutmstruct.zone = aoi.utmz;
                        cutmstruct.geoid = wgs84Ellipsoid;
                        cutmstruct = defaultm(cutmstruct);
                        [Xd,Yd]=meshgrid(aoi.x,aoi.y);
                        [LATd,LONd]=minvtran(cutmstruct,Xd(:),Yd(:));
                        
                        wutmstruct=defaultm('utm'); %<----The "wrong" utm zone structure.
                        wutmstruct.zone=sprintf('%d%s',utmz,aoi.utmz(end)); % Inherit the letter from the correct utm zone.
                        wutmstruct.geoid=wgs84Ellipsoid;
                        wutmstruct=defaultm(wutmstruct);
                        [Xdw,Ydw]=mfwdtran(wutmstruct,LATd,LONd);
                        
                        padcut=3*dx;
                        xdlim=[min(Xdw)-padcut max(Xdw)+padcut];
                        ydlim=[min(Ydw)-padcut max(Ydw)+padcut];
                        thesex=x>=xdlim(1)&x<=xdlim(2);
                        thesey=y>=ydlim(1)&y<=ydlim(2);
                        
                        % Now define the location of the centroids of the cut out pixels in the
                        % correct utm zone.
                        x=x(thesex); y=y(thesey);
                        [X,Y]=meshgrid(x,y);
                        [LAT,LON]=minvtran(wutmstruct,X,Y);
                        [X,Y]=mfwdtran(cutmstruct,LAT,LON);
                    else
                        thesex=x>=xn&x<=xm;
                        thesey=y>=yn&y<=ym;
                        x=x(thesex); y=y(thesey);
                        [X,Y]=meshgrid(x,y);
                    end
                    nr=size(X,1);
                    nc=size(X,2);
                    ris(s).X=X;
                    ris(s).Y=Y;
                    ris(s).r=zeros(nr,nc,Nb,'uint8'); % Landsat 8 is in uint16 by default.
                    ris(s).sf=ones(Nb,1);
                    
                end
                
                
                % Save the image for the given band to the output structure
                % along with a band dependent scale factor.
                if fis<4
                    bis=4-f; % RGB format (for simple visualization).
                else
                    bis=f;
                end
                if bis<7 % An acutal band, rescale from uint16 to uint8.
                    
                    % Fill values will have ALL bands equal to 0.01 for
                    % surface reflectance.
                    scaled=double(temp(thesey,thesex)).*1e-4;
                    scaled(scaled<0)=0.01; 
                    scaled(scaled>1)=1; % Saturation flags stored in own layer.
                    ris(s).sf(bis)=max(scaled(:))./255;
                    ris(s).r(:,:,bis)=uint8(scaled./ris(s).sf(bis));
                elseif  bis==7 % pixel QA "band"
                    %% Simplify the pixel QA:
                    % Simplify qa bits. Not interested in cirrus or fill values, so can
                    % create an unsigned 8 bit integer. Note, matlab starts counting (bits too)
                    % at 1 not 0.
                    qa=temp(thesey,thesex);
                    qa=bitset(qa,9,0); qa=bitset(qa,10,0); % Remove cirrus bits.
                    tocc=bitget(qa,11,'uint16'); % Terrain occlusion bit.
                    qa=bitset(qa,1,0); % Remove the "fill" flag in bit 1 (i.e. bit 0).
                    qa(tocc==1)=bitset(qa(tocc==1),1,1); % Set the first bit to the terrain occlusion bit.
                    ris(s).r(:,:,bis)=qa;
                elseif bis==8 % Pixel saturation "band"
                    %% Simplify the saturation QA, all bands > 7 (bits>7 (i.e. 8 in matlab)) are irrelevant.
                    sat=temp(thesey,thesex);
                    for bit=9:12
                        sat=bitset(sat,bit,0);
                    end
                    % Bits 0 through 7 (i.e. 1 through 8 in matlab) represent: data fill
                    % flag (bit 0) and band saturation flags (bits 1-7).
                    ris(s).r(:,:,bis)=sat;
                end % <- End band (f) loop.
                
                ris(s).t=t;
                ris(s).order={'R'; 'G'; 'B'; 'NIR'; 'SWIR1'; 'SWIR2'; 'pixel_QA'; 'saturation_QA'};
                ris(s).SE=searchval(5);
                ris(s).utmz=aoi.utmz;
                ris(s).delta=dx;
                ris(s).ulx=single(round(ulx)); % ULC identifiers for UNCROPPED scene.
                ris(s).uly=single(round(uly));
                
                
                
            end  %<- End "band" (f) loop.
            
        end % <- End scene (s) loop.
            
        end % <- End read Landsat_8 SR function.
    
    
    %% Landsat 7 surface reflectance
    function ris=r_Landsat_7_SR(indir,aoi)
        
        cont=dir(indir); %  <-- Contents of the read directory
        cont=cont(3:end);
        Nb=8;
        
        % Search strings  in the Landsat 8 metadata.
        searchstr=cell(5,1); searchval=zeros(size(searchstr));
        searchstr{1}='LANDSAT_SCENE_ID';
        searchstr{2}='DATE_ACQUIRED';
        searchstr{3}='SCENE_CENTER_TIME';
        searchstr{4}='SUN_ELEVATION';
        searchstr{5}='UTM_ZONE';
        
        % 'Band' file names [the last 2 'bands' are not bands but QA and
        % saturation layers].
        bandstr=cell(Nb,1);
        bandstr{1}='band1.tif';
        bandstr{2}='band2.tif';
        bandstr{3}='band3.tif';
        bandstr{4}='band4.tif';
        bandstr{5}='band5.tif';
        bandstr{6}='band7.tif';
        bandstr{7}='pixel';
        bandstr{8}='radsat';
        
        if exist([indir '1'],'dir') % Multiple scenes.
            ns=numel(cont);
        else
            ns=1;
        end
        
        sourcecell=cell(ns,1);
        for s=1:ns
            sourcecell{s}='Landsat-7 SR';
        end
        ris=struct('source',sourcecell); % <-- Allocate the output structure.
        
        
        % Loop over scenes.
        for s=1:ns
            %  Define the read directory for the given scene.
            if ns==1
                readdir=indir;
            else
                readdir=[indir sprintf('%d/',s)];
            end
            
            filesare=dir(readdir);
            filesare=filesare(3:end);
            
            % Read in metadata for the given scene.
            nf=numel(filesare); ok=0; f=1;
            while ~ok
                fis=filesare(f).name;
                snis=strsplit(fis,'_');
                if strcmp(snis(end),'MTL.txt')
                    ok=1;
                elseif f==nf
                    error('Couldn''nt find MTL file');
                else
                    f=f+1;
                end
            end
            mis=[readdir fis];
            fid=fopen(mis);
            t=[];
            str=1;
            while str<=numel(searchstr)
                
                line=fgetl(fid);
                lines=strsplit(line,'=');
                if any(strcmp(searchstr{str},strtrim(lines{1})))
                    if str==1% Scene ID
                        sceneid=strtrim(lines{1});
                    elseif str<=3 % Time stamp
                        t=[t lines{end}];
                        if str==3
                            t=datenum([t(1:11) ' ' t(14:21)],'yyyy-mm-dd HH:MM:SS');
                        end
                    else
                        searchval(str)=str2double(lines{end});
                    end
                    str=str+1;
                end
            end
            
            % Loop over bands
            for f=1:8
                fprintf('\n Landsat 7 SR Band %d \n',f);
                % Find the relevant file.
                nf=numel(filesare); ok=0; k=1;
                while ~ok
                    fis=filesare(k).name;
                    snis=strsplit(fis,'_');
                    if strcmp(snis(end),bandstr{f})||strcmp(snis(end-1),bandstr{f})
                        ok=1;
                    elseif k==nf
                        error('Couldn''t find band file %d',f);
                    else
                        k=k+1;
                    end
                end
                
                
                % Read in the file.
                [temp,sr]=geotiffread([readdir filesare(k).name]); % k is the actual file index.
                
                % Retrieve georeferencing information for first band (same for
                % all).
                if f==1
                    dx=sr.CellExtentInWorldX;
                    dy=sr.CellExtentInWorldY;
                    nr=sr.RasterSize(1);
                    nc=sr.RasterSize(2);
                    ulx=sr.XWorldLimits(1);
                    uly=sr.YWorldLimits(1);
                    x=double((ulx+(dx.*(0:(nc-1))))');
                    y=double(flipud((uly+(dy.*(0:(nr-1))))'));
                    
                    
                    %  Define your aoi.
                    xn=min(aoi.x)-3.*dx; xm=max(aoi.x)+3.*dx;
                    yn=min(aoi.y)-3.*dy; ym=max(aoi.y)+3.*dy;
                    
                    utmz=searchval(end); % UTM zone of scene
                    targetutmz=str2double(aoi.utmz(1:2));
                    if utmz~=targetutmz % If wrong utm zone.
                        disp('wrong utm zone');
                        
                        % First define a padded AOI in the "wrong" utm zone to speed up the
                        % interpolation/utm zone conversion.
                        cutmstruct = defaultm('utm'); %<--- The "correct" utm zone structure.
                        cutmstruct.zone = aoi.utmz;
                        cutmstruct.geoid = wgs84Ellipsoid;
                        cutmstruct = defaultm(cutmstruct);
                        [Xd,Yd]=meshgrid(aoi.x,aoi.y);
                        [LATd,LONd]=minvtran(cutmstruct,Xd(:),Yd(:));
                        
                        wutmstruct=defaultm('utm'); %<----The "wrong" utm zone structure.
                        wutmstruct.zone=sprintf('%d%s',utmz,aoi.utmz(end)); % Inherit the letter from the correct utm zone.
                        wutmstruct.geoid=wgs84Ellipsoid;
                        wutmstruct=defaultm(wutmstruct);
                        [Xdw,Ydw]=mfwdtran(wutmstruct,LATd,LONd);
                        
                        padcut=3*dx;
                        xdlim=[min(Xdw)-padcut max(Xdw)+padcut];
                        ydlim=[min(Ydw)-padcut max(Ydw)+padcut];
                        thesex=x>=xdlim(1)&x<=xdlim(2);
                        thesey=y>=ydlim(1)&y<=ydlim(2);
                        
                        % Now define the location of the centroids of the cut out pixels in the
                        % correct utm zone.
                        x=x(thesex); y=y(thesey);
                        [X,Y]=meshgrid(x,y);
                        [LAT,LON]=minvtran(wutmstruct,X,Y);
                        [X,Y]=mfwdtran(cutmstruct,LAT,LON);
                    else
                        thesex=x>=xn&x<=xm;
                        thesey=y>=yn&y<=ym;
                        x=x(thesex); y=y(thesey);
                        [X,Y]=meshgrid(x,y);
                    end
                    nr=size(X,1);
                    nc=size(X,2);
                    ris(s).X=X;
                    ris(s).Y=Y;
                    ris(s).r=zeros(nr,nc,Nb,'uint8'); % Landsat 8 is in uint16 by default.
                    ris(s).sf=ones(Nb,1);
                    
                end
                
                
                % Save the image for the given band to the output structure
                % along with a band dependent scale factor.
                if fis<4
                    bis=4-f; % RGB format (for simple visualization).
                else
                    bis=f;
                end
                if bis<7 % An acutal band, rescale from uint16 to uint8.
                    
                    % Fill values will have ALL bands equal to 0.01 for
                    % surface reflectance.
                    scaled=double(temp(thesey,thesex)).*1e-4;
                    scaled(scaled<0)=0.01;
                    scaled(scaled>1)=1; % Saturation flags stored in own layer.
                    ris(s).sf(bis)=max(scaled(:))./255;
                    ris(s).r(:,:,bis)=uint8(scaled./ris(s).sf(bis));
                elseif  bis==7 % pixel QA "band"
                    %% Simplify the pixel QA:
                    qa=uint8(temp(thesey,thesex)); % All bits higher than 7 (8 in matlab) are unused.
                    ris(s).r(:,:,bis)=qa;
                elseif bis==8 % Pixel saturation "band"
                    sat=uint8(temp(thesey,thesex)); % All bits higher than 7 (8 in matlab) are unused.
                    ris(s).r(:,:,bis)=sat;
                end % <- End band (f) loop.
                
                ris(s).t=t;
                ris(s).order={'R'; 'G'; 'B'; 'NIR'; 'SWIR1'; 'SWIR2'; 'pixel_QA'; 'saturation_QA'};
                ris(s).SE=searchval(5);
                ris(s).utmz=aoi.utmz;
                ris(s).delta=dx;
                ris(s).ulx=single(round(ulx)); % ULC identifiers for UNCROPPED scene.
                ris(s).uly=single(round(uly));
                
                
                
            end  %<- End "band" (f) loop.
            
        end % <- End scene (s) loop.
        
    end % <- End read Landsat_8 SR function.
    
    
    
    
    %% Sentinel-2 SR
    function ris=r_Sentinel_2_SR(indir,aoi)
        
        cont=dir(indir); %  <-- Contents of the read directory
        cont=cont(3:end);
        
        sfs2=1e4; % Sentinel-2 scale factor.
        sats2=65535; % Saturated pixel value.
        
        if exist([indir '1'],'dir') % Multiple scenes.
            ns=numel(cont);
        else
            ns=1;
        end
        
        sourcecell=cell(ns,1);
        for s=1:ns
            sourcecell{s}='Sentinel-2';
        end
        ris=struct('source',sourcecell); % <-- Allocate the output structure.
        
        
        % Loop over scenes.
        for s=1:ns
            %  Define the read directory for the given scene.
            if ns==1
                readdir=indir;
            else
                readdir=[indir sprintf('%d/',s)];
            end
            
            filesare=dir(readdir);
            filesare=filesare(3:end);
            
            
            % Read in metadata for the given scene.
            % Changes position depending on the type of Sen2Cor run.
            nf=numel(filesare);
            for f=1:nf
                fis=filesare(f).name;
                fis=strsplit(fis,'.');
                if strcmp(fis{end},'xml')
                    mis=filesare(f).name;
                end
            end
            out=parseXML([readdir mis]);
            sourceis=out.Children(2).Children(2).Children.Data;
            sourceis=sourceis(1:3)
            ris(s).source=sourceis;
            these=out.Children(4).Children(2);
            utmz=these.Children(2).Children.Data;
            utmz=strsplit(utmz,' ');
            utmz=utmz{end};
            nr=str2double(these.Children(6).Children(2).Children.Data);
            nc=str2double(these.Children(6).Children(4).Children.Data);
            ulx=str2double(these.Children(12).Children(2).Children.Data);
            uly=str2double(these.Children(12).Children(4).Children.Data);
            dx=str2double(these.Children(12).Children(6).Children.Data);
            dy=str2double(these.Children(12).Children(8).Children.Data);
            x=(ulx+(dx.*(0:(nc-1))))';
            y=(uly+(dy.*(0:(nr-1))))';
            
            % Sometimes the position of sensing time changes in the xml file...
            if strcmp(out.Children(2).Children(8).Name,'SENSING_TIME')
                t=out.Children(2).Children(8).Children.Data;
            else
                t=out.Children(2).Children(10).Children.Data;
            end
            t
            t=datenum([t(1:10) ' ' t(12:19)]);
            these=out.Children(4).Children(4).Children(4);
            SE=90-str2double(these.Children(2).Children.Data);
            
            % Read in the scene classification file. 
            % Changes position depending on the type of Sen2Cor run.
            nf=numel(filesare);
            for f=1:nf
                fis=filesare(f).name;
                fis=strsplit(fis,'_');
                if any(strcmp(fis,'SCL'))
                    sclis=filesare(f).name;
                end
            end
            scl=imread([readdir sclis]);
            % Resample to 10 meters.
            scl=imresize(scl,[nr,nc],'nearest');
            
            % Loop over bands (files in the subdirectory)
            seeks={'B02'; 'B03'; 'B04'; 'B08'; 'B11'; 'B12'};
            
            
            for f=1:6
                disp(sprintf('Sentinel-2 SR Band %d',f));
                
                % Find the right file...
                for k=1:numel(filesare)
                    fis=filesare(k).name;
                    fis=strsplit(fis,'_');
                    if any(strcmp(fis,seeks{f}))
                        kis=k;
                    end
                end
                fis=filesare(kis).name;
                
                
                temp=single(imread([readdir fis]));
                
                if f==1
                %  Define your aoi. 
                    xn=min(aoi.x)-3.*dx; xm=max(aoi.x)+3.*dx;
                    yn=min(aoi.y)-3.*dy; ym=max(aoi.y)+3.*dy;
                    
                    % Check utm zone.
                    utmz=str2double(utmz(1:2)); % UTM zone of scene
                    targetutmz=str2double(aoi.utmz(1:2));
                    if utmz~=targetutmz % If wrong utm zone.
                        disp('wrong utm zone');
                        
                        % First define a padded AOI in the "wrong" utm zone to speed up the
                        % interpolation/utm zone conversion.
                        cutmstruct = defaultm('utm'); %<--- The "correct" utm zone structure.
                        cutmstruct.zone = aoi.utmz;
                        cutmstruct.geoid = wgs84Ellipsoid;
                        cutmstruct = defaultm(cutmstruct);
                        [Xd,Yd]=meshgrid(aoi.x,aoi.y);
                        [LATd,LONd]=minvtran(cutmstruct,Xd(:),Yd(:));
                        
                        wutmstruct=defaultm('utm'); %<----The "wrong" utm zone structure.
                        wutmstruct.zone=sprintf('%d%s',utmz,aoi.utmz(end)); % Inherit the letter from the correct utm zone.
                        wutmstruct.geoid=wgs84Ellipsoid;
                        wutmstruct=defaultm(wutmstruct);
                        [Xdw,Ydw]=mfwdtran(wutmstruct,LATd,LONd);
                        
                        padcut=3*dx;
                        xdlim=[min(Xdw)-padcut max(Xdw)+padcut];
                        ydlim=[min(Ydw)-padcut max(Ydw)+padcut];
                        thesex=x>=xdlim(1)&x<=xdlim(2);
                        thesey=y>=ydlim(1)&y<=ydlim(2);
                        
                        % Now define the location of the centroids of the cut out pixels in the
                        % correct utm zone.
                        x=x(thesex); y=y(thesey);
                        [X,Y]=meshgrid(x,y);
                        [LAT,LON]=minvtran(wutmstruct,X,Y);
                        [X,Y]=mfwdtran(cutmstruct,LAT,LON);
                        
                    else
                        thesex=x>=xn&x<=xm;
                        thesey=y>=yn&y<=ym;
                        x=x(thesex); y=y(thesey);
                        [X,Y]=meshgrid(x,y);
                    end
                    ii=size(X,1);
                    jj=size(X,2);
                    ris(s).X=X;
                    ris(s).Y=Y;
                    ris(s).r=zeros(ii,jj,7,'uint8'); % Add extra band for scene classification.
                    ris(s).sf=ones(7,1);
            end
                
                
                % Save the image for the given band to the output structure
                % along with a band dependent scale factor.
                % to scale from uint8 compression to the measured
                % reflectance.
                if f<4
                    bis=4-f; % RGB
                else
                    bis=f; % NIR, SWIR1, SWIR2
                end
                if f>4 % % Interpolate SWIR (20 m) to 10 m grid using NN.
                    temp=imresize(temp,[nr,nc],'nearest');
                end
                temp=temp(thesey,thesex); %  Crop
                temp(temp==sats2)=max(temp(temp~=sats2)); % Set saturated values to maximum reflectance.
                scaled=double(temp)./sfs2;
                 % Fill values will have ALL bands equal to 0.01 for
                 % surface reflectance.
                scaled(scaled<=0)=0.01;
                scaled(scaled>1)=1;
                ris(s).sf(bis)=max(scaled(:))./255; % Maximize use of bit range;
                ris(s).r(:,:,bis)=uint8(scaled./ris(s).sf(bis));  
                %ris(s).sf(bis)=1/255;
                %ris(s).r(:,:,bis)=uint8(scaled.*255);  
                
            end % <- End band (f) loop.
                ris(s).r(:,:,end)=scl(thesey,thesex); % Save the scene classification.
                
                
                ris(s).t=t;
                ris(s).order={'R'; 'G'; 'B'; 'NIR'; 'SWIR1'; 'SWIR2'};
                ris(s).SE=SE;
                ris(s).utmz=aoi.utmz;
                ris(s).delta=dx;
                ris(s).ulx=single(round(ulx)); % ULC identifiers for UNCROPPED scene.
                ris(s).uly=single(round(uly));    
                
        end % <- End scene (s) loop.  
     end % <- End read Sentinel-2 function.

    
    
    
    


end








