% This creates TS plots and maps of water masses with CHL above desired
% thresholds
% (Dataset ID: erdMH1chlamday_R2022SQ)
% Written by E. Montes 
% October 2023.

% create working data directory;
% close all; clear all;
clearvars -except lat lon chla time lat1 lon1 % % CLASS P time lat lon lat1 lon1 

%make download and save directory
path = ('I:\My Drive\GDrive\OCED_AOML\WS_cruises\sat_analysis'); 
cd(path)
addpath('I:\My Drive\GDrive\software\matlab\m_map_toolbox\m_map\m_map');
addpath('I:\My Drive\GDrive\software\matlab\error_envelope\kakearney-boundedline-pkg-5d00182\boundedline');
addpath('I:\My Drive\GDrive\software\matlab\satellite_weekly_monthly_means\Dan_programs');

% % USE WITH SFER_stations.csv
sfer_sta_list = readtable("SFER_stations.csv");
site_coord = [sfer_sta_list.decLat sfer_sta_list.decLon];
sta_ids = sfer_sta_list.Station;

% Florida
n_lat = '29'; s_lat = '24'; e_lon = '-79.5'; w_lon = '-86';

% % read Florida shp
addpath('G:\My Drive\MBON\FKNMS\GIS_data\fknms_py2')
S = shaperead('fknms_py.shp');
tmp_coord = [S.X; S.Y]; 
xx = tmp_coord(1,:)';
yy = tmp_coord(2,:)';

% % LOAD BATHYMETRY
% Load bathymetry data from the provided GeoTIFF file (download a sgeotiff
% at https://download.gebco.net/. Use the same LAT/LON limits as above!)
    bathymetry_data = geotiffread('gebco_bathy.tif');
    bathymetry_data = flipud(bathymetry_data); % Flip the data over the x-axis

% %%%%%%%%%%%%%% USE FOR erdMH1chla1day_R2022SQ FOR DAILY AND MONTHLY TIME SERIES FROM COASTWATCH - 2002-PRESENT ################################ 
% Florida (large region)
n_lat = '29'; s_lat = '24'; e_lon = '-79.5'; w_lon = '-86';

start_time = '2003-01';
end_time = '2023-08';

options = weboptions('Timeout',1200); % this command is necessary in how matlab looks for bigger data files
% % Monthly data
urlfile = strcat('https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chlamday_R2022SQ.nc?chlor_a%5B(',start_time,'-16T00:00:00Z):1:(',end_time,'-16T00:00:00Z)][(',n_lat,'):1:(',s_lat,')][(',w_lon,'):1:(',e_lon,')]');
% 8-day data (lastest layer: 2022-12)
% urlfile = strcat('https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chla8day_R202SQ.nc?chlor_a%5B(',start_time,'-31T00:00:00Z):1:(',end_time,'-31T00:00:00Z)%5D%5B(',n_lat,'):1:(',s_lat,')%5D%5B(',w_lon,'):1:(',e_lon,')%5D');
% % Daily data (don't use a very broad time range)
% urlfile = strcat('https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chla1day_R2022SQ.nc?chlor_a%5B(',start_time,'-31T12:00:00Z):1:(',end_time,'-31T12:00:00Z)][(',n_lat,'):1:(',s_lat,')][(',w_lon,'):1:(',e_lon,')]');

websave('dummy.nc',urlfile,options); % note in matlab this saves the url exactly as the file that it is, not a text or html
lat=double(ncread('dummy.nc','latitude'));
lon=double(ncread('dummy.nc','longitude'));
chla=double(ncread('dummy.nc','chlor_a')); % Dataset ID: erdMH1chlamday_R2022SQ
time=double(ncread('dummy.nc', 'time'));
[lat1,lon1]=meshgrid(lat,lon);

% % MAP SST
%Define lat/lon min and max
latmin = min(lat);
latmax = max(lat);
lonmin = min(lon);
lonmax = max(lon);
LATLIMS=[latmin latmax];
LONLIMS=[lonmin lonmax];

figure();
    m_proj('equidistant cylindrical','lon',LONLIMS,'lat',LATLIMS);
    m_pcolor(lon1,lat1,chla(:, :, 2));
    shading flat;
    colormap(jet); 
    m_gshhs_i('patch',[0.5 0.5 0.5],'edgecolor','none');
    m_grid_v2('box','fancy','fontsize',12,'fontweight','bold');
    h = colorbar;
    scale = [0.1 20];
    caxis(scale);
    set(gca, 'colorscale', 'log');
    title('CHL', 'FontSize',14);
hold on
% Plot bathymetry contour line
    % Define the bathymetry level you want to contour (e.g., -100 meters)
    bathy_levels = [-50, -100, -200];
    % Create X and Y grids
    [lonGrid, latGrid] = meshgrid(linspace(lonmin, lonmax, size(bathymetry_data, 2)), linspace(latmin, latmax, size(bathymetry_data, 1)));
    m_contour(lonGrid, latGrid, bathymetry_data, [bathy_levels, bathy_levels], 'linecolor', [0.5 0.5 0.5], 'linewidth', 1.5);
% To plot stations
for n = 1:length(site_coord)
    [X, Y] = m_ll2xy(site_coord(n, 2), site_coord(n, 1));
    line(X, Y, 'marker', 'o', 'markersize', 10, 'linewidth', 0.5, 'color', 'k', 'MarkerFaceColor', 'white');
end
hold off

%time is seconds since Jan 01, 1970 in Coastwatch world
dnum=datenum([1970,1,1])+time./(60*60*24) ; % need to convert seconds to days
TIMEVEC=datevec(dnum);
doy=nan(length(time),1);
for i=1:length(doy);
    if TIMEVEC(i,2)==1;
        doy(i,1)=TIMEVEC(i,3);
    else;
        doy(i,1)=sum(eomday(TIMEVEC(i,1),1:TIMEVEC(i,2)-1),2)+TIMEVEC(i,3); 
        %eomday provides the total number days in a month and accounts for
        %leap years. This is an easy way to calculate the day of the year
        %so that we can calculate a series date for easy plotting
    end;
end;

%%create the series date
SDATE=TIMEVEC(:,1)+doy./366;

%%%%%%%%%%%%% CREATE MONTHLY SST CLIMATOLOGIES AND ANOMALIES FOR FLORIDA %%%%%%%%%%%%%%%%%
% Select year span to calculate climatologies
yr_min = 2003;
yr_max = 2012;
min_idx = find(TIMEVEC(:,1) == yr_min);
max_idx = find(TIMEVEC(:,1) == yr_max);
first = min_idx(1); 
last = max_idx(end);
tmp_arr = chla(:,:, first:last);  % sst array used for climatology
tmp_timevec = TIMEVEC(first:last,:);

zeros_mo = 12; 
zeros_lat = length(lat);
zeros_lon = length(lon);
mo_clim = zeros(zeros_lon,zeros_lat,zeros_mo); 
for q=1:12
    mo_idx = find(tmp_timevec(:,2) == q);
    if mo_idx > 0
        mo_chla = tmp_arr(:,:,mo_idx);
        mo_chla_mean = squeeze(nanmean(mo_chla,3));
        mo_clim(:,:,q) = mo_chla_mean;
    else
    end
end

%%%%%%%%%%%%% THIS GENERATES MONTHLY ANOMALIES %%%%%%%%%%%%%%%%%
yr_val = 2005; % selected year for anomaly calculation
mo_val = 09; % selected month for anomaly calculation
time_idx = find(TIMEVEC(:,1) == yr_val & TIMEVEC(:,2) == mo_val);
chla_mo = chla(:,:, time_idx); % monthly mean
clim_mo = squeeze(mo_clim(:,:, mo_val)); % selected monthly climatology 
chla_anom = chla_mo - clim_mo;

% Turn values above the threshold into NaN
threshold = 3;
chla_anom_filt = chla_anom;
chla_anom_filt(chla_anom_filt < threshold) = NaN;

addpath('I:\My Drive\GDrive\OCED_AOML\Omics\sat_analysis');
load('cmap');
figure();
    m_proj('equidistant cylindrical','lon',LONLIMS,'lat',LATLIMS);
    m_pcolor(lon1,lat1,chla_anom_filt);
    shading flat;
    colormap(cm); %%% Use for SST anomalies
    m_gshhs_i('patch',[0.5 0.5 0.5],'edgecolor','none');
    m_grid_v2('box','fancy','fontsize',12,'fontweight','bold');
    h=colorbar;
    scale = [-3 3];
    caxis(scale);
    title('CHL anomaly Sep 2005 (vs 2003-2012 baseline)', 'FontSize',14);
hold on
    % Plot bathymetry contour line
    % Define the bathymetry level you want to contour (e.g., -100 meters)
    bathy_levels = [-50, -100, -200];
    % Create X and Y grids
    [lonGrid, latGrid] = meshgrid(linspace(lonmin, lonmax, size(bathymetry_data, 2)), linspace(latmin, latmax, size(bathymetry_data, 1)));
    m_contour(lonGrid, latGrid, bathymetry_data, [bathy_levels, bathy_levels], 'linecolor', [0.5 0.5 0.5], 'linewidth', 1.5);
% % To plot stations
% for n = 1:length(site_coord)
%     [X, Y] = m_ll2xy(site_coord(n, 2), site_coord(n, 1));
%     line(X, Y, 'marker', 'o', 'markersize', 10, 'linewidth', 0.5, 'color', 'k', 'MarkerFaceColor', 'white');
% end
hold off

% % Compute a area time series of water mass with SST anomalies above or below certain threshold
threshold_anom_low = -3;
threshold_anom_high = 3;
threshold_low = 0.5;
threshold_high = 40;
pstep=mode(diff(lat)); % just determining the step as I am not certain it is in the metadata
km2p=cosd(lat1).*((pstep*110).^2);

for k=1:length(SDATE)
    
    % generate anomaly arrays
    chla_array = chla(:,:, k); % monthly mean
    month_val = TIMEVEC(k, 2);
    clim_array = squeeze(mo_clim(:,:, month_val)); % selected monthly climatology 
    chla_anom_tmp = chla_array - clim_array;

    % Calculate areas of low and high SST waters
    AREA_low_anom(k) = sum(km2p(chla_anom_tmp <= threshold_anom_low));
    AREA_high_anom(k) = sum(km2p(chla_anom_tmp >= threshold_anom_high));
    AREA_low_chla(k) = sum(km2p(chla_array <= threshold_low));
    AREA_high_chla(k) = sum(km2p(chla_array >= threshold_high));
end

period = [2003 2024];
figure()
subplot(4,1,1)
plot(SDATE, AREA_low_anom, 'b')
xlim(period)
subplot(4,1,2)
plot(SDATE, AREA_high_anom, 'r')
xlim(period)
subplot(4,1,3)
xlim(period)
plot(SDATE, AREA_low_chla, 'b', 'LineStyle','--')
xlim(period)
subplot(4,1,4)
plot(SDATE, AREA_high_chla, 'r', 'LineStyle','--')
xlim(period)

