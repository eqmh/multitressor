% This creates TS plots and maps of SST and CHL at desired locations 
% Written by E. Montes 
% January 2020.

% create working data directory;
close all; clear all;
% clearvars -except lat lon SST chla time lat1 lon1 % % CLASS P time lat lon lat1 lon1 

%make download and save directory
path = ('I:\My Drive\GDrive\proposals\2022_02_MultiStressor_NOAA\trends_analysis'); 
cd(path)
addpath('I:\My Drive\GDrive\software\matlab\m_map_toolbox\m_map\m_map');
addpath('I:\My Drive\GDrive\software\matlab\error_envelope\kakearney-boundedline-pkg-5d00182\boundedline');
addpath('I:\My Drive\GDrive\software\matlab\satellite_weekly_monthly_means\Dan_programs');
addpath('I:\My Drive\GDrive\software\matlab\anomaly');

% % % Set up the Import Options and import the data
% opts = delimitedTextImportOptions("NumVariables", 9);
% 
% % Specify range and delimiter
% opts.DataLines = [1, Inf];
% opts.Delimiter = "\t";
% 
% % Specify column names and types
% opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "LatitudeDecimal", "Var6", "Var7", "LongitudeDecimal", "Station"];
% opts.SelectedVariableNames = ["LatitudeDecimal", "LongitudeDecimal", "Station"];
% opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string"];
% 
% % Specify file level properties
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% 
% % Specify variable properties
% opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "LatitudeDecimal", "Var6", "Var7", "LongitudeDecimal", "Station"], "WhitespaceRule", "preserve");
% opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "LatitudeDecimal", "Var6", "Var7", "LongitudeDecimal", "Station"], "EmptyFieldRule", "auto");
% 
% % Import the data
% SFERsites = readmatrix("I:\My Drive\GDrive\OCED_AOML\WS_cruises\sat_analysis\SFER_sites.txt", opts);
% 
% %  % Clear temporary variables
% clear opts

% % % Add BG line extension sites
% bg_ext = readtable('bg_extension_line.csv');
% bg_ext_lat = bg_ext.latitude;
% bg_ext_lon = bg_ext.longitude;

% % % SFER sites (row 4673 is the first site sampled in 2022)
% sta_ids = unique(SFERsites(4673:end, 3));
% % % transform lat/lon from string to numeric
% SFERsites_sel = SFERsites(4673:end, :);
% for f=1:size(SFERsites_sel, 1)
%     dec_lat(f,1) = str2num(SFERsites_sel(f, 1));
%     dec_lon(f,1) = str2num(SFERsites_sel(f, 2));
% end
% 
% % % calculate mean/sd of station lat/lon
% for n=1:length(sta_ids)
%     logic_idx = strcmp(SFERsites_sel(:, 3), sta_ids(n));
%     row_idx = find(logic_idx);
%     lat_avg(n,1) = mean(dec_lat(row_idx));
%     lon_avg(n,1) = mean(dec_lon(row_idx));
%     lat_sd(n,1) = std(dec_lat(row_idx));
%     lon_std(n,1) = std(dec_lon(row_idx));
% end
% 
% lat_array = [lat_avg; bg_ext_lat];
% lon_array = [lon_avg; bg_ext_lon];
% site_coord = [lat_array lon_array];

% % USE WITH SFER_stations.csv
sfer_sta_list = readtable("SFER_stations.csv");
site_coord = [sfer_sta_list.decLat sfer_sta_list.decLon];
sta_ids = sfer_sta_list.Station;

% Florida
n_lat = '29'; s_lat = '24'; e_lon = '-79.5'; w_lon = '-86';

% % read Florida shp
addpath('I:\My Drive\GDrive\MBON\FKNMS\GIS_data\fknms_py2')
S = shaperead('fknms_py.shp');
tmp_coord = [S.X; S.Y]; 
xx = tmp_coord(1,:)';
yy = tmp_coord(2,:)';

% % LOAD BATHYMETRY
% Load bathymetry data from the provided GeoTIFF file (download a sgeotiff
% at https://download.gebco.net/. Use the same LAT/LON limits as above!)
    bathymetry_data = geotiffread('gebco_bathy.tif');
    bathymetry_data = flipud(bathymetry_data); % Flip the data over the x-axis
    
%%
% %%%%%%%%%%%%%% USE FOR CHL FROM ERDDAP - 2003-PRESENT ################################ 
%%% DEFINE PERIOD AND REGION%%%%
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

chla_avg = nanmean(chla, 3);

% % Index pixel values within complex polygons
% [IN]= inpolygon(lon1, lat1,xx, yy);
% chla_avg(~IN) = NaN;
    
%Define lat/lon min and max
latmin = min(lat);
latmax = max(lat);
lonmin = min(lon);
lonmax = max(lon);
LATLIMS=[latmin latmax];
LONLIMS=[lonmin lonmax];

figure;
    m_proj('equidistant cylindrical', 'lon', LONLIMS, 'lat', LATLIMS);
    m_pcolor(lon1, lat1, chla(:, :, 189));
    shading flat;
    colormap(jet);
    m_gshhs_i('patch', [0.5 0.5 0.5], 'edgecolor', 'none');
    m_grid_v2('box', 'fancy', 'fontsize', 12, 'fontweight', 'bold');
    h = colorbar;
    scale = [0.1 20];
    caxis(scale);
    set(gca, 'colorscale', 'log');
    title('CHL September 2018', 'FontSize', 16);
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

%%%%%%
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

%%% CALCULATES MEAN SST WITH SD WITHIN BOX FOR SELECTED STATIONS %%%%
% % selects sites to be analyzed
sta_ids_sel = {'MR'; 'LK'; 'WS'; 'CAL6'; '30'; 'TB5'; 'BG20'};
conc_chl_mean = [];
conc_chl_sd = [];
for f=1:length(sta_ids_sel)
    logic_idx = strcmp(sta_ids, sta_ids_sel(f));
    sta_idx = find(logic_idx);
    lat_pt = site_coord(sta_idx, 1);
    lon_pt = site_coord(sta_idx, 2);
    [loc_pix,loc_line]=latlon2pixline(lat_pt,lon_pt,lat,lon); % loc_pix=lon and loc_line=lat
    
    box=2;%must be an even number
    x1=loc_pix-box/2;
    x2=loc_pix+box/2;
    y1=loc_line-box/2;
    y2=loc_line+box/2;

    conc_chl = [];
        for j=1:size(chla,3)
            tmp_chl = squeeze(chla(:,:,j));
            tmp2_chl = tmp_chl(x1:x2, y1:y2);
            mean_chl = mean(tmp2_chl(~isnan(tmp2_chl)));
            std_chl = std(tmp2_chl(~isnan(tmp2_chl)));
            
            conc_chl = [conc_chl; mean_chl, std_chl];
        end

    conc_chl_mean = [conc_chl_mean, conc_chl(:, 1)]; 
    conc_chl_sd = [conc_chl_sd, conc_chl(:, 2)];
end

num_plots = length(sta_ids_sel);
period = [2003 2024];
baseline_top = zeros(length(conc_chl_mean),1)+2;
baseline_bottom = zeros(length(conc_chl_mean),1)+0;
% Make plots
figure()
for w=1:num_plots
    subplot(num_plots,1,w)
        ax(2) = subplot(num_plots,1,w);
        boundedline(SDATE(:), conc_chl_mean(:,w), conc_chl_sd(:,w), 'b', 'alpha', 'transparency', 0.10);
        hold on
        plot(SDATE, baseline_top, 'r')
        plot(SDATE, baseline_bottom, 'r')
        xlim(period)
        ylim([0 4])
        % set(gca,'fontsize',14);
        ylabel('Chl-a (mg m-3)')
%         title(site_names(w)) % % for MIR analysis
        title(sta_ids_sel(w), 'FontSize', 16)  % % for Outplant analysis
        hold off
end
 
% chl_outplant = [SDATE conc_chl_mean conc_chl_sd];
% writematrix(chl_outplant, 'chl_outplant_8day.csv');
% str_date = datestr(dnum, 'yyyy-mm-dd');
% dlmwrite('str_date2019_chl.txt', str_date, 'delimiter', '');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% %%%%%%%%%%%%%% USE FOR jplMURSST41 FOR DAILY AND MONTHLY TIME SERIES FROM COASTWATCH - 2002-PRESENT ################################ 
% Florida (large region)
n_lat = '29'; s_lat = '24'; e_lon = '-79.5'; w_lon = '-86';

start_time = '2003-01';
end_time = '2023-09';

% % DAILY
% urlfile = strcat('https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41.nc?analysed_sst[(',start_time,'T09:00:00Z):1:(',end_time,'T09:00:00Z)][(',s_lat,'):1:(',n_lat,')][(',w_lon,'):1:(',e_lon,')]');

% % MONTHLY
urlfile = strcat('https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41mday.nc?sst%5B(',start_time,'-16T00:00:00Z):1:(',end_time,'-16T00:00:00Z)%5D%5B(',s_lat,'):1:(',n_lat,')%5D%5B(',e_lon,'):1:(',w_lon,')%5D,nobs%5B(',start_time,'-16T00:00:00Z):1:(',end_time,'-16T00:00:00Z)%5D%5B(',s_lat,'):1:(',n_lat,')%5D%5B(',e_lon,'):1:(',w_lon,')%5D,mask%5B(',start_time,'-16T00:00:00Z):1:(',end_time,'-16T00:00:00Z)%5D%5B(',s_lat,'):1:(',n_lat,')%5D%5B(',e_lon,'):1:(',w_lon,')%5D');

options = weboptions('Timeout',1200); % this command is necessary in how matlab looks for bigger data files
websave('dummy.nc',urlfile,options); % note in matlab this saves the url exactly as the file that it is, not a text or html
lat=double(ncread('dummy.nc','latitude'));
lon=double(ncread('dummy.nc','longitude'));
% SST=double(ncread('dummy.nc','analysed_sst'));
SST=double(ncread('dummy.nc','sst'));
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
    m_pcolor(lon1,lat1,SST(:, :, end));
    shading flat;
    colormap(jet); 
    m_gshhs_i('patch',[0.5 0.5 0.5],'edgecolor','none');
    m_grid_v2('box','fancy','fontsize',12,'fontweight','bold');
    h = colorbar;
    scale = [25 32];
    caxis(scale);
    title('SST September 2023', 'FontSize',14);
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

%%% CALCULATES MEAN SST WITH SD WITHIN BOX FOR SELECTED STATIONS %%%%
% % selects sites to be analyzed
sta_ids_sel = {'MR'; 'LK'; 'WS'; 'CAL6'; '30'; 'TB5'; 'BG20'};
conc_sst_mean = [];
conc_sst_sd = [];
for f=1:length(sta_ids_sel)
    logic_idx = strcmp(sta_ids, sta_ids_sel(f));
    sta_idx = find(logic_idx);
    lat_pt = site_coord(sta_idx, 1);
    lon_pt = site_coord(sta_idx, 2);
    [loc_pix,loc_line]=latlon2pixline(lat_pt,lon_pt,lat,lon); % loc_pix=lon and loc_line=lat
    
    box=2;%must be an even number
    x1=loc_pix-box/2;
    x2=loc_pix+box/2;
    y1=loc_line-box/2;
    y2=loc_line+box/2;

    conc_sst = [];
        for j=1:size(SST,3)
            tmp_sst = squeeze(SST(:,:,j));
            tmp2_sst = tmp_sst(x1:x2, y1:y2);
            mean_sst = mean(tmp2_sst(~isnan(tmp2_sst)));
            std_sst = std(tmp2_sst(~isnan(tmp2_sst)));
            
            conc_sst = [conc_sst; mean_sst, std_sst];
        end

    conc_sst_mean = [conc_sst_mean, conc_sst(:, 1)]; 
    conc_sst_sd = [conc_sst_sd, conc_sst(:, 2)];
end

% % CREATE PLOTS
num_plots = length(sta_ids_sel);
period = [2003 2024];
baseline_top = zeros(length(conc_sst_mean),1)+2;
baseline_bottom = zeros(length(conc_sst_mean),1)+0;
% Make plots
figure()
for w=1:num_plots
    subplot(num_plots,1,w)
        ax(2) = subplot(num_plots,1,w);
        boundedline(SDATE(:), conc_sst_mean(:,w), conc_sst_sd(:,w), 'b', 'alpha', 'transparency', 0.10);
        hold on
        plot(SDATE, baseline_top, 'r')
        plot(SDATE, baseline_bottom, 'r')
        xlim(period)
        ylim([17 33])
        % set(gca,'fontsize',14);
        ylabel('SST degC')
%         title(site_names(w)) % % for MIR analysis
        title(sta_ids_sel(w), 'FontSize', 16)  % % for Outplant analysis
        hold off
end
 
% sst_outplant = mir_sst_all;
% writematrix(sst_outplant, 'sst_outplant_daily.csv');
% str_date = datestr(dnum, 'yyyy-mm-dd');
% dlmwrite('str_date2019_sst.txt', str_date, 'delimiter', '');

%%
%%%%%%%%%%%%% CREATE MONTHLY SST CLIMATOLOGIES AND ANOMALIES FOR FLORIDA %%%%%%%%%%%%%%%%%
% Select year span to calculate climatologies
yr_min = 2003;
yr_max = 2012;
min_idx = find(TIMEVEC(:,1) == yr_min);
max_idx = find(TIMEVEC(:,1) == yr_max);
first = min_idx(1); 
last = max_idx(end);
tmp_arr = SST(:,:, first:last);  % sst array used for climatology
tmp_timevec = TIMEVEC(first:last,:);

zeros_mo = 12; 
zeros_lat = length(lat);
zeros_lon = length(lon);
mo_clim = zeros(zeros_lon,zeros_lat,zeros_mo); 
for q=1:12
    mo_idx = find(tmp_timevec(:,2) == q);
    if mo_idx > 0
        mo_sst = tmp_arr(:,:,mo_idx);
        mo_sst_mean = squeeze(nanmean(mo_sst,3));
        mo_clim(:,:,q) = mo_sst_mean;
    else
    end
end

%%%%%%%%%%%%% THIS GENERATES MONTHLY ANOMALIES %%%%%%%%%%%%%%%%%
yr_val = 2010; % selected year for anomaly calculation
mo_val = 02; % selected month for anomaly calculation
time_idx = find(TIMEVEC(:,1) == yr_val & TIMEVEC(:,2) == mo_val);
sst_mo = SST(:,:, time_idx); % monthly mean
clim_mo = squeeze(mo_clim(:,:, mo_val)); % selected monthly climatology 
sst_anom = sst_mo - clim_mo;

addpath('I:\My Drive\GDrive\OCED_AOML\Omics\sat_analysis');
load('cmap');
figure();
    m_proj('equidistant cylindrical','lon',LONLIMS,'lat',LATLIMS);
    m_pcolor(lon1,lat1,sst_anom);
    shading flat;
    colormap(cm); %%% Use for SST anomalies
    m_gshhs_i('patch',[0.5 0.5 0.5],'edgecolor','none');
    m_grid_v2('box','fancy','fontsize',12,'fontweight','bold');
    h=colorbar;
    scale = [-2 2];
    caxis(scale);
    title('SST anomaly February 2010 (vs 2003-2012 baseline)', 'FontSize',14);
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

%%% CALCULATES MONTLY SST ANOM VALUES WITHIN BOX %%%%
% % Compute SST time series
sta_ids_sel = {'MR'; 'LK'; 'WS'; 'CAL6'; '30'; 'TB5'; 'BG20'};
conc_sst_mean = [];
conc_sst_sd = [];
for f=1:length(sta_ids_sel)
    logic_idx = strcmp(sta_ids, sta_ids_sel(f));
    sta_idx = find(logic_idx);
    lat_pt = site_coord(sta_idx, 1);
    lon_pt = site_coord(sta_idx, 2);
    [loc_pix,loc_line]=latlon2pixline(lat_pt,lon_pt,lat,lon); % loc_pix=lon and loc_line=lat
    
    box=10;%must be an even number
    x1=loc_pix-box/2;
    x2=loc_pix+box/2;
    y1=loc_line-box/2;
    y2=loc_line+box/2;

    conc_sst = [];
        for j=1:size(SST,3)
            tmp_sst = squeeze(SST(:,:,j));
            tmp2_sst = tmp_sst(x1:x2, y1:y2);
            mean_sst = mean(tmp2_sst(~isnan(tmp2_sst)));
            std_sst = std(tmp2_sst(~isnan(tmp2_sst)));
            
            conc_sst = [conc_sst; mean_sst, std_sst];
        end

    conc_sst_mean = [conc_sst_mean, conc_sst(:, 1)]; 
    conc_sst_sd = [conc_sst_sd, conc_sst(:, 2)];
end

% % Compute SST climatological value within each box
conc_clim_sst_val = [];
for f=1:length(sta_ids_sel)
    logic_idx = strcmp(sta_ids, sta_ids_sel(f));
    sta_idx = find(logic_idx);
    lat_pt = site_coord(sta_idx, 1);
    lon_pt = site_coord(sta_idx, 2);
    [loc_pix,loc_line]=latlon2pixline(lat_pt,lon_pt,lat,lon); % loc_pix=lon and loc_line=lat
    
    box=10;%must be an even number
    x1=loc_pix-box/2;
    x2=loc_pix+box/2;
    y1=loc_line-box/2;
    y2=loc_line+box/2;

    clim_sst_val = [];
    for k=1:12 
        clim_val = squeeze(mo_clim(:,:, k));
        tmp_clim_arr = clim_val(x1:x2, y1:y2);
        tmp_clim = mean(tmp_clim_arr(~isnan(tmp_clim_arr)));
        clim_sst_val = [clim_sst_val; tmp_clim];
    end
    conc_clim_sst_val = [conc_clim_sst_val, clim_sst_val];
end

%%% Compute SST anomaly within each box
conc_anom_sst_ts = [];
for f=1:length(sta_ids_sel)
    anom_sst_ts = [];
    for n=1:size(conc_sst_mean, 1)      
        mm_idx = TIMEVEC(n,2); 
        clim_val = conc_clim_sst_val(mm_idx,f);
        delta_sst = conc_sst_mean(n,f) - clim_val;
        anom_sst_ts = [anom_sst_ts; delta_sst];        
    end

    conc_anom_sst_ts = [conc_anom_sst_ts, anom_sst_ts];
end

%%% PLOT SST MEAN AND ANOMALY WITH SD WITHIN BOX %%%%
figure()
plot(SDATE, conc_sst_mean(:,:))
set(gca,'fontsize',14)
ylim([19 33])
xlim([2002 2024])
ylabel('Mean monthly SST (C)')
legend(sta_ids_sel)

%%%% Plot anomalies  %%%%%%
num_plots = length(sta_ids_sel);
period = [2003 2024];
baseline_zero = zeros(length(conc_anom_sst_ts),1);
figure()
for w=1:num_plots
    subplot(num_plots,1,w)
        anomaly(SDATE, conc_anom_sst_ts(:, w));
        xlim(period)
        ylim([-3 3])
        set(gca,'fontsize',10);
        ylabel('SST anom degC')
        title(sta_ids_sel(w))
        hold off
end

% sst_mir_anom = [SDATE conc_anom_sst_ts];
% writematrix(sst_mir_anom, 'sst_mir_anom2.csv');
% str_date = datestr(dnum, 'yyyy-mm-dd');
% dlmwrite('str_date2.txt', str_date, 'delimiter', '');

% % Plot SST anomalies for specific months
sel_mo_idx = find(TIMEVEC(:, 2) == 6);
sel_dates = SDATE(sel_mo_idx);
sel_sst_anom = conc_anom_sst_ts(sel_mo_idx, :);
figure()
for w=1:num_plots
    subplot(num_plots,1,w)
        bar(sel_dates, sel_sst_anom(:, w));
        xlim(period)
        ylim([-1.5 1.5])
        set(gca,'fontsize',12);
        ylabel('SST anomaly (C)')
        title(sta_ids_sel(w))
        hold off
end

% % COMPARE ANOMALIES WITH GHRSST ANOMALY PRODUCT
% Florida (large region)
n_lat = '29'; s_lat = '24'; e_lon = '-79.5'; w_lon = '-86';

start_time = '2010-02';
end_time = '2010-02';

% % MONTHLY ANOMALY (https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41anommday.graph)
urlfile = strcat('https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41anommday.nc?sstAnom%5B(',start_time,'-16T00:00:00Z):1:(',end_time,'-16T00:00:00Z)%5D%5B(',s_lat,'):1:(',n_lat,')%5D%5B(',e_lon,'):1:(',w_lon,')%5D');
options = weboptions('Timeout',1200); % this command is necessary in how matlab looks for bigger data files
websave('dummy.nc',urlfile,options); % note in matlab this saves the url exactly as the file that it is, not a text or html
lat=double(ncread('dummy.nc','latitude'));
lon=double(ncread('dummy.nc','longitude'));
SST_anom=double(ncread('dummy.nc','sstAnom'));
time=double(ncread('dummy.nc', 'time'));
[lat1,lon1]=meshgrid(lat,lon);

figure();
    m_proj('equidistant cylindrical','lon',LONLIMS,'lat',LATLIMS);
    m_pcolor(lon1,lat1,SST_anom);
    shading flat;
    colormap(cm); %%% Use for SST anomalies
    m_gshhs_i('patch',[0.5 0.5 0.5],'edgecolor','none');
    m_grid_v2('box','fancy','fontsize',12,'fontweight','bold');
    h=colorbar;
    scale = [-2 2];
    caxis(scale);

%%
%%%%%%%%%%%%% CREATE MONTHLY CHL CLIMATOLOGIES AND ANOMALIES FOR FLORIDA %%%%%%%%%%%%%%%%%
% % NOTE: DON'T RUN SST CODE ABOVE, JUST CHL-A SECTION!
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
mo_clim_chla = zeros(zeros_lon,zeros_lat,zeros_mo); 
for q=1:12
    mo_idx = find(tmp_timevec(:,2) == q);
    if mo_idx > 0
        mo_chla = tmp_arr(:,:,mo_idx);
        mo_chla_mean = squeeze(nanmean(mo_chla,3));
        mo_clim_chla(:,:,q) = mo_chla_mean;
    else
    end
end

%%%%%%%%%%%%% THIS GENERATES MONTHLY CHLA ANOMALIES %%%%%%%%%%%%%%%%%
yr_val = 2023; % selected year for anomaly calculation
mo_val = 08; % selected month for anomaly calculation
time_idx = find(TIMEVEC(:,1) == yr_val & TIMEVEC(:,2) == mo_val);
chla_mo = chla(:,:, time_idx); % monthly value
clim_mo_chla = squeeze(mo_clim_chla(:,:, mo_val)); % selected monthly climatology 
chla_anom = chla_mo - clim_mo_chla;

addpath('I:\My Drive\GDrive\OCED_AOML\Omics\sat_analysis');
load('cmap');
figure();
    m_proj('equidistant cylindrical','lon',LONLIMS,'lat',LATLIMS);
    m_pcolor(lon1,lat1,chla_anom);
    shading flat;
    colormap(cm); %%% Use cm for CHLA anomalies
    m_gshhs_i('patch',[0.5 0.5 0.5],'edgecolor','none');
    m_grid_v2('box','fancy','fontsize',12,'fontweight','bold');
    h=colorbar;
    scale = [-1 1];
    caxis(scale);
%     set(gca, 'colorscale', 'log')
    title('CHL anomaly September 2018 (vs 2003-2012 baseline)', 'FontSize',16);
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

%%% CALCULATES MONTLY CHL ANOM VALUES WITHIN BOX %%%%
% % Compute SST time series
sta_ids_sel = {'MR'; 'LK'; 'WS'; 'CAL6'; '30'; 'TB5'; 'BG20'};
conc_chla_mean = [];
conc_chla_sd = [];
for f=1:length(sta_ids_sel)
    logic_idx = strcmp(sta_ids, sta_ids_sel(f));
    sta_idx = find(logic_idx);
    lat_pt = site_coord(sta_idx, 1);
    lon_pt = site_coord(sta_idx, 2);
    [loc_pix,loc_line]=latlon2pixline(lat_pt,lon_pt,lat,lon); % loc_pix=lon and loc_line=lat
    
    box=2;%must be an even number
    x1=loc_pix-box/2;
    x2=loc_pix+box/2;
    y1=loc_line-box/2;
    y2=loc_line+box/2;

    conc_chla = [];
        for j=1:size(chla,3)
            tmp_chla = squeeze(chla(:,:,j));
            tmp2_chla = tmp_chla(x1:x2, y1:y2);
            mean_chla = mean(tmp2_chla(~isnan(tmp2_chla)));
            std_chla = std(tmp2_chla(~isnan(tmp2_chla)));
            
            conc_chla = [conc_chla; mean_chla, std_chla];
        end

    conc_chla_mean = [conc_chla_mean, conc_chla(:, 1)]; 
    conc_chla_sd = [conc_chla_sd, conc_chla(:, 2)];
end

% % Compute CHL climatological value within each box
conc_clim_chla_val = [];
for f=1:length(sta_ids_sel)
        logic_idx = strcmp(sta_ids, sta_ids_sel(f));
    sta_idx = find(logic_idx);
    lat_pt = site_coord(sta_idx, 1);
    lon_pt = site_coord(sta_idx, 2);
    [loc_pix,loc_line]=latlon2pixline(lat_pt,lon_pt,lat,lon); % loc_pix=lon and loc_line=lat
    
    box=2;%must be an even number
    x1=loc_pix-box/2;
    x2=loc_pix+box/2;
    y1=loc_line-box/2;
    y2=loc_line+box/2;

    clim_chla_val = [];
    for k=1:12 
        clim_val_chla = squeeze(mo_clim_chla(:,:, k));
        tmp_clim_arr = clim_val_chla(x1:x2, y1:y2);
        tmp_clim = mean(tmp_clim_arr(~isnan(tmp_clim_arr)));
        clim_chla_val = [clim_chla_val; tmp_clim];
    end
    conc_clim_chla_val = [conc_clim_chla_val, clim_chla_val];
end

%%% Compute CHL anomaly within each box
conc_anom_chla_ts = [];
for f=1:length(sta_ids_sel)
    anom_chla_ts = [];
    for n=1:size(conc_chla_mean, 1)      
        mm_idx = TIMEVEC(n,2); 
        clim_val_chla = conc_clim_chla_val(mm_idx,f);
        delta_chla = conc_chla_mean(n,f) - clim_val_chla;
        anom_chla_ts = [anom_chla_ts; delta_chla];        
    end

    conc_anom_chla_ts = [conc_anom_chla_ts, anom_chla_ts];
end

%%% PLOT CHL MEAN AND ANOMALY WITH SD WITHIN BOX %%%%
figure()
plot(SDATE, conc_chla_mean(:,:))
set(gca,'fontsize',14)
ylim([0 4])
xlim([2002 2024])
ylabel('Mean monthly Chl-a (mg m-3)')
legend(sta_ids_sel)

%%%% Plot anomalies  %%%%%%
num_plots = length(sta_ids_sel);
period = [2003 2024];
baseline_zero = zeros(length(conc_anom_chla_ts),1);
figure()
for w=1:num_plots
    subplot(num_plots,1,w)
        anomaly(SDATE, conc_anom_chla_ts(:, w));
        xlim(period)
        ylim([-1 1])
        set(gca,'fontsize',10);
        ylabel('Chla anom mg m-3')
        title(sta_ids_sel(w))
        hold off
end

% sst_mir_anom = [SDATE conc_anom_sst_ts];
% writematrix(sst_mir_anom, 'sst_mir_anom2.csv');
% str_date = datestr(dnum, 'yyyy-mm-dd');
% dlmwrite('str_date2.txt', str_date, 'delimiter', '');

% % Plot CHL anomalies for specific months
sel_mo_idx = find(TIMEVEC(:, 2) == 6);
sel_dates = SDATE(sel_mo_idx);
sel_chla_anom = conc_anom_chla_ts(sel_mo_idx, :);
figure()
for w=1:num_plots
    subplot(num_plots,1,w)
        bar(sel_dates, sel_chla_anom(:, w));
        xlim(period)
        ylim([-1.5 1.5])
        set(gca,'fontsize',12);
        ylabel('Chla anomaly (mg m-3)')
        title(sta_ids_sel(w))
        hold off
end




