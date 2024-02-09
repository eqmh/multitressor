% This script creates Mann-Kendall trend analysis using SST data and maps
% showing SST trends in selected regions of interest
% Written by E. Montes 
% October 2023.

% create working data directory;
close all; clear all;
% clearvars -except lat lon SST time lat1 lon1 % % CLASS P time lat lon lat1 lon1 

%make download and save directory
path = ('I:\My Drive\GDrive\proposals\2022_02_MultiStressor_NOAA\trends_analysis'); 
cd(path)
addpath('I:\My Drive\GDrive\software\matlab\m_map_toolbox\m_map\m_map');
addpath('I:\My Drive\GDrive\software\matlab\satellite_weekly_monthly_means\Dan_programs');
addpath('I:\My Drive\GDrive\proposals\2022_02_MultiStressor_NOAA\trends_analysis\cmocean'); 

%%
% % %%%%%%%%%%%%%%%%%%%%%%% FOR SST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Florida
n_lat = '29'; s_lat = '24'; e_lon = '-79.5'; w_lon = '-86';

start_time = '2003-01';
end_time = '2023-09';

% % Extract monthly SST 
% urlfile = strcat('https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41mday.nc?sst%5B(',start_time,'-16T00:00:00Z):1:(',end_time,'-16T00:00:00Z)%5D%5B(',s_lat,'):1:(',n_lat,')%5D%5B(',w_lon,'):1:(',e_lon,')%5D');

% Extract monthly SST anomaly
urlfile = strcat('https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41anommday.nc?sstAnom%5B(',start_time,'-16T00:00:00Z):1:(',end_time,'-16T00:00:00Z)%5D%5B(',s_lat,'):1:(',n_lat,')%5D%5B(',w_lon,'):1:(',e_lon,')%5D')

options = weboptions('Timeout',1200); % this command is necessary in how matlab looks for bigger data files
websave('dummy.nc',urlfile,options); % note in matlab this saves the url exactly as the file that it is, not a text or html
lat=double(ncread('dummy.nc','latitude'));
lon=double(ncread('dummy.nc','longitude'));
% SST=double(ncread('dummy.nc','sst'));
SST=double(ncread('dummy.nc','sstAnom'));
time=double(ncread('dummy.nc', 'time'));
[lat1,lon1]=meshgrid(lat,lon);

% % MAP 
%Define lat/lon min and max
latmin = min(lat);
latmax = max(lat);
lonmin = min(lon);
lonmax = max(lon);
LATLIMS=[latmin latmax];
LONLIMS=[lonmin lonmax];

selectet = SST(:,:,189);

figure();
    m_proj('equidistant cylindrical','lon',LONLIMS,'lat',LATLIMS);
    m_pcolor(lon1,lat1,selectet);
    shading flat;
    colormap(jet); 
    m_gshhs_i('patch',[0.5 0.5 0.5],'edgecolor','none');
    m_grid_v2('box','fancy','fontsize',12,'fontweight','bold');
    h = colorbar;
%     scale = [26 30];
%     caxis(scale);
    title('test', 'FontSize',16);

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

% % Generate per-pixel Mann-Kendall and linear regression statistics 
variable = SST;
for j=1:size(variable,1)
    for k=1:size(variable,2)
        pix_data = squeeze(variable(j,k,:));
            if all(isnan(mean(pix_data))) == 1
                mk = nan;
            elseif all(isnan(mean(pix_data))) == 0
        % the variable z is the Mann-Kendall trend statistic (MK).  
        % Larger absolute values of MK suggest stronger trends.
        [trend, z, p_value] = mann_kendall(pix_data);
                if p_value > 0.0001
                    mk = nan;
                elseif p_value < 0.0001
                    mk = z;
                end
         % Fit a linear model (polynomial of degree 1)
         p = polyfit(SDATE, pix_data, 1);
         slope = p(1); % Extract the slope and intercept from the polynomial coefficients
         residuals = pix_data - polyval(p, SDATE); 
         residualSE = sqrt(sum(residuals.^2) / (length(pix_data) - 2));
         t_statistic = slope / (residualSE / sqrt(sum((SDATE - mean(SDATE)).^2)));
         p_value_lst_sqr = 2 * (1 - tcdf(abs(t_statistic), length(pix_data) - 2));
                if p_value_lst_sqr > 0.001
                    lst_sqr = nan;
                elseif p_value_lst_sqr < 0.001
                    lst_sqr = slope;
                end
            end
        mk_map(j,k) = mk;
        lst_sqr_map(j,k) = lst_sqr;
    end
end

selected_map = lst_sqr_map;
cm3 = cmocean('thermal');
figure();
    m_proj('equidistant cylindrical','lon',LONLIMS,'lat',LATLIMS);
    m_pcolor(lon1,lat1,selected_map);
    shading flat;
    colormap(cm3);
    hold on
    [C, hContour] = m_contour(lon1, lat1, selected_map, [0.05,0.06], 'LineColor', 'r', 'LineWidth', 2);
    hold off
    m_gshhs_i('patch',[0.5 0.5 0.5],'edgecolor','none');
    m_grid_v2('box','fancy','fontsize',12,'fontweight','bold');
    h = colorbar;
    h.FontSize = 12; % Adjust the size as needed
    set(h, 'Ticks', [0.02:0.01:0.06], 'TickLabels', cellstr(num2str([0.02:0.01:0.06].')));
%     scale = [26 28];
%     caxis(scale);
    title('Monthly SST linear trend (p < 0.001)', 'FontSize',16);
    clabel(C, hContour, 'LabelSpacing', 100); % Adjust LabelSpacing as needed
    ylabel(h, 'Slope', 'FontSize', 14, 'Rotation', 90); % Change 'Your Colorbar Label' to your actual label

%%
% % %%%%%%%%%%%%%%%%%% FOR CHL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;

% Florida
n_lat = '29'; s_lat = '24'; e_lon = '-79.5'; w_lon = '-86';

start_time = '2003-01';
end_time = '2023-09';

% Extract monthly CHL
urlfile = strcat('https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1chlamday_R2022SQ.nc?chlor_a%5B(',start_time,'-16T00:00:00Z):1:(',end_time,'-16T00:00:00Z)][(',n_lat,'):1:(',s_lat,')][(',w_lon,'):1:(',e_lon,')]');

options = weboptions('Timeout',1200); % this command is necessary in how matlab looks for bigger data files
websave('dummy.nc',urlfile,options); % note in matlab this saves the url exactly as the file that it is, not a text or html
lat=double(ncread('dummy.nc','latitude'));
lon=double(ncread('dummy.nc','longitude'));
chla=double(ncread('dummy.nc','chlor_a')); % Dataset ID: erdMH1chlamday_R2022SQ
time=double(ncread('dummy.nc', 'time'));
[lat1,lon1]=meshgrid(lat,lon);

% % MAP 
%Define lat/lon min and max
latmin = min(lat);
latmax = max(lat);
lonmin = min(lon);
lonmax = max(lon);
LATLIMS=[latmin latmax];
LONLIMS=[lonmin lonmax];

selected = chla(:,:,189);

% figure;
%     m_proj('equidistant cylindrical', 'lon', LONLIMS, 'lat', LATLIMS);
%     m_pcolor(lon1, lat1, selected);
%     shading flat;
%     colormap(jet);
%     m_gshhs_i('patch', [0.5 0.5 0.5], 'edgecolor', 'none');
%     m_grid_v2('box', 'fancy', 'fontsize', 12, 'fontweight', 'bold');
%     h = colorbar;
%     scale = [0.1 20];
%     caxis(scale);
%     set(gca, 'colorscale', 'log');
%     title('CHL September 2018', 'FontSize', 16);

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

% Select year span to calculate climatologies
yr_min = 2003;
yr_max = 2012;
min_idx = find(TIMEVEC(:,1) == yr_min);
max_idx = find(TIMEVEC(:,1) == yr_max);
first = min_idx(1); 
last = max_idx(end);
tmp_arr = chla(:,:, first:last);  % sst array used for climatology
tmp_timevec = TIMEVEC(first:last,:);

% generate chl climatologies for the entire domain
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

% % Creates monthly CHL anomaly time series array for the entire domain
for t=1:size(chla,3)
    month_val = TIMEVEC(t,2);
    chla_mo = chla(:,:, t); % monthly array
    clim_mo_chla = squeeze(mo_clim_chla(:,:, month_val)); % selected monthly climatology 
    chla_anom(:,:,t) = chla_mo - clim_mo_chla;
end

% % Generate per-pixel Mann-Kendall statistics 
variable = chla_anom;
for j=1:size(variable,1)
    for k=1:size(variable,2)
        pix_data = squeeze(variable(j,k,:));
            if all(isnan(mean(pix_data))) == 1
                mk = nan;
            elseif all(isnan(mean(pix_data))) == 0
        % the variable z is indeed the Mann-Kendall trend statistic (MK).  
        % Larger absolute values of MK suggest stronger trends.
        [trend, z, p_value] = mann_kendall(pix_data);
                if p_value > 0.05
                    mk = nan;
                elseif p_value < 0.05
                    mk = z;
                end
         % Fit a linear model (polynomial of degree 1)
         p = polyfit(SDATE, pix_data, 1);
         slope = p(1); % Extract the slope and intercept from the polynomial coefficients
         residuals = pix_data - polyval(p, SDATE); 
         residualSE = sqrt(sum(residuals.^2) / (length(pix_data) - 2));
         t_statistic = slope / (residualSE / sqrt(sum((SDATE - mean(SDATE)).^2)));
         p_value_lst_sqr = 2 * (1 - tcdf(abs(t_statistic), length(pix_data) - 2));
                if p_value_lst_sqr > 0.05
                    lst_sqr = nan;
                elseif p_value_lst_sqr < 0.05
                    lst_sqr = slope;
                end
            end
        mk_map(j,k) = mk;
        lst_sqr_map(j,k) = lst_sqr;
    end
end

addpath('I:\My Drive\GDrive\OCED_AOML\Omics\sat_analysis');
load('cmap');
figure();
    m_proj('equidistant cylindrical','lon',LONLIMS,'lat',LATLIMS);
    m_pcolor(lon1,lat1,mk_map);
    shading flat;
    colormap(cm);
%     hold on
%     [C, hContour] = m_contour(lon1, lat1, mk_map, [2.5,2.5], 'LineColor', 'r', 'LineWidth', 2);
%     hold off
    m_gshhs_i('patch',[0.5 0.5 0.5],'edgecolor','none');
    m_grid_v2('box','fancy','fontsize',12,'fontweight','bold');
    h = colorbar;
    h.FontSize = 12; % Adjust the size as needed
    scale = [-4 4];
    caxis(scale);
    title('Seasonal Mann-Kendall trend for Chl-a (p < 0.05)', 'FontSize',16);
%     clabel(C, hContour, 'LabelSpacing', 100); % Adjust LabelSpacing as needed
    ylabel(h, 'S', 'FontSize', 14, 'Rotation', 0); % Change 'Your Colorbar Label' to your actual label

% cm3 = cmocean('thermal');    
figure();
    m_proj('equidistant cylindrical','lon',LONLIMS,'lat',LATLIMS);
    m_pcolor(lon1,lat1,lst_sqr_map);
    shading flat;
    colormap(cm);
    hold on
    [C, hContour] = m_contour(lon1, lat1, lst_sqr_map, [-0.01,0.01], 'LineColor', 'k', 'LineWidth', 2);
    hold off
    m_gshhs_i('patch',[0.5 0.5 0.5],'edgecolor','none');
    m_grid_v2('box','fancy','fontsize',12,'fontweight','bold');
    h = colorbar;
    h.FontSize = 12; % Adjust the size as needed
    scale = [-0.03 0.03];
    caxis(scale);
    title('Monthly Chl-a linear trend (p < 0.05)', 'FontSize',16);
    clabel(C, hContour, 'LabelSpacing', 100); % Adjust LabelSpacing as needed
    ylabel(h, 'Slope', 'FontSize', 14, 'Rotation', 90)





