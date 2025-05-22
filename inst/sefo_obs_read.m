## Copyright (C) 2015-2016 Nir Krakauer
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

##read, fill in missing, and regrid observed values for one month

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>



function file = sefo_obs_read (year, month, opts)

  data_dir = getfield (opts, 'data_dir');
  re_download = getfield (opts, 're_download');
  re_process = getfield (opts, 're_process');  
  var_name = getfield (opts, 'var_name');
  obs_start_year = getfield (opts, 'obs_start_year');
  obs_dataset = getfield (opts, 'obs_dataset');
  fcst_lons = getfield (opts, 'fcst_lons');
  fcst_lats = getfield (opts, 'fcst_lats');
  
  file = [data_dir 'obs_regrid_nmme_' var_name '_' obs_dataset '_' num2str(year) '_' num2str(month, '%2.2i') '.mat'];

  if ~exist(file, "file") || re_process

    nlats = numel(fcst_lats);
    nlons = numel(fcst_lons);
    fcst_lon_lims = nan(nlons+1, 1);
    fcst_lon_lims(2:(end-1)) = (fcst_lons(1:(end-1)) + fcst_lons(2:end))/2;
    fcst_lon_lims(1) = (fcst_lons(1) + fcst_lons(end))/2 - 180;
    fcst_lon_lims(end) = fcst_lon_lims(1) + 360;
    fcst_lat_lims = nan(nlats+1, 1);
    fcst_lat_lims(2:(end-1)) = (fcst_lats(1:(end-1)) + fcst_lats(2:end))/2;
    fcst_lat_lims(1) = -90;
    fcst_lat_lims(end) = 90;
    fcst_lat_lims = sind (fcst_lat_lims);

    switch obs_dataset
      case 'cams'

      switch var_name
        case 'tref'
          obs_file_name = 't.long';
        otherwise
          error('invalid var_name')
      endswitch

      obs_file = [data_dir obs_file_name];

      if ~exist(obs_file, "file") || re_download
        data_url = ['ftp://ftp.cpc.ncep.noaa.gov/wd51yf/global_monthly/gridded_binary/' obs_file_name];
        urlwrite(data_url, obs_file);
      endif

      fid = fopen (obs_file, 'r', 'ieee-be');

      nlon = 720;
      nlat = 360;
      np = nlon*nlat;

      lons = ((1:nlon) - 0.5) / 2;
      lats = ((1:nlat) - 180.5) / 2;

      offset = (4*np + 8)*(12*(year - obs_start_year) + (month - 1)) + 4; #in bytes; records preceded and followed by int32
      fseek (fid, offset);

      arr = fread (fid, [nlon, nlat], 'single', 0, 'ieee-be');

      cams_nan = -999;
      arr(arr == cams_nan) = NaN;

      fclose (fid);

      obs_lons = lons;
      obs_lats = lats;
      nobs_lats = numel(obs_lats);
      nobs_lons = numel(obs_lons);
      obs_lon_lims = nan(nobs_lons+1, 1);
      obs_lon_lims(2:(end-1)) = (obs_lons(1:(end-1)) + obs_lons(2:end))/2;
      obs_lon_lims(1) = (obs_lons(1) + obs_lons(end))/2 - 180;
      obs_lon_lims(end) = obs_lon_lims(1) + 360;
      obs_lat_lims = nan(nobs_lats+1, 1);
      obs_lat_lims(2:(end-1)) = (obs_lats(1:(end-1)) + obs_lats(2:end))/2;
      obs_lat_lims(1) = -90;
      obs_lat_lims(end) = 90;
      obs_lat_lims = sind (obs_lat_lims);
   
      #fill in missing values
      arr_lat_mean = nan (nobs_lats, 1);
      for j = 1:nobs_lats
        non_missing_values = isfinite (arr(:, j));
        if any (non_missing_values)
          arr_lat_mean(j) = mean(arr(non_missing_values, j));
        endif
      endfor
      missing_values = isnan(arr_lat_mean);
      if any(missing_values)
        arr_lat_mean(missing_values) = interp1(lats(~missing_values), arr_lat_mean(~missing_values), lats(missing_values), 'nearest', 'extrap');
      endif
      
      arr_aug = [arr' arr' arr'];
      valid_aug = isfinite(arr_aug);
      for j = 1:nobs_lats
        missing_values = isnan(arr_aug(j, :));
        if all (missing_values)
          arr_aug(j, :) = arr_lat_mean(j);
        elseif any (missing_values)
          arr_aug(j, missing_values) = interp1((1:(3*nobs_lons))'(~missing_values), arr_aug(j, ~missing_values), (1:(3*nobs_lons))'(missing_values), 'nearest', 'extrap');
        endif
      endfor

      coords = {obs_lat_lims, [obs_lon_lims(1:(end-1))'-360 obs_lon_lims' obs_lon_lims(2:end)'+360]'};
      new_coords = {fcst_lat_lims, fcst_lon_lims};
      arr_regrid = regrid (coords, arr_aug, new_coords);
   
      arr_wgt = regrid (coords, valid_aug, new_coords); #fraction of each NMME gridcell area that actually has observations
   
      arr_regrid(arr_wgt < 0.5) = NaN;
 
    case "ncep"
      switch var_name
        case 'tref'
          var_name_ncep = 'air';  
          var_name_ncep_filename = 'air.2m';
          ncep_dir = 'surface_gauss';
        otherwise
          error('invalid var_name')
      endswitch
    
      filename =  [data_dir "lsmask.19294.nc"]; #land-sea mask
      if ~exist (filename, "file")
        urlwrite("ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface_gauss/lsmask.19294.nc", filename);
      endif   
      #lf = -ncread (filename, "lsmask"); #originally land:-1, ocean:0, no fractional values
      obs_lons = ncread (filename, "lon");
      obs_lats = flipud (ncread (filename, "lat"));   
      nobs_lats = numel(obs_lats);
      nobs_lons = numel(obs_lons);
      obs_lon_lims = nan(nobs_lons+1, 1);
      obs_lon_lims(2:(end-1)) = (obs_lons(1:(end-1)) + obs_lons(2:end))/2;
      obs_lon_lims(1) = (obs_lons(1) + obs_lons(end))/2 - 180;
      obs_lon_lims(end) = obs_lon_lims(1) + 360;
      obs_lat_lims = nan(nobs_lats+1, 1);
      obs_lat_lims(2:(end-1)) = (obs_lats(1:(end-1)) + obs_lats(2:end))/2;
      obs_lat_lims(1) = -90;
      obs_lat_lims(end) = 90;
      obs_lat_lims = sind (obs_lat_lims);
      ncep_filename = [var_name_ncep_filename '.mon.mean.nc'];
      obs_file = [data_dir ncep_filename];
      if ~exist (obs_file, "file")
        urlwrite(["ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/" ncep_dir "/" ncep_filename], obs_file);
      endif
      start_year = 1948;
      m_index = 12*(year - start_year) + month;
      arr = squeeze(ncread(obs_file, var_name_ncep, [1 1 m_index], [Inf Inf 1]))'; #original is 192*94*809
      arr = flipud(arr); #ascending latitudes

      arr_aug = [arr arr arr];
      coords = {obs_lat_lims, [obs_lon_lims(1:(end-1))'-360 obs_lon_lims' obs_lon_lims(2:end)'+360]'};
      new_coords = {fcst_lat_lims, fcst_lon_lims};
      arr_regrid = regrid (coords, arr_aug, new_coords);  
      
    case "best"
      source_file_name = ['Land_and_Ocean_LatLong1.nc'];
      if ~exist([data_dir source_file_name], 'file') || re_download
        data_url = ['http://berkeleyearth.lbl.gov/auto/Global/Gridded/' source_file_name];
        urlwrite(data_url, [data_dir source_file_name])  #download and save netcdf file 
      endif
      source_file = [data_dir source_file_name];
      lon_obs = ncread (source_file, 'longitude'); #-179.5:179.5
      lat_obs = ncread (source_file, 'latitude'); #-89.5:89.5
      t_obs = ncread (source_file, 'time');
      nx_obs = numel (lon_obs); #360
      ny_obs = numel (lat_obs); #180
      start_year = floor (t_obs(1)); #1850 [assumes it starts at the beginning of a year]
      mon = month + 12*(year - start_year);   
      arr = ncread (source_file, 'temperature', [1 1 mon], [nx_obs ny_obs 1]); #size: nx_obs*ny_obs
      
      #fill in any missing values with the mean for that latitude; interpolate to get values for any completely missing latitudes
      if any(isnan(arr(:)))
        lat_mean = mean(arr, 1);
        if any(isnan(lat_mean))
          s = isfinite (lat_mean);
          lat_mean(~s) = interp1 (lat_obs(s), lat_mean(s), lat_obs(~s), 'nearest', 'extrap');
        endif
        s = isnan (arr);
        arr(s) = (ones(nx_obs, 1) * lat_mean)(s);
      endif
      
      #regrid to match forecast (NMME)
      arr = arr';
      arr_aug = [arr arr arr];
      obs_lat_lims = sind([-90 lat_obs(1:(end-1))'+diff(lat_obs')/2 90])';
      obs_lon_lims = [lon_obs(end)-360 lon_obs' lon_obs(1)+360];
      obs_lon_lims =  obs_lon_lims(1:(end-1))'+diff(obs_lon_lims')/2;
      obs_lon_lims = [obs_lon_lims(1:(end-1))'-360 obs_lon_lims' obs_lon_lims(2:end)'+360]';
      coords = {obs_lat_lims, obs_lon_lims};
      new_coords = {fcst_lat_lims, fcst_lon_lims};
      arr_regrid = regrid (coords, arr_aug, new_coords); 

    case "bestland"
      #observations (5498 land cells only)
      source_file_name = ['Complete_TAVG_EqualArea.nc'];
      if ~exist([data_dir source_file_name], 'file') || re_download
        data_url = ['http://berkeleyearth.lbl.gov/auto/Global/Gridded/' source_file_name];
        urlwrite(data_url, [data_dir source_file_name])  #download and save netcdf file 
      endif
      source_file = [data_dir source_file_name];
      lon_obs = ncread (source_file, 'longitude');
      lat_obs = ncread (source_file, 'latitude');
      t_obs = ncread (source_file, 'time');
      n_obs = numel(lon_obs); #5498
      lon_obs(lon_obs < 0) += 360; #transform to 0-360 range
      start_year = floor(t_obs(1)); #1750
      mon = month + 12*(year - start_year);
      arr = ncread (source_file, 'temperature', [1 mon], [n_obs 1]);
  
      #closest forecast grid points
      n_fcst_lats = numel(fcst_lats);
      n_fcst_lons = numel(fcst_lons);
      fcst_inds = sub2ind([n_fcst_lats n_fcst_lons], round(interp1(fcst_lats, 1:n_fcst_lats, lat_obs, 'extrap')), round(interp1(fcst_lons, 1:n_fcst_lons, lon_obs, 'extrap')));
      
      arr_regrid = nan (n_fcst_lats, n_fcst_lons);
      arr_regrid(fcst_inds) = arr;
            
    case "gpcp" #precipitation
      source_file_name = ['precip.mon.mean.nc'];
      if ~exist([data_dir source_file_name], 'file') || re_download
        data_url = ['ftp://ftp.cdc.noaa.gov/Datasets/gpcp/' source_file_name];
        urlwrite(data_url, [data_dir source_file_name])  #download and save netcdf file 
      endif
      source_file = [data_dir source_file_name];
      lon_obs = ncread (source_file, 'lon');
      lat_obs = ncread (source_file, 'lat');
      nx_obs = numel (lon_obs); #144
      ny_obs = numel (lat_obs); #72
      start_year = 1979;
      mon = month + 12*(year - start_year);
      arr = ncread (source_file, 'precip', [1 1 mon], [nx_obs ny_obs 1]); #size: nx_obs*ny_obs
      
      #fill in any missing values (apparently there are none) with the mean for that latitude; interpolate to get values for any completely missing latitudes
      arr(arr < 0) = NaN;
      if any(isnan(arr(:)))
        lat_mean = mean(arr, 1);
        if any(isnan(lat_mean))
          s = isfinite (lat_mean);
          lat_mean(~s) = interp1 (lat_obs(s), lat_mean(s), lat_obs(~s), 'nearest', 'extrap');
        endif
        s = isnan (arr);
        arr(s) = (ones(nx_obs, 1) * lat_mean)(s);
      endif

      #regrid to match forecast (NMME)
      arr = arr';
      arr_aug = [arr arr arr];
      obs_lat_lims = sind([-90 lat_obs(1:(end-1))'+diff(lat_obs')/2 90])';
      obs_lon_lims = [lon_obs(end)-360 lon_obs' lon_obs(1)+360];
      obs_lon_lims =  obs_lon_lims(1:(end-1))'+diff(obs_lon_lims')/2;
      obs_lon_lims = [obs_lon_lims(1:(end-1))'-360 obs_lon_lims' obs_lon_lims(2:end)'+360]';
      coords = {obs_lat_lims, obs_lon_lims};
      new_coords = {fcst_lat_lims, fcst_lon_lims};
      arr_regrid = regrid (coords, arr_aug, new_coords);
     
    case 'hadsst' #5-degree resolution
      switch var_name
        case 'sst'
        otherwise
          error('invalid var_name for this data source')
      endswitch    
      vers = '3.1.1.0';
      source_file_name = ['HadSST.' vers '.median'];
      if ~exist([data_dir source_file_name '.nc'], 'file') || re_download
        data_url = ['https://www.metoffice.gov.uk/hadobs/hadsst3/data/HadSST.' vers '/netcdf/' source_file_name '_netcdf.zip'];
        urlwrite(data_url, [data_dir source_file_name '.zip'])  #download and save netcdf file 
        system (['unzip ' data_dir source_file_name '.zip -d ' data_dir]);
      endif          

      source_file = [data_dir source_file_name '.nc'];
      lon_obs = ncread (source_file, 'longitude');
      lat_obs = ncread (source_file, 'latitude');
      nx_obs = numel (lon_obs); #72
      ny_obs = numel (lat_obs); #36
      start_year = 1850;
      mon = month + 12*(year - start_year);
      arr = ncread (source_file, 'sst', [1 1 mon], [nx_obs ny_obs 1]); #size: nx_obs*ny_obs

      #fill in any missing values with the mean for that latitude; interpolate to get values for any completely missing latitudes
      if any(isnan(arr(:)))
        lat_mean = mean(arr, 1); #assumes that the nan package's mean is active
        if any(isnan(lat_mean))
          s = isfinite (lat_mean);
          lat_mean(~s) = interp1 (lat_obs(s), lat_mean(s), lat_obs(~s), 'nearest', 'extrap');
        endif
        s = isnan (arr);
        arr(s) = (ones(nx_obs, 1) * lat_mean)(s);
      endif

      #flip latitudes to be in increasing order
      lat_obs = flip (lat_obs, 1);
      arr = flip (arr, 2);

      #regrid to match forecast (NMME)
      arr = arr';
      arr_aug = [arr arr arr];
      obs_lat_lims = sind([-90 lat_obs(1:(end-1))'+diff(lat_obs')/2 90])';
      obs_lon_lims = [lon_obs(end)-360 lon_obs' lon_obs(1)+360];
      obs_lon_lims =  obs_lon_lims(1:(end-1))'+diff(obs_lon_lims')/2;
      obs_lon_lims = [obs_lon_lims(1:(end-1))'-360 obs_lon_lims' obs_lon_lims(2:end)'+360]';
      coords = {obs_lat_lims, obs_lon_lims};
      new_coords = {fcst_lat_lims, fcst_lon_lims};
      arr_regrid = regrid (coords, arr_aug, new_coords);
      
    otherwise
      error ([obs_dataset " not recognized"])
    endswitch
 
    save (file, "arr_regrid")  

  endif
endfunction

#colormap('default'); imagesc_withnan(lons, lats, arr'); axis xy; colorbar
#colormap('default'); imagesc_withnan(lons, lats, arr_aug(:, 721:1440)); axis xy; colorbar
#colormap('default'); imagesc_withnan(fcst_lons, fcst_lats, arr_wgt); axis xy; colorbar
#colormap('default'); imagesc_withnan(fcst_lons, fcst_lats, arr_regrid, [-2 2]); axis xy; colorbar; pcontinents2 
  #cf. with https://www.ncdc.noaa.gov/sotc/global/201401

%!demo
%! if !exist("data_dir", "var")
%!   data_dir = [pwd '/'];
%! endif
%! fcst_lons = (0:359)';
%! fcst_lats = (-90:90)';
%! opts = struct ('data_dir', data_dir, ...
%!                 're_download', false, ...
%!                 're_process', false, ...
%!                 'var_name', 'tref', ...
%!                 'fcst_lons', fcst_lons, ...
%!                 'fcst_lats', fcst_lats, ... 
%!                 'obs_start_year', 1957, ...  
%!                 'obs_dataset', 'best' ...                 
%!                 );
%! year = 2015; month = 12;
%! file = sefo_obs_read (year, month, opts);
%! load (file)
%! colormap("default"); imagesc_withnan(fcst_lons, fcst_lats, arr_regrid, [-5 5]); axis xy; pcontinents2
%! # map observed temperatures for an example month (Jan 2014) using the default colormap
%! # cf. with https://www.ncdc.noaa.gov/sotc/global/201401

