## Copyright (C) 2015 Nir Krakauer
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

##read 1 month of NMME ensemble forecast from a given GCM at a given lag, saving the ensemble average 

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>



%{
#example:
if !exist("data_dir", "var")
  data_dir = [pwd '/'];
endif
pkg load netcdf
opts = struct ('data_dir', data_dir, ...
                're_download', false, ...
                're_process', false, ...
                'var_name', 'tref', ...
                'fcst_lons', (0:359)', ...
                'fcst_lats', (-90:90)' ...                
                );
year = 2013; month = 11; lag = 0; gcm = 'COLA-RSMAS-CCSM3';
file = sefo_fcst_read (year, month, lag, gcm, opts);
load (file);
v = vals_mean;
a = [min(v(:)) mean(v(:)) max(v(:))]
#ll = quantile(abs(v(:)), 0.95); lims = [-ll ll];
lims = [a(1) a(end)];
colormap("default"); imagesc (lons, lats, v', lims); axis xy; pcontinents2; colorbar
%}

function file = sefo_fcst_read (year, month, lag, gcm, opts)

  data_dir = getfield (opts, 'data_dir');
  re_download = getfield (opts, 're_download');
  re_process = getfield (opts, 're_process');  
  var_name = getfield (opts, 'var_name');
  fcst_lons = getfield (opts, 'fcst_lons');
  fcst_lats = getfield (opts, 'fcst_lats');


  file = [data_dir 'nmme_' var_name '_' gcm '_' num2str(year) '_' num2str(month, '%2.2i') '_lag_' num2str(lag, '%2.2i') '.mat'];
  
  if ~exist(file, "file") || re_process

    if (strcmp(gcm, 'NCEP-CFSv2') && (lag > 9)) || (strcmp(gcm, 'NASA-GMAO-062012') && (lag > 8)) #these models do not have forecasts for the longer lags, so the IRI website misleadingly returns the forecast for the same forecast month with the longest available lag instead
        vals_mean = lats = lons = [];
        nr = 0;
        warning (['empty ' file])
        save (file, "vals_mean", "vals_sd", "lats", "lons", "nr")
        return
    endif

    source_file = [data_dir 'nmme_' var_name '_' gcm '_' num2str(year) '_' num2str(month, '%2.2i') '_lag_' num2str(lag, '%2.2i') '.nc'];

    if ~exist(source_file, "file") || re_download

      #when the forecast was made (i.e. lag months before the target month)
      fcst_month = month - lag;
      fcst_year = year;
      if fcst_month < 1
        fcst_year -= 1;
        fcst_month += 12;
      endif

      month_names = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
    
      if strcmp (gcm(1:3), 'CMC') 
        if fcst_year < 2011
          extra_string = '/.HINDCAST';
        else
          extra_string = '/.FORECAST';
        endif
      elseif strcmp (gcm, 'NCEP-CFSv2')
        if fcst_year < 2011 || (fcst_year == 2011 && fcst_month < 4)
          extra_string = '/.HINDCAST/.PENTAD_SAMPLES';
        else
          extra_string = '/.FORECAST/.PENTAD_SAMPLES';
        endif
      else
        extra_string = '';
      endif

      data_url = ['http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/.' gcm extra_string '/.MONTHLY/.' var_name '/L/%28' num2str(lag) '%29%28' num2str(lag+1) '%29RANGEEDGES/S/%280000%201%20' month_names{fcst_month} '%20' num2str(fcst_year) '%29%280000%2015%20' month_names{fcst_month} '%20' num2str(fcst_year) '%29RANGEEDGES/data.cdf'];
      urlwrite(data_url, source_file);  %download and save netcdf file
     
    endif
      
      try     
        vals = ncread (source_file, var_name);
        vals = squeeze (vals); #size:nlons x nlats x nruns; type:single
        lons = ncread (source_file, 'X');
        lats = ncread (source_file, 'Y');    
      catch
        warning(['unable to read ' source_file])
        vals = [];
      end_try_catch
        
       if !isempty (vals) 
        if any (vals(:) > 1E20)
          warning(['values above 1E20 in ' source_file])
          vals(vals > 1E20) = NaN;
        endif 
      
        if any (vals(:) < 0) #can't have negative absolute temperature or precipitation
          #warning(['values below 0 in ' source_file])
          vals(vals < 0) = NaN;
        endif

        if strcmp (gcm, 'NCEP-CFSv2') && strcmp (var_name, 'tref') && any (vals(:) == 0) #NCEP-CFSv2 has spurious zeros for the lag-10 forecast of 2013/03
          warning(['values of 0 in ' source_file])
          vals(vals == 0) = NaN;
        endif

        if any (lons != fcst_lons)
          warning(['different longitudes in ' source_file])
        endif
        
        if any (lats != fcst_lats)
          warning(['different latitudes in ' source_file])
        endif        

        if any (((sum(isfinite(vals), 3)(:))) == 0)
          if any (isfinite(vals(:))) #assume that model is only missing the first and last (polar) latitudes, so fill those in [this is the case with GFDL]
            vals(:, 1, :) = vals(:, 2, :);
            vals(:, end, :) = vals(:, end-1, :);
            if any (((sum(isfinite(vals), 3)(:))) == 0) #check that we were successful at filling in
              warning(['unsuccessful at filling in ' source_file])
            endif
          else
            vals = [];
          endif
        endif
      endif     

      if isempty(vals)
        vals_mean = lats = lons = [];
        nr = 0;
        warning (['empty ' file])
      else
        vals_mean = squeeze (mean (vals, 3));
        vals_sd = squeeze (std (vals, [], 3));
        nr = sum (isfinite(vals(1, 1, :))(:)); #number of ensemble members
        vals_mean = single (vals_mean);
        vals_sd = single (vals_sd);
      endif
      
      save (file, "vals_mean", "vals_sd", "lats", "lons", "nr")
  endif

endfunction

%!demo
%! if !exist("data_dir", "var")
%!   data_dir = [pwd '/'];
%! endif
%! opts = struct ('data_dir', data_dir, ...
%!                 're_download', true, ...
%!                 're_process', true, ...
%!                 'var_name', 'tref', ...
%!                 'fcst_lons', (0:359)', ...
%!                 'fcst_lats', (-90:90)' ...                
%!                 );
%! year = 2008; month = 9; lag = 11; gcm = 'GFDL-CM2p5-FLOR-A06';
%! file = sefo_fcst_read (year, month, lag, gcm, opts);
%! load (file)
%! colormap("default"); imagesc_withnan (lons, lats, vals_mean', [min(vals_mean(:)) max(vals_mean(:))]); axis xy; pcontinents2; colorbar  
%! # map predicted temperatures for an example month
