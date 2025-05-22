## Copyright (C) 2016 Nir Krakauer
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

##predict values using specified forecast methods
##
##outputs a cell array of files with the predictions, one file per method

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>



function files = sefo_predict (predict_year, predict_month, lag, methods, opts)

  warning ("off", "Octave:broadcast", "global");

  data_dir = getfield (opts, 'data_dir');
  re_download = getfield (opts, 're_download');
  re_process = getfield (opts, 're_process');  
  var_name = getfield (opts, 'var_name');
  fcst_lons = getfield (opts, 'fcst_lons');
  fcst_lats = getfield (opts, 'fcst_lats');
  fcst_start_year = getfield (opts, 'fcst_start_year');
  obs_start_year = getfield (opts, 'obs_start_year');
  obs_end_year = getfield (opts, 'obs_end_year');
  clim_years = getfield (opts, 'clim_years');
  obs_dataset = getfield (opts, 'obs_dataset');
  p_spline = getfield (opts, 'p_spline');
  dof_spline = getfield (opts, 'dof_spline');
  lamda = getfield (opts, 'lamda');    

  n_months_combine = numel (lag);

  if iscell (methods)
    methods_use = methods;
  else
    methods_use = {methods};
  endif

  n_methods = numel(methods_use);
  mus = sigmas = nus = [];

  for method_index = 1:n_methods

    method = methods_use{method_index};

    file = [data_dir 'predict_' var_name '_' method '_' obs_dataset];
    for m = 1:n_months_combine
      file = [file '_' num2str(predict_year(m)) '_' num2str(predict_month(m), '%2.2i') '_lag_' num2str(lag(m), '%2.2i')];
    endfor
    file = [file '.mat'];

    run_function = true;
    if exist(file, "file") && !re_process
      if strcmp (method, "multimethod") #rerun unless the methods used were the same
        load (file, "methods")
        if isequal (methods, methods_use)
          run_function = false;
        endif
      else
        run_function = false;
      endif    
    endif

    function_already_run = false;

    if run_function

      if ~function_already_run
        if n_months_combine > 1
          if (numel(predict_year) != n_months_combine) || (numel(predict_month) != n_months_combine)
            error ('predict_year, predict_month, lag should have the same numbers of elements')
          endif
          n_days_per_month = days_per_month (predict_month, predict_year);
          w = n_days_per_month / sum(n_days_per_month);
          fcsts_mean_combine = fcsts_mean_trend_combine = obs_combine = obs_trend_combine = 0;
          for m = 1:n_months_combine
            file_obs = sefo_obs_assemble (predict_year(m), predict_month(m), opts);
            load (file_obs)
            if m == 1
              ny_obs = size (obs, 2);
            endif
            obs_combine += w(m)*obs(:, (end-ny_obs+1):end);
            obs_trend_combine += w(m)*obs_trend(:, (end-ny_obs):end);
            file_fcst = sefo_fcst_assemble (predict_year(m), predict_month(m), lag(m), opts);
            load (file_fcst)
            if m == 1
              ny_fcsts = size (fcsts_mean, 2);
            endif
            fcsts_mean_combine += w(m)*fcsts_mean(:, (end-ny_fcsts+1):end, :);
            fcsts_mean_trend_combine += w(m)*fcsts_mean_trend(:, (end-ny_fcsts+1):end, :);
          endfor
          obs = obs_combine; clear obs_combine
          obs_trend = obs_trend_combine; clear obs_trend_combine
          fcsts_mean = fcsts_mean_combine; clear fcsts_mean_combine
          fcsts_mean_trend = fcsts_mean_trend_combine; clear fcsts_mean_trend_combine
        else
          file_obs = sefo_obs_assemble (predict_year, predict_month, opts);
          load (file_obs)
          ny_obs = size (obs, 2);
          file_fcst = sefo_fcst_assemble (predict_year, predict_month, lag, opts);
          load (file_fcst)
        endif
        
        #areas represented by each grid cell
        nlats = numel(fcst_lats);
        nlons = numel(fcst_lons);
        nmme_lat_lims = nan(nlats+1, 1);
        nmme_lat_lims(2:(end-1)) = (fcst_lats(1:(end-1)) + fcst_lats(2:end))/2;
        nmme_lat_lims(1) = -90;
        nmme_lat_lims(end) = 90;
        if strcmp (obs_dataset, 'bestland') #all grid points represent equal areas
          area_wgts = ones (nlats, nlons);
        else
          area_wgts = diff(sind(nmme_lat_lims)) * ones(1, nlons);
        endif
        #only consider grid cells with available historical observations
        if exist("obs_inds", "var") && ~isempty(obs_inds)
          fcsts_mean = fcsts_mean(obs_inds, :, :);
          fcsts_mean_trend = fcsts_mean_trend(obs_inds, :, :);
          area_wgts = area_wgts(obs_inds);
        else
          obs_inds = [];
        endif
            
        area_wgts = area_wgts(:);
        area_wgts /= sum(area_wgts);    
        np = size (fcsts_mean, 1);
        #obs_global = area_wgts' * obs;
        #fcsts_mean_global = area_wgts' * mean(fcsts_mean, 3);
        
        #it is assumed that fcst_start_year >= obs_start_year and that forecasts are available through predict_year
        if (obs_end_year - obs_start_year + 1) != ny_obs
          obs_end_year = obs_start_year + ny_obs - 1;
          #warning (['sefo_predict: input obs_end_year not consistent with saved sefo_obs_assemble output, using ' num2str(obs_end_year)])
        endif
        ny_fcsts = size (fcsts_mean, 2) - max(1, predict_year(end)-obs_end_year); #years with both forecasts and observations
        n_gcms = size (fcsts_mean, 3);

        obs_offset = fcst_start_year - obs_start_year;
        
        fcsts_mean = double (fcsts_mean);
        fcsts_mean_trend = double (fcsts_mean_trend);

        #for precipitation, transform observations and forecasts to approximate standard normal distributions
        switch var_name
          case {'prec'}
            obs_trend = sefo_normalization (obs, struct ("data_to_transform", obs_trend));
            obs = sefo_normalization (obs);
            for f = 1:n_gcms
              fcsts_mean_trend(:, :, f) = sefo_normalization (fcsts_mean(:, :, f), struct ("data_to_transform", fcsts_mean_trend(:, :, f)));
              fcsts_mean(:, :, f) = sefo_normalization (fcsts_mean(:, :, f));
            endfor
        endswitch

        options.w = area_wgts;
        options.dof_spline = dof_spline;
        options.lamda = lamda;
        options.fcst_trend = fcsts_mean_trend;
        options.obs_trend = obs_trend;
        options.t_clim = clim_years - obs_start_year + 1;
        if strcmp (method, "MA")
          options.timescale = 30;
        elseif length (method) >= 2 && strcmp (method(1:2), "EW")
          options.timescale = 15;
        endif

        function_already_run = true;
      endif

      options.mus = mus;
      options.sigmas = sigmas;
      options.nus = nus;

      [mu, sigma, nu] = sefo_prediction (fcsts_mean, obs, method, options);

      methods = methods_use;
      save (file, "mu", "sigma", "nu", "methods", "area_wgts")
    else
      load (file, "mu", "sigma", "nu")
    endif

    if method_index == 1
      np = numel(mu);
      mus = sigmas = nus = nan (np, n_methods);
    endif
    mus(:, method_index) = mu;
    sigmas(:, method_index) = sigma;
    nus(:, method_index) = nu;

    files{method_index} = file;

  endfor

  
endfunction

%!demo
%! if !exist("data_dir", "var")
%!   data_dir = [pwd '/'];
%! endif
%! methods = {'linear', 'rawmodel'};
%! opts = struct ('data_dir', data_dir, ...
%!                 're_download', false, ...
%!                 're_process', false, ...
%!                 'var_name', 'tref', ...
%!                 'fcst_lons', (0:359)', ...
%!                 'fcst_lats', (-90:90)', ... 
%!                 'fcst_start_year', 1982, ...
%!                 'obs_start_year', 1957, ...  
%!                 'obs_end_year', 2015, ...                  
%!                 'obs_dataset', 'best', ...              
%!                 'p_spline', 3E-4, ...
%!                 'dof_spline', 2.5, ...
%!                 'lamda', 10, ...
%!                 'predict_adj', true, ...
%!                 'predict_adj_years', 1994:1998, ...
%!                 'predict_adj_months', 1:12, ...                
%!                 'methods', methods, ...
%!                 'clim_years', 1981:2010, ...   
%!                 'gcms', {{'CMC1-CanCM3', 'CMC2-CanCM4', 'GFDL-CM2p1-aer04', 'COLA-RSMAS-CCSM3', 'NASA-GMAO-062012'}} ...
%!                 );
%! predict_year = 2015; predict_month = 1; lag = 1;
%! files = sefo_predict (predict_year, predict_month, lag, methods, opts)
