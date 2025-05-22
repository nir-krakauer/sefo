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

##compare prediction with observations

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>


function file = sefo_verify (years, months, lag, opts)

  data_dir = getfield (opts, 'data_dir');
  re_download = getfield (opts, 're_download');
  re_process = getfield (opts, 're_process');  
  var_name = getfield (opts, 'var_name');
  fcst_lons = getfield (opts, 'fcst_lons');
  fcst_lats = getfield (opts, 'fcst_lats');
  fcst_start_year = getfield (opts, 'fcst_start_year');
  obs_start_year = getfield (opts, 'obs_start_year');   
  obs_dataset = getfield (opts, 'obs_dataset');
  methods = getfield (opts, 'methods');
  predict_adj = getfield (opts, 'predict_adj');
  predict_adj_years = getfield (opts, 'predict_adj_years');
  predict_adj_months = getfield (opts, 'predict_adj_months');
  clim_years = getfield (opts, 'clim_years');

  n_months_combine = numel (lag);
  n_methods = numel(methods);
  n_years = numel(years);
  n_months = numel(months);
  if size(years) == size(months)
    nt = n_years;
  elseif n_years > 1 && n_months == 1
    months = months * ones(size(years));
    nt = n_years;
  elseif n_months > 1 && n_years == 1 
    years = years * ones(size(months));
    nt = n_months;
  else
    years = ones(n_months, 1) * years(:)';
    months = months(:) * ones(1, n_years);
    nt = n_years*n_months;
  endif
  if nt == 1
    dates = [num2str(years(1)) '_' num2str(months(1), '%2.2i')];
  else
    dates = [num2str(years(1)) '_' num2str(months(1), '%2.2i') '_' num2str(years(end)) '_' num2str(months(end), '%2.2i')];
  endif

  if isfield(opts, 'subset_verify') && !isempty(opts.subset_verify)
    subset_tag = ['_' opts.subset_verify];
  else
    subset_tag = '';
  endif

  if predict_adj
    adj_tag = '_adj';
    factors_adj = nan(n_methods, 2);
  else
    adj_tag = '';
    factors_adj = [];
  endif

  rank_calc = !isempty (clim_years);
  if rank_calc
    rank_tag = '_rnk';
  else
    rank_tag = '';
  endif

  file = [data_dir 'verify_' var_name '_' obs_dataset '_' dates '_lag' num2str(lag, '_%2.2i') adj_tag rank_tag subset_tag '.mat'];

  disp ([datestr(clock) "\n" 'sefo_verify: initialized'])

  if ~exist(file, "file") || re_process

    obss = area_wgtss = [];
    for t = 1:nt
      year = years(t);
      month = months(t);
      predict_month = month + (0:(n_months_combine-1))';
      predict_year = year * ones(n_months_combine, 1);
      inds_next_year = find(predict_month > 12);
      predict_year(inds_next_year) = predict_year(inds_next_year) + 1;
      predict_month(inds_next_year) = predict_month(inds_next_year) - 12;
      n_days_per_month = days_per_month (predict_month, predict_year);
      w = n_days_per_month / sum(n_days_per_month);
      obs = 0;
      for m = 1:n_months_combine
        file_obs = sefo_obs_read (predict_year(m), predict_month(m), opts);
        load (file_obs)
        if ~exist("obs_inds", "var") || isempty(obs_inds)
          obs_inds = isfinite (arr_regrid);
        endif
        obs += w(m)*arr_regrid(obs_inds);
        clear arr_regrid
      endfor
      if t == 1
        np = numel (obs);
        obss = nan (np, nt, "single");
      endif
      obss(:, t) = obs;

      #for precipitation, transform the verifying observations to approximate standard normal distribution based on past observations
      switch var_name
          case {'prec'}
            obs_hist = 0;
            for m = 1:n_months_combine
              file_obs = sefo_obs_assemble (predict_year(m), predict_month(m), opts);
              load (file_obs, "obs")
              if m == 1
                ny_obs = size (obs, 2);
              endif
              obs_hist += w(m)*obs(:, (end-ny_obs+1):end);
            endfor
            obss(:, t) = sefo_normalization (obs_hist, struct ("data_to_transform", obss(:, t)));
            disp ([datestr(clock) "\n" 'sefo_verify: normalized observations'])
      endswitch
      if rank_calc
      
        #get climatology
        clim_opts = opts;
        clim_opts.obs_start_year = min(clim_years);
        clim_obs_file = sefo_obs_assemble (max(clim_years)+1, month, clim_opts); ##TODO: get to work for n_months_combine > 1
        load (clim_obs_file, "obs");
        
        switch var_name
            case {'prec'}
              obs = sefo_normalization (obs_hist, struct ("data_to_transform", obs));
              clear obs_hist
        endswitch

        %{
        #based on actual climatology ranks -- not very robust to almost-repeated values during climatology period
        obs = sort (obs, 2);
        R = ranks([obss(:, t) obs], 2);
        n_cats = size (R, 2);
        obs_cat = R(:, 1); #rank of current year's observations in the climatology [1 to n_clim_years+1]
        thresh_obs_cat = nan (np, 2);
        inds = floor(obs_cat) > 1;
        thresh_obs_cat(inds, 1) = obs(sub2ind([np n_cats], (1:np)'(inds), (floor(obs_cat)-1)(inds)));
        thresh_obs_cat(!inds, 1) = -Inf;
        inds = ceil(obs_cat) < n_cats;
        thresh_obs_cat(inds, 2) = obs(sub2ind([np n_cats], (1:np)'(inds), (ceil(obs_cat))(inds)));
        thresh_obs_cat(!inds, 2) = Inf;
        %}
        
        #based on assumed normality over the climatology period
        %n_cats = 1 + size (obs, 2);
        n_cats = 100;
        if t == 1
          obs_cats = zeros (np, nt, "uint16");
          thresh_obs_cats = nan (np, nt, 2, "single");
          thresh_cats = nan (np, nt, n_cats-1, "single");
        endif
        thresh_z = tinv ((1:(n_cats-1))'/n_cats, n_cats-1);
        clim_mean = mean(obs, 2);
        clim_sd = std(obs, [], 2);
        obs_cat = ceil (n_cats * tcdf ((obss(:, t) - clim_mean) ./ clim_sd, n_cats-1));
        obs_cat(obs_cat == 0) = 1; #cases where the z score is low enough that tcdf underflows to 0
        thresh_obs_cat = nan (np, 2);
        inds = (obs_cat == 1);
        thresh_obs_cat(!inds, 1) = clim_mean(!inds) + thresh_z(obs_cat(!inds) - 1) .* clim_sd(!inds);
        thresh_obs_cat(inds, 1) = -Inf;
        inds = (obs_cat == n_cats);
        thresh_obs_cat(!inds, 2) = clim_mean(!inds) + thresh_z(obs_cat(!inds)) .* clim_sd(!inds);
        thresh_obs_cat(inds, 2) = Inf;
                
        obs_cats(:, t) = obs_cat;
        thresh_obs_cats(:, t, :) = thresh_obs_cat;
        thresh_cats(:, t, :) = clim_mean + clim_sd * thresh_z';
        disp ([datestr(clock) "\n" 'sefo_verify: computed observation ranks'])
      endif
      disp ([datestr(clock) "\n" 'sefo_verify: compiled observations, t = ' num2str(t)])
    endfor
    
    obs = obss;
    clear obss
    if rank_calc
      freq_obs = histc (obs_cats(:), (1:n_cats)') / numel(obs_cats); #frequency of each rank in the observations
      freq_for = ll_for = nan (n_cats, n_methods);
    else
      freq_obs = freq_for = ll_for = [];
    endif
    disp ([datestr(clock) "\n" 'sefo_verify: computed observation frequencies'])
    
    methods_in = methods;
    rmse = nll = crps = predent = bias = ks = nan(1, n_methods);
    rmse_map = bias_map = nll_map = predent_map = nan(np, n_methods, "single");
    rmse_ts = bias_ts = nll_ts = predent_ts = nan(nt, n_methods, "single");
    predict_mu = predict_sigma = predict_nu = area_wgtss = nan(np, nt, "single");
    nq = 100;
    hists = nan (nq, n_methods);

    for m = 1:n_methods
      method = methods_in{m};
      disp([datestr(clock) "\n" 'sefo_verify: started method ' method])

      for t = 1:nt
        year = years(t);
        month = months(t);
        predict_month = month + (0:(n_months_combine-1))';
        predict_year = year * ones(n_months_combine, 1);
        inds_next_year = find(predict_month > 12);
        predict_year(inds_next_year) = predict_year(inds_next_year) + 1;
        predict_month(inds_next_year) = predict_month(inds_next_year) - 12;
        file_predict = sefo_predict (predict_year, predict_month, lag, method, opts){1};
        load (file_predict, "mu", "sigma", "nu", "area_wgts")
        predict_mu(:, t) = mu(:);
        predict_sigma(:, t) = sigma(:);
        predict_nu(:, t) = nu(:);
        if m == 1
          area_wgts = area_wgts / mean(double(area_wgts));
          area_wgtss(:, t) = area_wgts;
        endif
        disp ([datestr(clock) "\n" 'sefo_verify: compiled predictions, t = ' num2str(t)])
      endfor
      
      disp ([datestr(clock) "\n" 'sefo_verify: collected predictions'])

      if m == 1
        #area weighting and land/ocean mask
        w = mean (double (area_wgtss), 2);
        if !isempty(subset_tag)
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
          #get a land mask from BEST and regrid it to match forecast (NMME)
          source_file_name = ['Land_and_Ocean_LatLong1.nc'];
          if ~exist([data_dir source_file_name], 'file') || re_download
            data_url = ['http://berkeleyearth.lbl.gov/auto/Global/Gridded/' source_file_name];
            urlwrite(data_url, [data_dir source_file_name])  #download and save netcdf file 
          endif
          source_file = [data_dir source_file_name];
          lon_obs = ncread (source_file, 'longitude'); #-179.5:179.5
          lat_obs = ncread (source_file, 'latitude'); #-89.5:89.5
          land_mask = ncread (source_file, 'land_mask');
          arr = land_mask';
          arr = [arr arr arr];
          obs_lat_lims = sind([-90 lat_obs(1:(end-1))'+diff(lat_obs')/2 90])';
          obs_lon_lims = [lon_obs(end)-360 lon_obs' lon_obs(1)+360];
          obs_lon_lims =  obs_lon_lims(1:(end-1))'+diff(obs_lon_lims')/2;
          obs_lon_lims = [obs_lon_lims(1:(end-1))'-360 obs_lon_lims' obs_lon_lims(2:end)'+360]';
          coords = {obs_lat_lims, obs_lon_lims};
          new_coords = {fcst_lat_lims, fcst_lon_lims};
          land_mask_regrid = regrid (coords, arr, new_coords);
          #subset land or ocean
          switch subset_tag
            case '_land'
              w .*= land_mask_regrid(:);
            case '_ocean'
              w .*= (1 - land_mask_regrid(:));
            case ''
            otherwise
              error (['Undefined subset: ' subset_tag(2:end)])
          endswitch
        endif
        options.w = double (w);
        options.nq = 100;
        if rank_calc
          options.n_cats = n_cats;
          options.obs_cat = obs_cats;
          options.thresh_obs_cat = thresh_obs_cats;
        endif
      endif

      if predict_adj
        adj_file = sefo_adj (predict_adj_years, predict_adj_months, lag, method, opts);
        load(adj_file)
        factors_adj(m, :) = exp (facs_optim);
        predict_sigma = predict_sigma * factors_adj(m, 1);
        predict_nu = predict_nu / factors_adj(m, 2);
      endif
      
      if rank_calc #mean forecast probability of each climatology category
        thresh_probs = diff (cat (3, zeros(np, nt, 1), tcdf ((thresh_cats - predict_mu) ./ predict_sigma, repmat(predict_nu, [1 1 n_cats-1])), ones(np, nt, 1)), [], 3);
        freq_for(:, m) = mean (mean (thresh_probs, 1, double(w)), 2);
        ll_for(:, m) = mean (mean (log(thresh_probs), 1, double(w)), 2);
        clear thresh_probs
      endif
      
      disp ([datestr(clock) "\n" 'sefo_verify: calculating performance measures'])
      
      [struct_out] = sefo_verification (double(predict_mu), double(predict_sigma), double(predict_nu), double(obs), options);
      save_names = {'rmse', 'bias', 'nll', 'crps', 'predent', 'ks', 'kp', 'rmse_map', 'bias_map', 'nll_map', 'predent_map', ...
    'rmse_ts', 'bias_ts', 'nll_ts', 'predent_ts', 'hists', 'nll_cats', 'prob_cats'};
      for i = 1:numel(save_names)
        eval([save_names{i} '(:, m) = struct_out.' save_names{i} ';']);
      endfor

    endfor
    
    save (file, "rmse", "bias", "nll", "crps", "predent", "ks", "hists", "factors_adj", "rmse_map", "bias_map", "nll_map", "predent_map", "rmse_ts", "bias_ts", "nll_ts", "predent_ts", "nll_cats", "prob_cats", "freq_obs", "freq_for", "ll_for")
  endif
  
endfunction


%!demo
%! if !exist("data_dir", "var")
%!   data_dir = [pwd '/'];
%! endif
%! if !exist("methods", "var")
%!   methods = {'rawmodel', 'rawmodel_linear', 'EW-a'};
%! endif
%! opts = struct ('data_dir', data_dir, ...
%!                 're_download', false, ...
%!                 're_process', false, ...
%!                 'var_name', 'prec', ...
%!                 'fcst_lons', (0:359)', ...
%!                 'fcst_lats', (-90:90)', ... 
%!                 'fcst_start_year', 1982, ...
%!                 'obs_start_year', 1957, ...  
%!                 'obs_end_year', 2015, ...
%!                 'obs_dataset', 'gpcp', ...
%!                 'p_spline', 1E-4, ...
%!                 'dof_spline', 2.5, ...
%!                 'lamda', 10, ...
%!                 'predict_adj', false, ...
%!                 'predict_adj_years', 1994:1998, ...
%!                 'predict_adj_months', 1:12, ...                
%!                 'methods', {methods}, ...  
%!                 'clim_years', 1981:2010, ...   
%!                 'gcms', {{'CMC1-CanCM3', 'CMC2-CanCM4', 'GFDL-CM2p1-aer04', 'COLA-RSMAS-CCSM3', 'NASA-GMAO-062012'}} ...
%!                 );
%! years = 2012; months = 1; lag = 1;
%! file = sefo_verify (years, months, lag, opts)
%! load(file)
%! methods, rmse, bias, predent, nll, ks, freq_obs
%{
methods = 
{
  [1,1] = rawmodel
  [1,2] = rawmodel_linear
  [1,3] = trend
  [1,4] = model_trend
  [1,5] = multimethod
  
opts = struct ('data_dir', data_dir, ...
                 're_download', false, ...
                 're_process', false, ...
                 'var_name', 'tref', ...
                 'fcst_lons', (0:359)', ...
                 'fcst_lats', (-90:90)', ...
                 'fcst_start_year', 1982, ...
                 'obs_start_year', 1979, ...
                 'obs_end_year', 2015, ...
                 'obs_dataset', 'gpcp', ...
                 'p_spline', 1E-4, ...
                 'dof_spline', 2.5, ...
                 'lamda', 10, ...
                 'predict_adj', false, ...
                 'predict_adj_years', 1994:1998, ...
                 'predict_adj_months', 1:12, ...
                 'methods', {{'EW-a', 'MMT'}}, ...
                 'clim_years', 1981:2010, ...
                 'gcms', {{'CMC1-CanCM3', 'CMC2-CanCM4', 'COLA-RSMAS-CCSM4', 'GFDL-CM2p1-aer04', 'GFDL-CM2p5-FLOR-A06', 'GFDL-CM2p5-FLOR-B01', 'NASA-GMAO-062012', 'NCEP-CFSv2'}} ...
                 );
years = 2016; months = 1; lag = 1;
tic; file = sefo_verify (years, months, lag, opts); t = toc;
load(file)
rmse, bias, predent, nll, ks, freq_obs

}
rmse =

   0.83574   0.97047   1.14472   0.99649   0.98982

bias =

  -9.4099e-02  -1.7746e-02   8.5556e-04  -8.0852e-02  -7.9320e-02

predent =

   0.76428   0.84310   0.92438   0.77071   0.69417

nll =

   0.86605   0.96329   1.29008   0.98046   0.95606

ks =

   0.086916   0.084878   0.054888   0.081961   0.086308
  
  
  
%}

