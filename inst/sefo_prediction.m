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

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{mu}, @var{sigma}, @var{nu}] =} sefo_prediction (@var{fcst}, @var{obs}, @var{method}, @var{options})
## Compute a predictive distribution using past data and forecasts.
##
## It's assumed that forecasts from @var{k} different models are available for @var{nf} instances (time periods) at @var{m} sites,
## and that past observations are also available for @var{no} instances at the same @var{m} sites.
##
## The various methods implemented all are derived on the assumption of normally distributed data/forecast errors, 
## leading to predictive t distributions. Where distributions are not normal, pre-transforming the data and forecasts may be useful.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{fcst} (dimensions: @var{m} by @var{nf} by @var{k}) should contain the forecast values
## @item
## @var{obs} (dimensions: @var{m} by @var{no}) should contain the observations
## @item
## @var{method} (string), which prediction method to use
## @item
## @var{options} (structure), optional additional arguments. Default values are used for any that are unspecified or empty.
## @end itemize
##
## @subheading @var{options} fields
##
## @itemize @bullet
## @item
## @var{t_fcst} (dimensions: @var{nf} by 1), times (instances) at which forecasts are available. Default: @code{(1:nf)'}
## @item
## @var{t_obs} (dimensions: @var{no} by 1), instances at which observations are available. Default: @code{((nf-no):(nf-1))'}
## @item
## @var{t_clim} (dimensions: @var{nc} by 1), observation instances that should be used to compute climatology. Default: @code{(1:no)'}
## @item
## @var{t_pred} (scalar), instance at which the predictions should be made. Default: @code{nf}
## @item
## @var{w} (dimensions: @var{m} by 1), site weights to use for computing global parameters. Default: @code{ones(m, 1)} (equal weight)
## @item
## @var{timescale} (scalar), consider for prediction only observations and forecasts that are recent compared to this
## (used in moving average -type methods). Default: 30
## @item
## @var{dof_spline} (scalar), degrees of freedom to assume that are used by trend fits (in methods that employ them). Default: 2.5
## @item
## @var{lamda} (scalar), parameter controlling the degree of smoothing of parameters across sites (in methods that employ it). Default: 10
## @item
## @var{fcst_trend} (dimensions: @var{m} by @var{nf} by @var{k}) trend component of forecasts (used in trend-based methods).
## Default: @code{repmat(mean(fcst, 2), [1 nf 1])}
## @item
## @var{obs_trend} (dimensions: @var{m} by @var{no+1}) trend component of observations 
## (with the last column being the estimated trend at the prediction time). Default: @code{repmat(mean(obs, 2), [1 no+1 1])}
## @item
## @var{mus, sigmas, nus} (dimensions of each: @var{m} by @var{n_methods}) predictions from other methods to combine for @code{multimethods}. Default: empty
## @end itemize
##
## @subheading Implemented methods
##
## @itemize @bullet
## @item
## @var{history} or @var{C}: Climatology (assumes stationarity) -- average and standard deviation of all observations (or those with indices given in @var{t_clim}).
## @item
## @var{linear} or @var{T}: Regression model using all observations and linear trend.
## @item
## @var{quadratic}: Regression model using all observations and quadratic trend.
## @item
## @var{ma} or @var{MA}: Climatology based only on observations no older than @var{timescale}.
## @item
## @var{ma_globcorr}: @var{ma} with a correction for the effect of recent trends based on a spline fit to the multi-site mean.
## @item
## @var{spline}: Smoothing spline fit.
## @item
## @var{spline_scaled}: Smoothing spline fit of time series normalized by their standard deviation -- useful if variance differs across sites.
## @item
## @var{ewma} or @var{EW}: Climatology based on past observations being exponentially downweighted with scale @var{timescale}.
## @item
## @var{ewma_globcorr} or @var{EW-a}: @var{ewma} with a correction for the effect of recent trends based on a spline fit to the multi-site mean.
## @item
## @var{ewma_linear}: Regression with linear trend, with past observations being exponentially downweighted with scale @var{timescale}.
## @item
## @var{ewma_linear_globcorr}: @var{ewma_linear} with a correction for the effect of recent trends based on a spline fit to the multi-site mean.
## @item
## @var{ewma_model}: Regression with mean forecast, with past observations being exponentially downweighted with scale @var{timescale}.
## @item
## @var{ewma_linear_model_reg}: Regression with linear trend and mean forecast, with past observations being exponentially downweighted with scale @var{timescale} and parameters smoothed toward their across-site means.
## @item
## @var{ewma_linear_multimodel_reg}: Regression with linear trend and forecasts, with past observations being exponentially downweighted with scale @var{timescale} and parameters smoothed toward their across-site means.
## @item
## @var{ewma_linear_model}: Regression with linear trend and mean forecast, with past observations being exponentially downweighted with scale @var{timescale}.
## @item
## @var{model_linear_reg}: Regression model using linear trend and mean forecast, with parameters smoothed toward their across-site means.
## @item
## @var{multimodel_linear_reg}: Regression model using linear trend and forecasts, with parameters smoothed toward their across-site means.
## @item
## @var{rawmodel} or @var{M}: Mean forecast, with offset determined by regression model.
## @item
## @var{rawmodel_globlinear} or @var{MT}: Mean forecast, with offset determined by regression model and a globally determined linear trend.
## @item
## @var{rawmodel_scaled}: Regression model using mean forecast.
## @item
## @var{rawmodel_globscaled} or @var{MS}: Mean Forecast with globally determined weight, offset determined by regression model.
## @item
## @var{rawmodel_globscaled_globlinear} or @var{MST}: Mean Forecast with globally determined weight and linear trend, offset determined by regression model.
## @item
## @var{rawmultimodel_glob} or @var{MM}: Forecasts with globally determined weights, offset determined by regression model.
## @item
## @var{rawmultimodel_glob_globlinear} or @var{MMT}: Forecasts and linear trend with globally determined weights, offset determined by regression model.
## @item
## @var{rawmultimodel}: Regression model using forecasts.
## @item
## @var{rawmodel_pooled}: Regression model using mean forecast, with partial pooling of coefficients.
## @item
## @var{rawmodel_linear_pooled}: Regression model using mean forecast and linear trend, with partial pooling of coefficients.
## @item
## @var{rawmultimodel_pooled}: Regression model using forecasts, with partial pooling of coefficients.
## @item
## @var{rawmodel_linear}: Mean forecast, with offset and linear trend determined by regression model.
## @item
## @var{rawmodel_globcorr}: Mean forecast, with offset and globally determined linear trend.
## @item
## @var{trend}: Forecasts the trend, with variability from previous residuals.
## @item
## @var{model_trend}: Regression model using detrended mean forecast to predict detrended observations.
## @item
## @var{multimodel_trend}: Regression model using detrended forecasts to predict detrended observations.
## @item
## @var{globalmultimodel_trend}: Regression model using detrended forecasts, with globally determined weights, to predict detrended observations.
## @item
## @var{regmultimodel_trend}: Regression model using detrended forecasts to predict detrended observations, with parameters smoothed toward their across-site means.
## @item
## @var{multimethod}: For each grid point, select the prediction method the claims most confidence (lowest entropy of predictive distribution).
## @end itemize
##
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{mu} (dimensions: @var{m} by 1), means of the predictive t distributions
## @item
## @var{sigma} (dimensions: @var{m} by 1), scales of the predictive t distributions
## @item
## @var{nu} (dimensions: @var{m} by 1, or scalar), degrees of freedom of the predictive t distributions
## @end itemize
##
## @subheading Example
##
## @example
## @group
## m = 10; nf = 20; no = 30; k = 3;
## obs = randn(m, no+1); coeffs = 1 + randn(1, 1, k);
## fcst = repmat(obs(:, (no-nf+2):end), [1 1 k]) .* coeffs + randn(m, nf, k);
## #obs = randn(m, no+1) + 0.1*(ones(m, 1)*(1:no+1)); #to add a linear trend not accounted for in the forecasts
## method = 'rawmodel_globlinear'; [mu, sigma, nu] = sefo_prediction (fcst, obs(:, 1:(end-1)), method);
## [rms(mu - obs(:, end)) rms(mean(obs(:, 1:(end-1)), 2) - obs(:, end))]
## @end group
## @end example
##
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>

#methods_list = {'history', 'linear', 'quadratic', 'ma', 'ma_globcorr', 'spline', 'spline_scaled', 'ewma', 'ewma_globcorr', 'ewma_linear', 'ewma_linear_globcorr', 'ewma_model', 'ewma_linear_model_reg', 'ewma_linear_model', 'model_linear_reg', 'multimodel_linear_reg', 'rawmodel', 'rawmodel_scaled', 'rawmodel_globscaled', 'rawmultimodel_glob', 'rawmultimodel_glob_linear', 'rawmultimodel', 'rawmodel_pooled', 'rawmodel_linear_pooled', 'rawmultimodel_pooled', 'rawmodel_linear', 'rawmodel_globcorr', 'trend', 'model_trend', 'multimodel_trend', 'globalmultimodel_trend', 'regmultimodel_trend', 'multimethod'};

function [mu, sigma, nu] = sefo_prediction (fcst, obs, method, options)

  warning ("off", "Octave:broadcast", "global");

  [m nf k] = size(fcst);
  no = size(obs, 2);

  #read in options or use default values

  options_default = struct('t_fcst', (1:nf)', ...
                          't_obs', ((nf-no):(nf-1))', ...
                          't_clim', (1:no)', ...
                          't_pred', nf, ...
                          'w', ones(m, 1), ...
                          'timescale', 30, ...
                          'dof_spline', 2.5, ...
                          'lamda', 10, ...
                          'fcst_trend', repmat(mean(fcst, 2), [1 nf 1]), ...
                          'obs_trend', repmat(mean(obs, 2), [1 no+1 1]), ...
                          'mus', [], ...
                          'sigmas', [], ...
                          'nus', []);

  option_names = fieldnames(options_default);
  if nargin > 3
    for i = 1:numel(option_names)
      option_name = option_names{i};
      if isfield(options, option_name) && !isempty(options.(option_name))
        options_default.(option_name) = options.(option_name);
      endif
    endfor
    clear options
  endif
  for i = 1:numel(option_names)
    option_name = option_names{i};
    eval([option_name " = getfield (options_default, '" option_name "');"])
  endfor
  clear options_default option_names

    switch method
      case {'history', 'C'}
        [mu, sigma, nu] = ma_params (obs(:, t_clim)', numel(t_clim));
        mu = mu';
        sigma = sigma';
      case {'linear', 'T'}
        t_mean = mean(t_obs);
        t_obs = t_obs - t_mean;
        t_pred = t_pred - t_mean;
        [mu, sigma, nu] = regress_params ([ones(no, 1) t_obs], obs', [1 t_pred]);
        mu = mu';
        sigma = sigma';
      case 'quadratic'
        t_mean = mean(t_obs);
        t_obs = t_obs - t_mean;
        t_pred = t_pred - t_mean;
        [mu, sigma, nu] = regress_params ([ones(no, 1) t_obs t_obs.^2], obs', [1 t_pred t_pred.^2]);
        mu = mu';
        sigma = sigma';
     case {'ma', 'MA'}
        t_inds = abs(t_pred - t_obs) <= timescale;
        obs = obs(:, t_inds);
        [mu, sigma, nu] = ma_params (obs', no);
        mu = mu';
        sigma = sigma';
      case 'ma_globcorr'
        obs_global = (w' * obs) / sum(w);
        
        t_inds = abs(t_pred - t_obs) <= timescale;
        obs = obs(:, t_inds);
        [mu, sigma, nu] = ma_params (obs', no);
              
        #global correction to mu
        [mu_global, sigma_global, nu_global] = ma_params (obs_global(t_inds)', no);
        [mu_global_spline] = csaps_sel (t_obs, obs_global', t_pred);
        
        mu = mu' + (mu_global_spline - mu_global);
        sigma = sigma';
      case 'spline'
        [mu, sigma, nu, p] = spline_params (t_obs, obs', t_pred);
      case 'spline_scaled'
        s = std (obs, [], 2);
        obs ./= s;
        [mu, sigma, nu, p] = spline_params (t_obs, obs', t_pred);
        mu .*= s;
        sigma .*= s;
      case {'ewma', 'EW'}
        d = 0;
        gamma = 1 ./ (1 + timescale .^ 2);
        
        v = t_pred - t_obs;
        R = min(v, v');
        [P, D] = eig (R); #presumably faster to pre-calculate this
        mu = sigma = nu = nan (m, 1);
        for p = 1:m
          [mu(p), sigma(p), nu(p)] = ew_params (t_obs, obs(p, :)', t_pred, gamma, d, {P, D});
        endfor

      case {'ewma_globcorr', 'EW-a'}
        d = 0;
        gamma = 1 ./ (1 + timescale .^ 2);
        
        v = t_pred - t_obs;
        R = min(v, v');
        [P, D] = eig (R); #presumably faster to pre-calculate this
        mu = sigma = nu = nan (m, 1);
        for p = 1:m
          [mu(p), sigma(p), nu(p)] = ew_params (t_obs, obs(p, :)', t_pred, gamma, d, {P, D});
        endfor
              
        #global correction to mu
        obs_global = (w' * obs) / sum(w);
        [mu_global, sigma_global, nu_global] = ew_params (t_obs, obs_global', t_pred, gamma, d, {P, D});
        [mu_global_spline] = csaps_sel (t_obs, obs_global', t_pred);
        
        mu = mu + (mu_global_spline - mu_global);
        
      case 'ewma_linear'
        d = 1;
        gamma = 1 ./ (1 + timescale .^ 2);
        
        t_mean = mean(t_obs);
        v = t_pred - t_obs;
        t_obs = t_obs - t_mean;
        t_pred = t_pred - t_mean;
        R = min(v, v');
        [P, D] = eig (R); #presumably faster to pre-calculate this
        mu = sigma = nu = nan (m, 1);
        for p = 1:m
          [mu(p), sigma(p), nu(p)] = ew_params (t_obs, obs(p, :)', t_pred, gamma, d, {P, D});
        endfor

      case 'ewma_linear_globcorr'
        d = 1;
        gamma = 1 ./ (1 + timescale .^ 2);
        
        t_mean = mean(t_obs);
        v = t_pred - t_obs;
        t_obs = t_obs - t_mean;
        t_pred = t_pred - t_mean;
        R = min(v, v');
        [P, D] = eig (R); #presumably faster to pre-calculate this
        mu = sigma = nu = nan (m, 1);
        for p = 1:m
          [mu(p), sigma(p), nu(p)] = ew_params (t_obs, obs(p, :)', t_pred, gamma, d, {P, D});
        endfor
                      
        #global correction to mu
        obs_global = (w' * obs) / sum(w);
        [mu_global, sigma_global, nu_global] = ew_params (t_obs, obs_global', t_pred, gamma, d, {P, D});
        [mu_global_spline] = csaps_sel (t_obs, obs_global', t_pred);
        
        mu = mu + (mu_global_spline - mu_global);

      case 'ewma_model'

        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        
        gamma = 1 ./ (1 + timescale .^ 2);
        
        v = t_pred - t;
        R = min(v, v');
        [P, D] = eig (R); #presumably faster to pre-calculate this
        mu = sigma = nu = nan (m, 1);
        fcst_mean = mean(fcst, 3) - mean(fcst_pred, 3);
        for p = 1:m
          X = [ones(n, 1) fcst_mean(p, :)'];
          [mu(p), sigma(p), nu(p)] = ew_params ([], obs(p, :)', [], gamma, X, {P, D});
        endfor

      case 'ewma_linear_model_reg'
        
        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        
        gamma = 1 ./ (1 + timescale .^ 2);
        
        v = t_pred - t;
        R = min(v, v');

        C = (1-gamma)*speye(n) + gamma*R;
        L = [zeros(1, 2) speye(1)];

        fcst_mean = mean(fcst, 3) - mean(fcst_pred, 3);

        W = (w * ones(1, n))(:);
        c = ([ones(m*n, 1) (ones(m, 1)*v')(:) center(fcst_mean, 2)(:)].*W) \ (obs(:).*W); #global-mean value for best model regression coefficient
        beta0 = [0; 0; c(end)];
        Xstar = [1 0 0];
        
        mu = sigma = nu = nan (m, 1);

        for p = 1:m
          X = [ones(n, 1) v fcst_mean(p, :)'];
          y = obs(p, :)';
          [mu(p), sigma(p), nu(p)] = regress_reg_gen_aicc_params (X, y, Xstar, C, L, beta0);
        endfor
      
      case 'ewma_linear_multimodel_reg'

        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);

        gamma = 1 ./ (1 + timescale .^ 2);
        
        v = t_pred - t;
        R = min(v, v');

        C = (1-gamma)*speye(n) + gamma*R;
        L = [zeros(k, 2) speye(k)];


        fcst = fcst - fcst_pred;

        W = (w * ones(1, n))(:);
        c = ([ones(m*n, 1) (ones(m, 1)*v')(:) reshape(center(fcst, 2), m*n, k)].*W) \ (obs(:).*W); #global-mean value for best model regression coefficients
        beta0 = [0; 0; c((end-k+1):end)];
        Xstar = [1 0 zeros(1, k)];
        
        mu = sigma = nu = nan (m, 1);

        for p = 1:m
          X = [ones(n, 1) v squeeze(fcst(p, :, :))];
          y = obs(p, :)';
                   
          [mu(p), sigma(p), nu(p)] = regress_reg_gen_aicc_params (X, y, Xstar, C, L, beta0);
        endfor

      case 'ewma_linear_model'

        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
      
        gamma = 1 ./ (1 + timescale .^ 2);
        
        v = t_pred - t;
        R = min(v, v');
        [P, D] = eig (R);
        mu = sigma = nu = nan (m, 1);

        fcst_mean = mean(fcst, 3) - mean(fcst_pred, 3);

        for p = 1:m
          X = [ones(n, 1) v fcst_mean(p, :)(:)];
          [mu(p), sigma(p), nu(p)] = ew_params ([], obs(p, :)', [], gamma, X, {P, D});
        endfor

      case 'model_linear_reg' #multimodel mean (regularized, with offset and linear trend)
      
        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        v = t_pred - t;
      
        fcst_mean = mean(fcst, 3) - mean(fcst_pred, 3);
        
        W = (w * ones(1, n))(:);
        
        c = ([ones(m*n, 1) (ones(m, 1)*v')(:) center(fcst_mean, 2)(:)].*W) \ (obs(:).*W); #global-mean value for best model regression coefficient
        
        mu = sigma = nu = nan (m, 1);
      
        for p = 1:m
          X = [ones(n, 1) v fcst_mean(p, :)'];
          Xstar = [1 0 0];
          y = obs(p, :)';
          C = speye (n);
          beta0 = [0 0 c(end)]';
          L = [zeros(1, 2) eye(1, 1)];
          [mu(p), sigma(p), nu(p)] = regress_reg_gen_aicc_params (X, y, Xstar, C, L, beta0);

        endfor

      case 'multimodel_linear_reg' #multimodel mean (regularized, with offset and linear trend)

        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        v = t_pred - t;
 
        fcst = fcst - fcst_pred;
      
        W = (w * ones(1, n))(:);
      
        c = ([ones(m*n, 1) (ones(m, 1)*v')(:) reshape(center(fcst, 2), m*n, k)].*W) \ (obs(:).*W); #global-mean value for best model regression coefficients
        
        mu = sigma = nu = nan (m, 1);
     
        for p = 1:m
          X = [ones(n, 1) v squeeze(fcst(p, :, :))];
          Xstar = [1 0 zeros(1, k)];
          y = obs(p, :)';
          C = speye (n);
          beta0 = [0; 0; c(3:end)];
          L = [zeros(k, 2) eye(k)];
          [mu(p), sigma(p), nu(p)] = regress_reg_gen_aicc_params (X, y, Xstar, C, L, beta0);
        endfor

     case {'rawmodel', 'M'} #multimodel mean (with offset) T ~ alpha_x + M
        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        fcst_mean = mean(fcst, 3) - mean(fcst_pred, 3);
        [mu, sigma, nu] = regress_params ([ones(n, 1)], (obs - fcst_mean)', 1);
        mu = mu';
        sigma = sigma';

     case {'rawmodel_globlinear', 'MT'} #multimodel mean (with offset) and global linear trend T ~ alpha_x + M + gamma*t
        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        fcst_mean = mean(fcst, 3) - mean(fcst_pred, 3);
        Y = center (obs - fcst_mean, 2);
        X = (ones(m, 1) * center(t)');
        W = (w * ones(1, n));
        b = (X .* W)(:) \ (Y .* W)(:);
        [mu, sigma, nu] = regress_params ([ones(n, 1)], (obs - fcst_mean)' - b*t, 1);
        mu = mu' + b*t_pred;
        sigma = sigma';

      case 'rawmodel_scaled' #T ~ alpha_x + beta_x*M

        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        fcst_mean = mean(fcst, 3) - mean(fcst_pred, 3);
        mu = sigma = nu = nan (m, 1);
        for p = 1:m
          [mu(p), sigma(p), nu(p)] = regress_params ([ones(n, 1) fcst_mean(p, :)'], obs(p, :)', [1 0]);
        endfor
     
      case {'rawmodel_globscaled', 'MS'} #T ~ alpha_x + beta*M
        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        fcst_mean = mean(fcst, 3) - mean(fcst_pred, 3);
        W = (w * ones(1, n));
        b = (center(fcst_mean, 2) .* W)(:) \ (center(obs, 2) .* W)(:);
        [mu, sigma, nu] = regress_params ([ones(n, 1)], (obs - b*fcst_mean)', 1);
        mu = mu';
        sigma = sigma';

      case {'rawmodel_globscaled_globlinear', 'MST'} #T ~ alpha_x + beta*M + gamma*t
        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        fcst_mean = mean(fcst, 3) - mean(fcst_pred, 3);
        X = center(fcst_mean, 2);
        X(:, :, 2) = (ones(m, 1) * center(t)');
        Y = center (obs, 2);
        W = (w * ones(1, n));
        b = reshape(X .* W, m*n, 2) \ (Y .* W)(:);
        [mu, sigma, nu] = regress_params ([ones(n, 1)], (obs - b(1)*fcst_mean)' - b(2)*t, 1);
        mu = mu' + b(2)*t_pred;
        sigma = sigma';

      case {'rawmultimodel_glob', 'MM'} #T ~ alpha_x + sum(beta_i*M_i)
        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        fcst = fcst - fcst_pred;
        W = (w * ones(1, n));
        
        b = reshape(center(fcst, 2) .* W, m*n, k) \ (center(obs, 2) .* W)(:);
        [mu, sigma, nu] = regress_params ([ones(n, 1)], (obs - sum(fcst .* reshape(b, [1 1 k]), 3))', 1);
        mu = mu';
        sigma = sigma';

      case {'rawmultimodel_glob_globlinear', 'MMT'} #T ~ alpha_x + sum(beta_i*M_i) + gamma*t
        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        fcst = fcst - fcst_pred;
        W = (w * ones(1, n));
        
        b = [reshape(center(fcst, 2) .* W, m*n, k) (w*center(t)')(:)] \ (center(obs, 2) .* W)(:);
        [mu, sigma, nu] = regress_params ([ones(n, 1)], (obs - sum(fcst .* reshape(b(1:end-1), [1 1 k]), 3) - b(end)*ones(m, 1)*t')' , 1);
        mu = mu' + b(end)*t_pred;
        sigma = sigma';

      case 'rawmultimodel' #T ~ alpha_x + sum(beta_xi*M_i)

        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        fcst = fcst - fcst_pred;

        mu = sigma = nu = nan (m, 1);
        for p = 1:m
          [mu(p), sigma(p), nu(p)] = regress_params ([ones(n, 1) squeeze(fcst(p, :, :))], obs(p, :)', [1 zeros(1, k)]);
        endfor

      case 'rawmodel_pooled' #T ~ beta_x1 + beta_x2*M, beta_x ~ beta_0

        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        F = center(mean([fcst fcst_pred], 3), 2)';
        [mu, sigma, nu, beta0, S0] = regress_params_pooled ([ones(n, 1, m) reshape(F(1:(end-1), :), [n 1 m])], obs', [ones(1, 1, m) reshape(F(end, :), [1 1 m])]);

      case 'rawmodel_linear_pooled' #T ~ beta_x1 + beta_x2*t + beta_x3*M, beta_x ~ beta_0

        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        F = center(mean([fcst fcst_pred], 3), 2)';
        [mu, sigma, nu, beta0, S0] = regress_params_pooled ([ones(n, 1, m) repmat(t, [1 1 m]) reshape(F(1:(end-1), :), [n 1 m])], obs', [ones(1, 1, m) t_pred*ones(1, 1, m) reshape(F(end, :), [1 1 m])]);
      
      case 'rawmultimodel_pooled' #T ~ beta_x1 + sum(beta_xi*M_i), beta_x ~ beta_0

        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        
        F = permute(center([fcst fcst_pred], 2), [2 3 1]);
        
        [mu, sigma, nu, beta0, S0] = regress_params_pooled ([ones(n, 1, m) F(1:(end-1), :, :)], obs', [ones(1, 1, m) F(end, :, :)]);
            
      case 'rawmodel_linear' #multimodel mean (with offset and linear trend) T ~ alpha_x + M + gamma_x*t
        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        fcst_mean = mean(fcst, 3) - mean(fcst_pred, 3);
        [mu, sigma, nu] = regress_params ([ones(n, 1) t], (obs - fcst_mean)', [1 t_pred]);
        mu = mu';
        sigma = sigma';
      
      case 'rawmodel_globcorr' #multimodel mean (with offset and globally-estimated linear trend)

        fcst_pred = fcst(:, t_fcst == t_pred, :);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        obs = obs(:, i_obs, :);
        n = numel(t);
        fcst_mean = mean(fcst, 3) - mean(fcst_pred, 3);
        
        obs_global = (w' * obs) / sum(w);
        fcst_mean_global = sum(reshape(w, [m 1]) .* fcst_mean, 1) / sum(w);

        b = [ones(n, 1) center(t)] \ (obs_global - fcst_mean_global)';
            
        [mu, sigma, nu] = regress_params ([ones(n, 1)], (obs - fcst_mean)' - b(2)*t, 1);
        mu = mu' + b(2)*t_pred;
        sigma = sigma';

      case 'trend' #trend of past observations
        [mu, sigma, nu] = ma_params ((obs - obs_trend(:, 1:(end-1)))', no);
        mu = mu' + obs_trend(:, end); #actually, at least if the observations are equally spaced, mu will always be 0 up to roundoff error 
        sigma = sigma' * sqrt(nu/(no-dof_spline)); #adjust standard deviation and effective degrees of freedom in an attempt to account for spline fitting of the historical observations
        nu = no - dof_spline;
        
      case 'model_trend' #trend plus detrended model mean as a predictor

        fcst_pred = fcst(:, t_fcst == t_pred, :);
        fcst_trend_pred = fcst_trend(:, t_fcst == t_pred, :);
        obs_trend_pred = obs_trend(:, end);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        fcst_trend = fcst_trend(:, i_fcst, :);
        obs = obs(:, i_obs);
        obs_trend = obs_trend(:, i_obs);
        n = numel(t);
        fcst_mean_detrend = (mean(fcst, 3) - mean(fcst_pred, 3)) - (mean(fcst_trend, 3) - mean(fcst_trend_pred, 3));
        obs_detrend = obs - (obs_trend - obs_trend_pred);
      
        mu = sigma = nu = nan (m, 1);
        for p = 1:m
          [mu(p), sigma(p), nu(p)] = regress_params ([ones(n, 1) fcst_mean_detrend(p, :)'], obs_detrend(p, :)', [1 0]);
        endfor
        
      case 'multimodel_trend' #trend plus detrended models as predictors
        fcst_pred = fcst(:, t_fcst == t_pred, :);
        fcst_trend_pred = fcst_trend(:, t_fcst == t_pred, :);
        obs_trend_pred = obs_trend(:, end);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        fcst_trend = fcst_trend(:, i_fcst, :);
        obs = obs(:, i_obs);
        obs_trend = obs_trend(:, i_obs);
        n = numel(t);
        fcst_detrend = (fcst - fcst_pred) - (fcst_trend - fcst_trend_pred);
        obs_detrend = obs - (obs_trend - obs_trend_pred);
      
        mu = sigma = nu = nan (m, 1);
        for p = 1:m
            [mu(p), sigma(p), nu(p)] = regress_params ([ones(n, 1) squeeze(fcst_detrend(p, :, :))], obs_detrend(p, :)', [1 zeros(1, k)]);
        endfor
           
      case 'globalmultimodel_trend' #trend plus detrended models as predictors, with globally constant weights
        fcst_pred = fcst(:, t_fcst == t_pred, :);
        fcst_trend_pred = fcst_trend(:, t_fcst == t_pred, :);
        obs_trend_pred = obs_trend(:, end);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        fcst_trend = fcst_trend(:, i_fcst, :);
        obs = obs(:, i_obs);
        obs_trend = obs_trend(:, i_obs);
        n = numel(t);
        fcst_detrend = (fcst - fcst_pred) - (fcst_trend - fcst_trend_pred);
        obs_detrend = obs - (obs_trend - obs_trend_pred);
      
        #globally optimal model weights
        X = reshape(fcst_detrend, n*m, k);
        X = [X ones(n*m, 1)];
        y = obs_detrend(:);
        W = (w * ones(1, n))(:);
        c = (X.*W) \ (y.*W);
      
        [mu, sigma, nu] = ma_params ((obs_detrend - reshape(X*c, m, n))', n);
        mu = mu'  + [squeeze(fcst_pred) ones(m, 1)]*c;
        #mu = mu' + c(end);
        sigma = sigma' * sqrt(nu/(n-dof_spline));
        nu = n - dof_spline;
        
      case 'regmultimodel_trend' #trend plus detrended models as predictors, regularized
      
        fcst_pred = fcst(:, t_fcst == t_pred, :);
        fcst_trend_pred = fcst_trend(:, t_fcst == t_pred, :);
        obs_trend_pred = obs_trend(:, end);
        [t, i_fcst, i_obs] = intersect (t_fcst, t_obs);
        fcst = fcst(:, i_fcst, :);
        fcst_trend = fcst_trend(:, i_fcst, :);
        obs = obs(:, i_obs);
        obs_trend = obs_trend(:, i_obs);
        n = numel(t);
        fcst_detrend = (fcst - fcst_pred) - (fcst_trend - fcst_trend_pred);
        obs_detrend = obs - (obs_trend - obs_trend_pred);
      
        #globally optimal model weights
        X = reshape(fcst_detrend, n*m, k);
        X = [X ones(n*m, 1)];
        y = obs_detrend(:);
        W = (w * ones(1, n))(:);
        c = (X.*W) \ (y.*W);
      
        mu = sigma = nu = nan (m, 1);
        for p = 1:m
            [mu(p), sigma(p), nu(p)] = regress_params ([[ones(n, 1) squeeze(fcst_detrend(p, :, :))]' [zeros(k, 1) lamda*eye(k)]']', [obs_detrend(p, :) lamda*c(1:k)']', [1 zeros(1, k)]);
        endfor
        
      case 'multimethod' #for each grid point, select the method the claims most confidence (lowest entropy of PDF)
        if isempty(mus) || isempty(sigmas) || isempty(nus)
          warning('sefo_prediction: no results from other methods given to multimethod, returning empty fields')
          mu = sigma = nu = [];
          return
        endif
        n_methods = size (mus, 2);
        
        for m = 1:n_methods
          mu = mus(:, m);
          sigma = sigmas(:, m);
          nu = nus(:, m);
          if m == 1
            mu_met = mu;
            sigma_met = sigma;
            nu_met = nu;
            ent_met = t_entropy (nu) + log(sigma);
          else
            ent = t_entropy (nu) + log(sigma);
            ii = ent < ent_met;
            if any(ii)
              mu_met(ii) = mu(ii);
              sigma_met(ii) = sigma(ii);
              nu_met(ii) = nu(ii);
              ent_met(ii) = ent(ii);
            endif
          endif
        endfor
        mu = mu_met;
        sigma = sigma_met;
        nu = nu_met;
               
      otherwise
        error(['method ' method ' not implemented; see function help text for the current list'])
    endswitch


endfunction


%{
m = 1000; nf = 20; no = 30; k = 3;
obs = randn(m, no+1); coeffs = 1 + randn(1, 1, k);
fcst = repmat(obs(:, (no-nf+2):end), [1 1 k]) .* coeffs + randn(m, nf, k);
method = 'rawmultimodel_glob'; [mu, sigma, nu] = sefo_prediction (fcst, obs(:, 1:(end-1)), method);
n_cats = 10; thresh_obs_cat = [-Inf; quantile(obs(:, 1:end-1)(:), (1:(n_cats-1))'/n_cats); Inf]; obs_cat = reshape(interp1(thresh_obs_cat, (1:n_cats+1)', obs(:, end), "previous"), [m 1]);
options.obs_cat = obs_cat; options.thresh_obs_cat = thresh_obs_cat; options.n_cats = n_cats;
[struct_out] = sefo_verification (mu, sigma, nu, obs(:, end), options);


m = 10000; nf = 30; no = 30; k = 1;
obs = randn(m, no+1); coeffs = 1 + randn(1, 1, k);
fcst = nan(m, nf, k);
method = 'history'; [mu, sigma, nu] = sefo_prediction (fcst, obs(:, 1:(end-1)), method);
n_cats = 10; thresh_obs_cat = [-Inf; quantile(obs(:, 1:end-1)(:), (1:(n_cats-1))'/n_cats); Inf]; obs_cat = reshape(interp1(thresh_obs_cat, (1:n_cats+1)', obs(:, end), "previous"), [m 1]);
options.obs_cat = obs_cat; options.thresh_obs_cat = thresh_obs_cat; options.n_cats = n_cats;
[struct_out] = sefo_verification (mu, sigma, nu, obs(:, end), options);
%}
