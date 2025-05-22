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
## @deftypefn {Function File} {[@var{struct_out}] =} sefo_verification (@var{pred_mu}, @var{pred_sigma}, @var{pred_nu}, @var{obs}, @var{options})
## Compare t-distribution predictions with data.
##
## It's assumed that forecasts using @var{k} different methods are available for @var{n} instances (time periods) at @var{m} sites,
## and that observations are also available for the same @var{n} instances at the same @var{m} sites.
##
## The various methods implemented all are derived on the assumption of normally distributed data/forecast errors, 
## leading to predictive t distributions. Where distributions are not normal, pre-transforming the data and forecasts may be useful.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{fcst_mu}, @var{fcst_sigma}, @var{fcst_nu} (dimensions: @var{m} by @var{n} by @var{k}), respectively
## the mean, scale parameter, and degrees of freedom of the forecast distributions
## @item
## @var{obs} (dimensions: @var{m} by @var{n}), the observations
## @item
## @var{options} (structure), optional additional arguments. Default values are used for any that are unspecified or empty.
## @end itemize
##
## @subheading @var{options} fields
##
## @itemize @bullet
## @item
## @var{w} (dimensions: @var{m} by 1), site weights to use for computing global means. Default: @code{ones(m, 1)} (equal weight)
## @item
## @var{nq} (positive integer), number of rank bins to compute histogram for. Default: 100
## @item
## @var{obs_cat} (positive integers, dimensions: @var{m} by @var{n}). Category into which each observed value falls, relative to some set of thresholds. Default: empty.
## @item
## @var{thresh_obs_cat} (dimensions: @var{m} by @var{n} by 2, or @var{n_cats} by 2 or @var{n_cats+1} by 1). Upper and lower bounds for inclusion in each observation's category. Default: empty.
## @item
## @var{n_cats} (positive integer). Number of categories of observations used. Default: max(obs_cat(:)).
## @end itemize
##
## @subheading Fields of output structure
##
## @itemize @bullet
## @item
## @var{rmse}, @var{bias}, @var{nll}, @var{predent}, @var{ks}, @var{kp} (dimensions: @var{k} by 1):
## Respectively prediction root mean square error, prediction mean bias, mean negative log likelihood of observations, mean entropy of predictive distribution,
## and Kolmogorov-Smirnov statistic and p value of locations of observations in the predicted CDF compared to a uniform distribution
## (all averaged across sites and times).
## @item
## @var{rmse_ts}, @var{bias_ts}, @var{nll_ts}, @var{predent_ts} (dimensions: @var{n} by @var{k}):
## Respectively prediction root mean square error, mean negative log likelihood of observations, mean entropy of predictive distribution
## (all averaged across sites).
## @item
## @var{rmse_map}, @var{bias_map}, @var{nll_map}, @var{predent_map} (dimensions: @var{m} by @var{k}):
## Respectively prediction root mean square error, mean negative log likelihood of observations, mean entropy of predictive distribution
## (all averaged across times).
## @item
## @var{hists} (dimensions: @var{nq} by @var{k}): Histogram of how many observations are in each centile (etc.) of the prediction CDF.
## @item
## @var{nll_cats} (dimensions: @var{n_cats} by @var{k}): Mean negative log likelihood conditional on the observations falling in a particular category.
## @item
## @var{prob_cats} (dimensions: @var{n_cats} by @var{k}): Mean correct-category forecast probability conditional on the observations falling in a particular category.
## [@var{nll_cats}, @var{prob_cats} are only computed if @var{obs_cat} and @var{thresh_obs_cat} are specified, otherwise they're empty.]
## @end itemize
##
## @subheading Example
## m = 10; n = 20; k = 3;
## obs = randn (m, n); pred_mu = obs + randn (m, n, k); pred_sigma = ones (m, n, k); pred_nu = 10 * ones(m, n, k);
## n_cats = 10; thresh_obs_cat = [-Inf; norminv((1:(n_cats-1))'/n_cats); Inf]; obs_cat = reshape(interp1(thresh_obs_cat, (1:n_cats+1)', obs(:), "previous"), [m n]);
## options.obs_cat = obs_cat; options.thresh_obs_cat = thresh_obs_cat; options.n_cats = n_cats;
## [struct_out] = sefo_verification (pred_mu, pred_sigma, pred_nu, obs, options);
## @example
## @group
## @end group
## @end example
##
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>

function [struct_out] = sefo_verification (pred_mu, pred_sigma, pred_nu, obs, options)

  [m n k] = size (pred_mu);

  #read in options or use default values

  options_default = struct('w', ones (m, 1), ...
                          'nq', 100, ...
                          'obs_cat', [], ...
                          'thresh_obs_cat', [], ...
                          'n_cats', []);

  option_names = fieldnames(options_default);
  if nargin > 4
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

  if !isempty(thresh_obs_cat)
    if isempty(n_cats)
      n_cats = max (obs_cat(:));
    endif
    thresh_nd = ndims (thresh_obs_cat);
    if thresh_nd == 2
      if size(thresh_obs_cat, 2) == 1
        thresh_obs_cat_orig = thresh_obs_cat;
        thresh_obs_cat = nan (n_cats, 2);
        thresh_obs_cat(:, 1) = thresh_obs_cat_orig(1:(end-1));
        thresh_obs_cat(:, 2) = thresh_obs_cat_orig(2:end);
      endif
      thresh_obs_cat_orig = thresh_obs_cat;
      thresh_obs_cat = nan (m, n, 2);
      thresh_obs_cat(:, :, 1) = reshape (thresh_obs_cat_orig(obs_cat, 1), m, n);
      thresh_obs_cat(:, :, 2) = reshape (thresh_obs_cat_orig(obs_cat, 2), m, n);
    endif
    nll_cats = prob_cats = nan (n_cats, k);
  else
    nll_cats = prob_cats = [];
  endif

  rmse = nll_mean = crps = predent_mean = bias = ks = kp = nan(k, 1);
  rmse_map = bias_map = nll_map = predent_map = nan(m, k);
  rmse_ts = bias_ts = nll_ts = predent_ts = nan(n, k);
  hists = nan (nq, k);

  w = w / sum(w);
  wgts = w * (ones(1, n)/n);

  for method = 1:k
    
    mu = pred_mu(:, :, method) - obs;
    sigma = pred_sigma(:, :, method);
    nu = pred_nu(:, :, method);

    rmse(method) = sqrt(mean(mu(:) .^ 2, 1, wgts(:)));
    bias(method) = mean(mu(:), 1, wgts(:));
    rmse_map(:, method) = rms(mu, 2);
    bias_map(:, method) = mean(mu, 2);
    rmse_ts(:, method) = rms(mu, 1, w); #time series of global means
    bias_ts(:, method) = mean(mu, 1, w);
    #mean negative log likelihood
    nll = -(tpdf_log(double(mu) ./ double(sigma), double(nu))) + log(double(sigma));
    nll_mean(method) = mean(nll(:), 1, wgts(:));
    nll_map(:, method) = mean(nll, 2);
    nll_ts(:, method) = mean(nll, 1, [], w);
    clear nll
    #continuous ranked probability score
    crps(method) = mean((double(sigma) .* t_crps (double(mu) ./ double(sigma), double(nu)))(:), 1, wgts(:));
    #forecast probability by observed category
    if !isempty(thresh_obs_cat)
      prob_cat = tcdf((thresh_obs_cat(:, :, 2) - pred_mu(:, :, method))./sigma, nu) - tcdf((thresh_obs_cat(:, :, 1) - pred_mu(:, :, method))./sigma, nu);
      for cat = 1:n_cats
        if any (obs_cat(:) == cat)
          prob_cats(cat, method) = mean(prob_cat(obs_cat == cat), [], wgts(obs_cat == cat));
          nll_cats(cat, method) = mean(-log(prob_cat(obs_cat == cat)), [], wgts(obs_cat == cat));
        endif
      endfor
      clear prob_cat
    endif
    #entropy of predictive distribution (should have same expectation as nll if forecasts are well calibrated)
    predent = t_entropy (double(nu)) + log(double(sigma));
    predent_mean(method) = mean (predent(:), 1, wgts(:));
    predent_map(:, method) = mean (predent, 2);
    predent_ts(:, method) = mean (predent, 1, w);
    clear predent
    #compare locations of observations to predicted CDF to uniform distribution
    cdf = tcdf (-double(mu) ./ double(sigma), double(nu));
    [pval kval] = kolmogorov_smirnov_test (cdf(:), "unif", 0, 1);
    ks(method) = kval / sqrt(n*m);
    kp(method) = pval;
    #histogram
    [hists(:, method), ~] = hist (ceil(nq*cdf(:)), 1:nq);
    clear cdf
  endfor

  nll = nll_mean;
  predent = predent_mean;
  save_names = {'rmse', 'bias', 'nll', 'crps', 'predent', 'ks', 'kp', 'rmse_map', 'bias_map', 'nll_map', 'predent_map', ...
    'rmse_ts', 'bias_ts', 'nll_ts', 'predent_ts', 'hists', 'nll_cats', 'prob_cats'};
  for i = 1:numel(save_names)
    if exist(save_names{i}, "var")
      eval(['struct_out.' save_names{i} ' = ' save_names{i} ';']);
    endif
  endfor

endfunction
