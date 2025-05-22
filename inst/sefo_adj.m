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

##fit multiplicative adjustment factors for the spread and dof parameters in the prediction t distribution based on flattening the histogram of observation CDF locations

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>




function file = sefo_adj (years, months, lag, method, opts)

  data_dir = getfield (opts, 'data_dir');
  re_download = getfield (opts, 're_download');
  re_process = getfield (opts, 're_process');  
  var_name = getfield (opts, 'var_name');  
  obs_dataset = getfield (opts, 'obs_dataset');

  n_months_combine = numel (lag);
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
    #error ("the inputs months, years are not of common size")
  endif
  if nt == 1
    dates = [num2str(years(1)) '_' num2str(months(1), '%2.2i')];
  else
    dates = [num2str(years(1)) '_' num2str(months(1), '%2.2i') '_' num2str(years(end)) '_' num2str(months(end), '%2.2i')];
  endif
  
  file = [data_dir 'adj_' var_name '_' obs_dataset '_' method '_' dates '_lag' num2str(lag, '_%2.2i') '.mat'];

  if ~exist(file, "file") || re_process
    
    #obtain all requested months and methods
    predict_mu = predict_sigma = predict_nu = obss = area_wgtss = [];
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
        if exist("obs_inds", "var") && ~isempty(obs_inds)    
          obs += w(m)*arr_regrid(obs_inds);
        else
          obs += w(m)*arr_regrid(:);
        endif
        clear arr_regrid
      endfor


      file_predict = sefo_predict (predict_year, predict_month, lag, method, opts);
      load (file_predict)
      predict_mu = [predict_mu mu(:)'];
      predict_sigma = [predict_sigma sigma(:)'];
      if (isscalar (nu));
        nu = nu * ones(numel(mu), 1);
      endif
      predict_nu = [predict_nu nu(:)'];
      area_wgtss = [area_wgtss area_wgts(:)'];
      obss = [obss obs(:)'];      
    endfor
    predict_mu = predict_mu';
    predict_sigma = predict_sigma';
    predict_nu = predict_nu';
    area_wgts = area_wgtss'; clear area_wgtss
    obs = obss'; #clear obss
    
    np = size (predict_mu, 1);
    nq = 100;
    
    area_wgts /= mean(area_wgts);
    
    #optimal adjustment factors
    hist_misfit = @(facs) sumsq(hist(ceil(nq*tcdf((obs - predict_mu) ./ (exp(facs(1)) * predict_sigma), predict_nu/exp(facs(2)))), 1:nq) - np/nq);
    [facs_optim, fval] = fminsearch(hist_misfit, [0 0.3]);

    
    save (file, "facs_optim")

  endif
  
endfunction

%{
rmse, nll_mean, ks, plot(hists) #general underestimation of frequency of extremes

m = 8; fac = 1.15;
nll_new = -(tpdf_log((obs - predict_mu(:, m)) ./ (fac * predict_sigma(:, m)), predict_nu(:, m)));
nll_mean_new = mean(nll_new .* area_wgts)';
cdf_new = tcdf((obs - predict_mu(:, m)) ./ (fac * predict_sigma(:, m)), predict_nu(:, m));
[pval kval] = kolmogorov_smirnov_test (cdf_new, "unif", 0, 1);
ks_new = kval / sqrt(numel(cdf_new));
nq = 100; [hists_new, ~] = hist (ceil(nq*cdf_new), 1:nq);

m = 8; fac = 1; dfm = 1;
nll_new = -(tpdf_log((obs - predict_mu(:, m)) ./ (fac * predict_sigma(:, m)), predict_nu(:, m)/dfm));
nll_mean_new = mean(nll_new .* area_wgts)';
cdf_new = tcdf((obs - predict_mu(:, m)) ./ (fac * predict_sigma(:, m)), predict_nu(:, m)/dfm);
[pval kval] = kolmogorov_smirnov_test (cdf_new, "unif", 0, 1);
ks_new = kval / sqrt(numel(cdf_new));
nq = 100; [hists_new, ~] = hist (ceil(nq*cdf_new), 1:nq);
nll_mean_new, ks_new, plot(hists_new)

np = numel(obs);

#try to fit adjustment parameters to flatten the rank histogram
np = numel(obs); nq = 100; m = 8;
hist_misfit = @(facs) sumsq(hist(ceil(nq*tcdf((obs - predict_mu(:, m)) ./ (exp(facs(1)) * predict_sigma(:, m)), predict_nu(:, m)/exp(facs(2)))), 1:nq) - np/nq);
[facs_optim, fval] = fminsearch(hist_misfit, [0 0.3], optimset("Display", "iter"));
  f = exp(facs_optim); m = 8; fac = f(1); dfm = f(2);
#NCEP 2000-2014; lag 2
#[0 0.3]: -7.9037e+09    -2.3616e+09
#0.0505150823595708   0.2485224721352517
#1.05181272705312   1.28212963423613
#1.052, 1.282: 1.4800, 0.018474 (some dipole near 0.5) 
#BEST 2000-2014; lag 2
#values: [0 0]: 1.0088e+10, [0 0.3]: 5.0573e+09
#0.0210853944751948   0.6130409552355541, 464929348
#exp(): 1.02130926208228   1.84603658652249
#BEST 1994-1998; lag 2
#1.5914e+09, 1.0077e+09
#0.0467334493871906   0.3113666221506303
#1.04784266874666   1.36528967487869
#1.048, 1.365: 1.5437, 0.025662 [more high extremes, less low extremes though]
#1, 1: 1.5956, 0.034362 [lots more high extremes, slightly more low extremes]
#BEST 2000-2014; lag 1
#6.4589e+09 4.3405e+08
#0.0256450667859212   0.6362132023858944
#1.02597673062219   1.88931287050099
#1.026, 1.889: 1.50707, 0.0034523
#1.049, 1.402: 1.48513, 0.0043146
#BEST 1994-1998; lag 1
#1.6136e+09 9.4143e+08
#0.0480739119304260   0.3377673742791326
#1.04924820441882   1.40181436739524
#1.049, 1.402: 1.54727, 0.023254 [not fully corrected]
#1, 1: 1.60111, 0.034044

NCEP 2000-2014; lag 2
rmse
   1.5450
   1.4339
   1.4343
   1.4403
   1.4139
   1.4138
   1.5034
   1.3981
   1.3992
nll_mean
   1.5911
   1.4818
   1.5111
   1.4894
   1.5505
   1.5581
   1.5602
   1.5345
   1.5504
ks
   0.181298
   0.053433
   0.053445
   0.017149
   0.032864
   0.030473
   0.029228
   0.025691
   0.029242

BEST 2000-2014; lag 2
rmse
   1.2502
   1.1808
   1.1686
   1.2104
   1.1874
   1.1792
   1.2533
   1.1661 -- still underestimates extremes
   1.1671
nll_mean
   1.5560
   1.4405
   1.4955
   1.4829
   1.5093
   1.5513
   1.5511
   1.5274
   1.5434
ks
   0.183053
   0.049832
   0.044637
   0.017167
   0.015137
   0.015626
   0.014124
   0.011429
   0.014965
BEST 1994-1998; lag 1
   1.2197
   1.2266
   1.2010
   1.3174
   1.2357
   1.2111
   1.3497
   1.1931
   1.1951   
nll_mean
   1.6193
   1.5540
   1.5209
   1.5262
   1.6332
   1.6166
   1.6326
   1.6011
   1.6327
ks
   0.1194629
   0.0556438
   0.0313242
   0.0091301
   0.0532904
   0.0366622
   0.0333513
   0.0340436
   0.0420969   
BEST 2000-2014; lag 1   
   1.2502
   1.1808
   1.2098
   1.2639
   1.1874
   1.1647
   1.2392
   1.1644
   1.1637
nll_mean
   1.5560
   1.4405
   1.4964
   1.4870
   1.5093
   1.5548
   1.5519
   1.5381
   1.5530
ks
   0.183053
   0.049832
   0.043752
   0.012088
   0.015137
   0.014802
   0.013757
   0.012225
   0.016384
   
   
m = 8
  fac: nll_mean_new, ks_new, plot(hists_new)
  1.1: 1.4307, 0.014686 -- too little at edges (except at extremes)
  1.05: 1.4759, 0.0062853 -- same, but less so
  1.15: 1.3908, 0.024534 -- strong inverse U, extremes also damped
  conclusion: shape is wrong (tails not fat enough), increasing SD insufficient
fac 1, dfm
  1.15: 1.5263, 0.010794 -- tails still too strong
  1.5: 1.5244, 0.0093573 -- not much improved
  2: 1.5235, 0.0078588 -- esp. lower extremes are flattened
  3: 1.5255, 0.0058174 -- now reversed
  2.5: 1.5240, 0.0066090 -- still reversed
1.1, 2: 1.4334, 0.018895 -- reversed
1.1, 1.5: 1.4314, 0.016664 -- still reversed
1.05, 1.5: 1.4750, 0.0071112 -- still somewhat reversed, except for extremes
1.04, 2: 1.4848, 0.0067042 -- now reversed at extremes
1.021, 1.846: 1.5026, 0.0053883 -- not great, but quantitatively pretty good (% mismatches small; high 99 percentile, low 1)
1, 1: 1.5274, 0.011429
1.048, 1.365: 1.4771, 0.0065489 -- actually not bad -- so, parameters are tranferable (optimized 1993-1998 for 2000-2014
%}

%!demo
%! if !exist("data_dir", "var")
%!   data_dir = [pwd '/'];
%! endif
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
%!                 'methods', {{'history', 'linear', 'rawmodel', 'rawmodel_linear', 'trend', 'model_trend', 'multimodel_trend', 'globalmultimodel_trend', 'regmultimodel_trend'}}, ...   
%!                 'gcms', {{'CMC1-CanCM3', 'CMC2-CanCM4', 'GFDL-CM2p1-aer04', 'COLA-RSMAS-CCSM3', 'NASA-GMAO-062012'}} ...       
%!                 );
%! years = 2006; months = 6; lag = 1; method = 'globalmultimodel_trend';
%! file = sefo_adj (years, months, lag, method, opts); load(file)    
%! facs_optim   
%! #-0.086882   0.629417
