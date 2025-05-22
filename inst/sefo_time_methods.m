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

##run various prediction methods and record computing speed of each

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>





function file = sefo_time_methods (predict_year, predict_month, lag, opts)

  data_dir = getfield (opts, 'data_dir');
  re_process = getfield (opts, 're_process');
  var_name = getfield (opts, 'var_name');
  obs_dataset = getfield (opts, 'obs_dataset');
  methods = getfield (opts, 'methods');
    
  n_methods = numel(methods);

  
  file = [data_dir 'time_methods_' var_name '_' obs_dataset '_' num2str(predict_year) '_' num2str(predict_month, '%2.2i') '_lag_' num2str(lag, '%2.2i') '.mat'];

  if ~exist(file, "file") || re_process
  
      timings = nan (n_methods, 1);  
      for m = 1:n_methods
        method = methods{m};
        tic;        
          file_predict = sefo_predict (predict_year, predict_month, lag, method, opts);
        timings(m) = toc;
     endfor
    
    save (file, "methods", "timings")

  endif
  
endfunction

%{
methods{timings > 1E3}
ans = ewma_linear_model_reg
ans = ewma_linear_multimodel_reg
ans = model_linear_reg
ans = multimodel_linear_reg


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
%!                 'p_spline', 1E-4, ...
%!                 'dof_spline', 2.5, ...
%!                 'lamda', 10, ...
%!                 'predict_adj', false, ...
%!                 'predict_adj_years', 1994:1998, ...
%!                 'predict_adj_months', 1:12, ...                
%!                 'methods', {{'history', 'linear', 'quadratic', 'ma', 'ma_globcorr', 'ewma', 'ewma_globcorr', 'ewma_linear', 'ewma_linear_globcorr', 'ewma_model', 'ewma_linear_model_reg', 'ewma_linear_model', 'model_linear_reg', 'multimodel_linear_reg', 'rawmodel', 'rawmodel_scaled', 'rawmodel_globscaled', 'rawmultimodel_glob', 'rawmultimodel_glob_linear', 'rawmultimodel', 'rawmodel_linear', 'rawmodel_globcorr', 'trend', 'model_trend', 'multimodel_trend', 'globalmultimodel_trend', 'regmultimodel_trend'}}, ... 
%!                 'clim_years', 1981:2010, ...
%!                 'gcms', {{'CMC1-CanCM3', 'CMC2-CanCM4', 'GFDL-CM2p1-aer04', 'COLA-RSMAS-CCSM3', 'NASA-GMAO-062012'}} ... 
%!                 );               
%! predict_year = 2000; predict_month = 1; lag = 1;
%! file = sefo_time_methods (predict_year, predict_month, lag, opts)
%! load(file)
%! methods, timings


