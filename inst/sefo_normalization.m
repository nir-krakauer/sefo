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
## @deftypefn {Function File} {[@var{transformed}] =} sefo_normalization (@var{data}[, @var{options}])
## Empirical transformation to standard normal distribution for arbitrary data, using ranks.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{data} (dimensions: @var{m} by @var{n})
## Data to use for computing transforms. Each row is transformed seperately.
## @item
## @var{options} (structure).
## @end itemize
##
## @subheading @var{options} fields
##
## @itemize @bullet
## @item
## @var{data_to_transform} (dimensions: @var{m} by @var{nd} or 1 by @var{nd}), data to transform. Default: @var{data}.
## @item
## @var{quantiles_to_find} (dimensions: @var{m} by @var{nq} or 1 by @var{nq}), estimate values of given quantiles in original data units.
## @item
## @var{cap_extremes} (logical) If set, set extreme values returned from interpolation at ones close to the distribution of the data. Default: true.
## @end itemize
##
## @subheading Output
##
## @itemize @bullet
## @item
## @var{transformed} (dimensions: @var{m} by @var{nd} or @var{m} by @var{nq}):
## Transformed @var{data_to_transform} or values of @var{quantiles_to_find}}.
## @end itemize
##
## @subheading Example
## m = 10; n = 20;
## data = randn (m, n);
##
## options = []; options.quantiles_to_find = (1:9)/10; [transformed] = sefo_normalization (data, options);
##
## options = []; options.data_to_transform = randn(1, 9); [transformed] = sefo_normalization (data, options);
## @group
## @end group
## @end example
##
## @end deftypefn

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>

function [transformed] = sefo_normalization (data, options)

  [m n] = size (data);

  if (nargin > 1) && isfield(options, 'quantiles_to_find') && !isempty(options.quantiles_to_find)

    data = sort(data, 2);
    qd = norminv ((1:n)/(n+1));
    qf = norminv (options.quantiles_to_find);
    transformed = interp1 (qd', data', qf, 'linear', 'extrap')';

  else
    if (nargin == 1) || !isfield(options, 'data_to_transform') || isempty(options.data_to_transform)
      transformed = norminv (ranks(data, 2) / (n + 1));
    else
      data = sort (data, 2);
      qd = norminv (ranks(data, 2) / (n + 1));
      df = options.data_to_transform;
      nd = size(df, 2);

      if m > 1
        transformed = nan (m, nd);
      endif
      if size(df, 1) == 1
        for i = 1:m
          [d, inds] = unique (data(i, :));
          if numel (inds) > 1
            transformed(i, :) = interp1 (d', qd(i, inds)', df, 'linear', 'extrap')';
          else
            transformed(i, :) = 0;
          endif
        endfor
      else
        for i = 1:m
          [d, inds] = unique (data(i, :));
          if numel (inds) > 1
            transformed(i, :) = interp1 (d', qd(i, inds)', df(i, :), 'linear', 'extrap')';
          else
            transformed(i, :) = 0;
          endif
        endfor
      endif

      if !isfield(options, 'cap_extremes') || isempty(options.cap_extremes) || options.cap_extremes
        cap_vals = [1 -1] * norminv(1/(n+nd+1));
        transformed(transformed < cap_vals(1)) = cap_vals(1);
        transformed(transformed > cap_vals(2)) = cap_vals(2);
      endif

    endif
  endif



endfunction
