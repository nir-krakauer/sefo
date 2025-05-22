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

%hack to produce a scaled image plot using imagesc, with the difference that NaN and infinity values are plotted as white, not as the first color in the default color map
%
%works by setting white as the first color in the color map, and setting NaN values equal to below the minimum value in the array
%
%the function bluered sets the color map to a blue-red range that dichotomizes the top and bottom halves of the range better
%
%inelegant in that colormap will show the expanded color range
%
%call with same arguments as imagesc

## Author: Nir Krakauer <nkrakauer@ccny.cuny.edu>




function H = imagesc_withnan(varargin)

%attempt to find what the arguments passed are -- in particular, which argument is the data array, which we'll call A, and whether a range to plot is given
pos_limits = 0;
if (nargin == 1)
	pos = 1;
elseif (nargin == 2)
	if isvector(varargin{2}) && (numel(varargin{2}) == 2)
		pos = 1;	
		pos_limits = 2;
	else
		pos = 2;
	end
elseif (nargin == 3)
	if isvector(varargin{3}) && (numel(varargin{3}) == 2)
		pos = 2;
		pos_limits = 3;
	else
		pos = 3;
	end
elseif (nargin == 4)
	pos = 3;
	if (numel(varargin{4}) == 2)
		pos_limits = 4;
	end
else	
	error('too many input arguments (4 maximum -- see imagesc help)')
end

A = varargin{pos};

%determine the maximum and minimum color map step that will be used in imagesc
h = colormap;
ncolors = size(h, 1); %number of levels in color map

if pos_limits
	maxval = varargin{pos_limits}(2);
	minval = varargin{pos_limits}(1);
else
	maxval = max(A(:));
	minval = min(A(:));
end
minval_thresh = minval + 0.5*(maxval-minval)/ncolors;
A(A < minval_thresh) = minval_thresh; %ensure that finite values that are below the range are plotted as the first color, not as white
#A(A < minval) = minval;
#A(A > maxval) = maxval;
#A = 1.5 + round((A-minval) * (ncolors-1)/(maxval-minval));

if any(isnan(A(:)))
	%modify the color map to include white for NaN values
	nanval = minval - (maxval - minval)/ncolors;
	A(~isfinite(A)) = nanval;
	#A(~isfinite(A)) = 0;
	h = [[1 1 1]; h]; %add white at beginning of colormap for NaN values
	colormap(h);

	%insert the redone data array into the input arguments for imagesc
	varargin{pos} = A;
	if pos_limits
		varargin{pos_limits} = [nanval maxval];
		#varargin{pos_limits} = [0 ncolors];
		disp(varargin{pos_limits})
	end
end

%finally, we're ready to call imagesc
H = imagesc(varargin{:});

