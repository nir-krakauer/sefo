## Copyright (C) 2014-2016 Nir Krakauer
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

%given the month and year, return the number of days in the month


function n_days = days_per_month(month, year)

	dpm = [31 28 31 30 31 30 31 31 30 31 30 31]';
	
	n_days = dpm(month) + (month == 2) .* is_leap_year(year);
