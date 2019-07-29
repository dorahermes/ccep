function F = LinePlusGauss(x,P,t)

% function for fitting a line + gaussian

%     Copyright (C) 2014  D Hermes
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% F= P - (x(1)*sqrt(2*pi)*normpdf(t,x(2),x(3)));
F= P - (x(1)*sqrt(2*pi)*normpdf(t,x(2),x(3)) + x(4));

% x(1)  = a*x(3), where x(3) makes the gaussian have an amplitude of 1
% --> a = the amplitude of the curve a = x(1)/x(3)

% x1: amplitude
% x2: mean
% x3: width