function F = ns_in2out(b,P,x,bb_in)

% function for fitting stuff

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

% this is a multidimensional function, with x and bb_in being vectors

% polynomial option:
% F= P - (b(4)./(1+bb_in) .* (b(1).*x.^2 + b(2).*x + b(3)));

% sqrt option:
% F= P - (b(1)./(b(2)+bb_in) .* x.^b(3));

% sqrt option:
F= P - (b(1)./(1+bb_in) .* x.^b(2) + bb_in*b(3));
