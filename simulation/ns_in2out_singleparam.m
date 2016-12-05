function F = ns_in2out_singleparam(b,P,x)

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

% sqrt option:
% F= P - (b(1) .* x.^b(2));

F= P - ( b(1) .* (log10((b(2)+x)./b(2))) );


% invert is: x = b(2).*(10.^(y./bb(1)) - 1);
