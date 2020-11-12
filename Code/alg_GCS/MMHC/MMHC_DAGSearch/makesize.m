function x = makesize(x)
%MAKESIZE Make size vector.
%   X = ORG.MENSXMACHINA.ARRAY.MAKESIZE(X) corrects numeric row vector X of
%   nonnegative integers (or Infs) to make a size vector. A size vector has
%   length >= 2.
%
%   Example:
%
%       import org.mensxmachina.array.makesize;
%
%       makesize(zeros(1, 0))
%       makesize(2)
%       makesize([3 2 4])

% Copyright 2010-2012 Mens X Machina
% 
% This file is part of Mens X Machina Common Toolbox.
% 
% Mens X Machina Common Toolbox is free software: you can redistribute it
% and/or modify it under the terms of the GNU General Public License
% alished by the Free Software Foundation, either version 3 of the License,
% or (at your option) any later version.
% 
% Mens X Machina Common Toolbox is distributed in the hope that it will be
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with Mens X Machina Common Toolbox. If not, see
% <http://www.gnu.org/licenses/>.

% no validation

if isempty(x)
    x = [1 1];
elseif length(x) == 1
    x = [x 1];
end

end