classdef sepsetidentifieddata < event.EventData
%SEPSETIDENTIFIEDDATA Sepset-identified data.
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.SEPSETIDENTIFIEDDATA is the class
%   of sepset-identified data. Sepset-identified data are passed to
%   listeners of sepsetIdentified events.

% Copyright 2010-2012 Mens X Machina
% 
% This file is part of Mens X Machina Probabilistic Graphical Model
% Toolbox.
% 
% Mens X Machina Probabilistic Graphical Model Toolbox is free software:
% you can redistribute it and/or modify it under the terms of the GNU
% General Public License alished by the Free Software Foundation, either
% version 3 of the License, or (at your option) any later version.
% 
% Mens X Machina Probabilistic Graphical Model Toolbox is distributed in
% the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
% the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
% PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with Mens X Machina Probabilistic Graphical Model Toolbox. If not, see
% <http://www.gnu.org/licenses/>.

properties(SetAccess = immutable)
    
i % linear index of variable X -- a numeric positive integer
j % linear index of variable Y -- a numeric positive integer
k % linear indices of the sepset variables -- a numeric row vector of positive integers

end

methods

% constructor

function Obj = sepsetidentifieddata(i, j, k)
%SEPSETIDENTIFIEDDATA Create sepset-identified data.
%   OBJ = ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.SEPSETIDENTIFIEDDATA(I, J, K)
%   creates sepset-identified data with linear index I and J of variable X
%   and Y, respectivelly, and linear indices K of the sepset. I and J are
%   numeric positive integers. K is a numeric row vector of positive
%   integers.
    
    % (no validation)
    
    % set properties
    Obj.i = i;
    Obj.j = j;
    Obj.k = k;

end

end

end