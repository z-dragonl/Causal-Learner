classdef(Sealed) citperformeddata < event.EventData
%CITPERFORMEDDATA Conditional-independence-test-performed data.
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.CIT.CITPERFORMEDDATA is the class
%   of conditional-independence-test-performed data.
%   Conditional-independence-test-performed data are passed to listeners of
%   citPerformed events.

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
j % linear index of variable X -- a numeric positive integer
k % linear indices of set Z of variables -- a numeric row vector of positive integers    
    
pValue % p-value -- either NaN or a numeric scalar in [0, 1]
stat % statistic -- either NaN or a numeric real scalar

end

methods

% constructor

function Obj = citperformeddata(i, j, k, p, t)
%CITPERFORMEDDATA Create conditional-independence-test-performed data.
%   OBJ = ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.CIT.CITPERFORMEDDATA(I, J, K,
%   S, P, T) creates conditional-independence-test-performed data with
%   linear index I and J of variable X and Y, respectivelly, linear indices
%   K of set Z of variables, p-value P, and statistic T. I and J are
%   numeric positive integes. K is a numeric row vector of positive
%   integers. S is a numeric integer in range [1, 3]. P is a numeric scalar
%   in [0, 1]. T is a numeric real scalar.
    
    % (no validation)
    
    % set properties
    Obj.i = i;
    Obj.j = j;
    Obj.k = k;
    Obj.pValue = p;
    Obj.stat = t;

end

end

end