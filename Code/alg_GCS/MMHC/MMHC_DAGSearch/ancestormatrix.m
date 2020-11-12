function AM = ancestormatrix(G)
%ANCESTORMATRIX Ancestor matrix.
%   AM = ORG.MENSXMACHINA.GRAPH.ANCESTORMATRIX(G) returns the ancestor
%   matrix of directed acyclic graph (DAG) G. G is a M-by-M sparse matrix
%   representing a directed acyclic graph. Each nonzero element in G
%   denotes an edge in the graph. AM is a M-by-M sparse matrix. Each
%   nonzero element in AM denotes an ancestor-descendant relationship.
%
%   ORG.MENSXMACHINA.GRAPH.ANCESTORMATRIX does not check if G is acyclic.
%
%   Example:
%
%       import org.mensxmachina.graph.ancestormatrix;
% 
%       % create a DAG
%       G = sparse([1 1 2 3], [2 3 4 4], ones(1, 4), 4, 4);
%
%       % view DAG
%       view(biograph(G));    
% 
%       % create ancestor matrix
%       am = ancestormatrix(G);
%
%       % view matrix
%       view(biograph(am)); 
%
%   See also ORG.MENSXMACHINA.GRAPH.ANCESTORS.

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

% References:
% [1] Giudici,P. and Castelo,R. (2003) Improving Markov chain Monte Carlo
%     model search for data mining. Machine Learning, 50, 127-158.

import ancestors;

assert(issparse(G) && size(G, 1) == size(G, 2));

nn = size(G, 1);

% initialize ancestor graph
AM_i = [];
AM_j = [];

for j = 1:nn
    
    ancestors_j = ancestors(G, j);
    
    AM_i = [AM_i ancestors_j];
    AM_j = [AM_j repmat(j, size(ancestors_j))];
    
end

AM = sparse(AM_i, AM_j, true, nn, nn);

end