function an = ancestors(G, i)
%ANCESTORS Find ancestors of directed acyclic graph (DAG) node.
%   AN = ORG.MENSXMACHINA.GRAPH.ANCESTORS(G, I) returns the ancestors of
%   node I in directed acyclic graph G. G is an M-by-M sparse matrix. Each
%   nonzero element in G denotes an edge in the graph. I is the linear
%   index of a node in G. AN is a column vector of node linear indices in
%   G.
%
%   ORG.MENSXMACHINA.GRAPH.ANCESTORMATRIX does not check if G is acyclic.
%
%   Example:
%
%       import org.mensxmachina.pgm.ancestors;
%       
%       % create a directed graph with 4 nodes
%       G = sparse([1 1 2 3], [2 3 4 4], ones(1, 4), 4, 4);
%
%       % view graph
%       view(biograph(G));    
%
%       % find ancestors of node #4
%       an = ancestors(G, 4)

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

import ancestors;
        
assert(issparse(G) && size(G, 1) == size(G, 2));
validateattributes(i, {'numeric'}, {'integer', 'scalar', 'positive', '<=', size(G, 2)});

% find parents
an = find(G(:, i))';

for i=1:length(an) % for each parent

    % find its ancestors
    an = [ an ancestors(G, an(i)) ];

end

if ~isempty(an)

    % keep unique ancestors
    an = unique(an);

end

end