% S = find_subset_include(sub(k),sub); %  find those subsets that include sub(k)
function Idx = find_subset_include (s0,sub)

if(isempty(s0) | isempty(sub))
    Idx = ones(1,length(sub));
else
    Idx = zeros(1,length(sub));
    for i=1:length(sub)
        tmp = intersect(s0,sub{i}) ;
        if(~isempty(tmp))
            if(isequal(tmp, s0))
                Idx(i)=1;
            end
        end
    end
end

