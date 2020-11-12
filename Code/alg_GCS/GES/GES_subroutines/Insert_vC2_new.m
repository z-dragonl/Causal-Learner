
function sign = Insert_vC2_new(G,j,i,NAT) % validity test for condition 2 of Insert operator
% here G is CPDAG
% Use Depth-first-Search
start=j;
target=i;
% stack(1)=start; % initialize the stack
stack(1) = struct('value',start,'pa',[]);
sign=1;% If every semi-pathway contains a node in NAT, than sign=1;
count=1;

while(~isempty(stack))
    top=stack(1);
    stack = stack(2:end); % pop
    if(top.value==target) % if find the target, search that pathway to see whether NAT is in that pathway
        curr = top;
        ss=0;
        while(1)
            if(~isempty(curr.pa))
                if(find(NAT==curr.pa.value)) % contains a node in NAT
                    ss=1;
                    break;
                end
            else
                break;
            end
            curr=curr.pa;
        end
        if(~ss) % do not include NAT
            sign = 0;
            break;
        end
        
    else
        child = find(G(top.value,:)==1 | G(top.value,:)==-1);
        sign_child=ones(1,length(child));
        % check each child, whether it has appeared before in the same pathway
        for k=1:length(child)
            curr = top;
            while(1)
                if(~isempty(curr.pa))
                    if(curr.pa.value==child(k))
                        sign_child(k)=0; % has appeared in that path before
                        break;
                    end
                else
                    break;
                end
                curr=curr.pa;
            end
        end
        for k=1:length(sign_child)
            if(sign_child(k))
                stack = [struct('value',child(k),'pa',top), stack];  % push
            end
        end
        
    end
    
end
