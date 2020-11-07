function sign = check_restrict(G,Mask)
sign = 1;
for i = 1:size(G,1)
    for j = 1:size(G,2)
        if(Mask(i,j)==0 & G(i,j)==1)
            sign = 0;
            break;
        end
    end
    if(~sign)
        break;
    end
end