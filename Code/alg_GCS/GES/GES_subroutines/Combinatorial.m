% sub = Combinatorial (T0); % find all the sbusets of T0
function sub = Combinatorial (T0)
sub={};
count=0;
if(isempty(T0))
    count=count+1;
    sub{count}=zeros(1,0)*[];%a 1x0 empty matrix
else
    if(length(T0)==1)
        sub{1}=zeros(1,0)*[];
        sub{2}=T0; % when T0 is a scale, it is a special case!!
    else
        for n=0:length(T0)
            com=nchoosek(T0,n);
            for i=1:size(com,1)
                count=count+1;
                sub{count}=com(i,:);
            end
        end
    end
end