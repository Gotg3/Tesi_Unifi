function limitIndexes=findLimits(sgn,thresholds)          
[~,m]=max(sgn);
limitIndexes=zeros(1,2);
% for j=1:length(thresholds)
    for i=m:-1:1
        if sgn(i)>thresholds
            limitIndexes(1,1)=i;
        else
           
        end
          
    end
    for i=m:length(sgn)
        if sgn(i)>thresholds
            limitIndexes(1,2)=i;
        else
           
        end
       
    end
% end
