% ------------------------------------------------------------------------------
% Create random sample accordibng to given distribution.
%
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function y = randdf(size1,df1,flag)
    if numel(size1) == 1
        n=1;
        m=size1(1);
    elseif numel(size1) == 2
        n=size1(1);m=size1(2);
    else
        return
    end
    
    all=n*m;
    if size(df1,1) ~= 2
        return
    end
    
    if sum(df1(1,:)<0) > 0
        return
    end
    
    if size(df1,2) < 2
        return
    end
    
    if strcmp(flag,'pdf')
        df1(1,:) = cumsum(df1(1,:))/sum(df1(1,:));
    elseif strcmp(flag,'cdf')
        if sum(diff(df1(1,:))<0)>0
            return
        end
        
        df1(1,:) = df1(1,:)/df1(1,end);
    else
        return
    end
    df1(1,:) = df1(1,:) + [1:size(df1,2)]*eps;
    temp = rand(1,all);
    temp = interp1(df1(1,:),df1(2,:),temp,'linear','extrap');
    y = reshape(temp,n,m);
