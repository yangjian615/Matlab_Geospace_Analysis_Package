function [IND] = N_nearest(x,Y,N,mode)
% INPUT: 
% x ---- a double: reference point in a vector (e.g. a time element of the primary field)
% Y ---- a column vector: the search space (e.g. secondary time field)
% N ---- a positive integer: the N nearest points in y[:] from x[i]
% mode -- for linear interpolation, and the special case N = 2, such that
%         the joined time t[j] has to fit between the two closest time elements T[i] in the search space
%         T[i] < t[j] < T[i+1]
% OUTPUT:
% ind -- a Nx1 vector of indeces of N elements in y which are nearest to an element x

numY = length(Y);
temp = Y - x*ones(numY,1); % difference
no_diff = find(temp < eps(x) & temp > -eps(x)); 


if (N == 2) && mode
    if no_diff
        IND = no_diff;
    else
        more_than_0 = temp(temp>0);
        less_than_0 = temp(temp<0);
        
        if ~isempty(more_than_0) && ~isempty(less_than_0)
            IND(1) = find(temp == -min(-less_than_0), 1); % find the index of a negative element in Y which is closest to x
            IND(2) = find(temp == min(more_than_0),1); % find the index of a positive element in Y which is closest to x
        else 
            IND = N_nearest(x,Y,N,0); % mode 0 = logical FALSE
        end
          
    end
else
    [yy,ind] = sort(abs(temp)); % sort the pointwise absolute difference between Y[:] and x,and obtain the corresponding indeces
    IND = ind(1:N);
end

end