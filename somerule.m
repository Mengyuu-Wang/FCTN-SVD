function [R] = somerule(rc,i,n)
    R = ones(1,n-i);
    temp = factor(rc);
    temp = sort(temp,"descend");
    if length(temp) <= n-i
        R(1:length(temp)) = temp;
    else
        R(1:n-i) = temp(1:n-i);
        temp = temp(n-i+1:end);
        while length(temp) > n-i
            R(1:n-i) = R(1:n-i) .* temp(1:n-i);
            temp = temp(n-i+1:end);
        end
        R(1:length(temp)) = R(1:length(temp)) .* temp;
    end
end


%temp = factor(rc);
%        temp = sort(temp,"descend");
%        if length(temp) <= n-1
%            R(1,2:length(temp)+1) = temp;
%        else
%            R(1,2:n-1) =  temp(1:n-2);
%            R(1,n) = prod(temp(n-1:end));
%        end

%temp = factor(rc);
%        temp = sort(temp,"descend");
%        if length(temp) <= n-i
%            R(i,i+1:length(temp)+i) = temp;
%        else
%            R(i,i+1:n-1) = temp(1:n-i-1);
%            R(i,n) = prod(temp(n-i:end));
%        end


