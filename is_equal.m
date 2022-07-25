function flag=is_equal(a,b)

% Compare two values to check if they are equal or not
%
% Return value "1" if "a" equals "b", otherwise return value "0"

    if a==b
        flag=1;
    end
    if a~=b
        flag=0;
    end
end