function m=calculate_m(w,u,v,sX,sY,b)

% Calculate equation (20) in Mahon, 1996
%
% References: 
% Mahon, Keith I. "The New “York” regression: Application of an improved 
%   statistical method to geochemistry." International Geology Review 38.4 
%   (1996): 293-303. (https://www.tandfonline.com/doi/abs/10.1080/
%   00206819709465336

count=size(sX,2);
m=zeros(1,count);
for i=1:count
    m(i)=0;
    for j=1:count
        m(i)=m(i)+w(j).*w(j).*(is_equal(i,j)-w(i)/sum(w)).*(b*b.*v(j).*...
            sX(j).*sX(j)+2*b.*u(j).*sY(j).*sY(j)-v(j).*sY(j).*sY(j));
    end
end
end