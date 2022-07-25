function sb=sigma_slope(b,data)

% Re-calculate the error of slope by equation (17) in Mahon, 1996
%
% References: 
% Mahon, Keith I. "The New “York” regression: Application of an improved 
%   statistical method to geochemistry." International Geology Review 38.4 
%   (1996): 293-303. (https://www.tandfonline.com/doi/abs/10.1080/
%   00206819709465336

    X=data.x;
    sX=data.dx;
    Y=data.y;
    sY=data.dy;
    w=1./(sY.*sY+b.*b.*sX.*sX);
    meanX=sum(w.*X)/sum(w);
    meanY=sum(w.*Y)/sum(w);
    u=X-meanX;
    v=Y-meanY;
    m=calculate_m(w,u,v,sX,sY,b);
    n=calculate_n(w,u,v,sX,sY,b);
    k=sum( w.*w.*(2.*b.*u.*v.*sX.*sX + (u.*u.*sY.*sY-v.*v.*sX.*sX) ) )+...
        4.*sum(-w.*w.*w.*b.*sX.*sX.* (b.*b.*u.*v.*sX.*sX + b.* (u.*u.*...
        sY.*sY-v.*v.*sX.*sX) - u.*v.*sY.*sY) );
    sb=sqrt( sum(m.*m.*sX.*sX + n.*n.*sY.*sY) ./k./k);
end