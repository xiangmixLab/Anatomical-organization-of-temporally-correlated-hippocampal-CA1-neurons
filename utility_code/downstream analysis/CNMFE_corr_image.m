function Cn=CNMFE_corr_image(Y)
    
%% reshape input
    [d1,d2,~]=size(Y);
    HY = reshape(Y, d1*d2, []);

%% substract median and get noise std
    HY = bsxfun(@minus, HY, median(HY, 2));
    Ysig = GetSn(HY);

%% denoise?
    HY_thr = HY;
    sig=3;
    HY_thr(bsxfun(@lt, HY_thr, Ysig*sig)) = 0;

%% CNMFE-style autocorr image
% the correlation is the multiplication of a pixel and its surrounding pixels, avged by pixel num
    Cn = correlation_image(HY_thr, [1,2], d1,d2); 