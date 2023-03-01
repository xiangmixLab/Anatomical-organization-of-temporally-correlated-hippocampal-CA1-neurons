
function velo_rate1=velo_rate_012021(velo,Fs,thresh,win_leng)

velo_rate=double(velo>thresh);
velo_rate1=velo_rate*0;
rate_win=win_leng*Fs; %10 sec window
for i=rate_win+1:length(velo_rate1)-rate_win
    velo_rate1(i)=sum(velo_rate(round(i-rate_win/2):round(i+rate_win/2)))/10; % (times_moving)/sec
end