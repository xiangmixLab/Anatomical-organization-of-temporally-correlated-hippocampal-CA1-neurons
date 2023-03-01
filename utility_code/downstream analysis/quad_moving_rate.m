
function quad_move_rate1=velo_rate_012021(quad,Fs,thresh)

quad_velo=abs(diff(double(quad)));
quad_move_rate=double(quad_velo>1);
quad_move_rate1=quad_move_rate*0;
rate_win=10*Fs; %10 sec window
for i=rate_win+1:length(quad)-rate_win
    quad_move_rate1(i)=sum(quad_move_rate(round(i-rate_win/2):round(i+rate_win/2)))/10; % (times_moving)/sec
end