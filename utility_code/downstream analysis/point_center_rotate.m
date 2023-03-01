function rotpoint=point_center_rotate(pt,theta)

bbox=[min(pt,[],1),max(pt,[],1)-min(pt,[],1)];
cen=bbox(1:2)+bbox(3:4)/2;

pt=pt-cen;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
rotpoint = (R*pt')';
rotpoint=rotpoint+cen;
rotpoint=rotpoint-min(rotpoint,[],1);
