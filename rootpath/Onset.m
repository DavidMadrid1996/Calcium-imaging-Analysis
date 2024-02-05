temporal = abs(unico(1)-CALCIUMroiTS.diff_perc03.times); 
temporal2 = find(temporal==min(temporal)); 

plot([CALCIUMroiTS.diff_perc03.times(temporal2) CALCIUM.AMPLITUDE_X.(myfield1)(j)],...
    [CALCIUM.baseline_amplitude.(myfield1)(j) CALCIUM.AMPLITUDE_Y.(myfield1)(j)],'LineWidth',7)

p = polyfit([CALCIUMroiTS.diff_perc03.times(temporal2:CALCIUM.AMPLITUDE_X.(myfield1)(j)/CALCIUM.srate)],...
    [-CALCIUMroiTS.diff_perc03.data(i,temporal2:CALCIUM.AMPLITUDE_X.(myfield1)(j)/CALCIUM.srate)],1);

x1 = linspace(CALCIUMroiTS.diff_perc03.times(temporal2),CALCIUMroiTS.diff_perc03.times(int32(CALCIUM.AMPLITUDE_X.(myfield1)(j)/CALCIUM.srate)));
y1 = polyval(p,x1);
plot(x1,y1,'LineWidth',7)

x = [CALCIUMroiTS.diff_perc03.times(temporal2) CALCIUM.baseline_amplitude.(myfield1)(j)]
y = [CALCIUMroiTS.diff_perc03.times(int32(CALCIUM.AMPLITUDE_X.(myfield1)(j)/CALCIUM.srate)) CALCIUM.AMPLITUDE_Y.(myfield1)(j)]

a = [x1(1) y1(1)]
b = [x1(end) y1(end)]


v_1 = [y(1),y(2),0] - [x(1),x(2),0];
theta1 = atan(v_1(2)/v_1(1))* 180/pi
v_2 = [b(1),b(2),0] - [a(1),a(2),0];
theta2 = atan(v_2(2)/v_2(1))* 180/pi; 
theta = abs(theta1-theta2)



v_1 = (v_1);
v_2 = (v_2);

Theta = atan2(norm(cross(v_1, v_2)), dot(v_1, v_2))* 180/pi

difference = (atan((y(2)-x(2))/(y(1)-x(1))) - atan((b(2)-a(2))/(b(1)-a(1)))) * 180/pi;

