val_data = readmatrix('a0_valueset.csv');

figure;
plot3(val_data(:,1),val_data(:,2),val_data(:,4),'.');
title('Size of value set');
xlabel('p (prime)');
ylabel('d (degree)');
zlabel('size of value set');

figure;
plot(gcd(val_data(:,1)-1,val_data(:,2)-1),val_data(:,4),'.');
title('gcd(d-1,p-1) vs size of value set');
xlabel('gcd(d-1,p-1)');
ylabel('size of value set');