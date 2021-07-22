disc_data = readmatrix('discrepancy_data.csv');
n0disc_data = readmatrix('discrepancy_data_n0.csv');

figure;
plot3(disc_data(:,1),disc_data(:,3),disc_data(:,4),'.');
title('Discrepancy of cx^d+x over F_p');
xlabel('p (prime)');
ylabel('d (degree)');
zlabel('discrepancy of value sequence');

figure;
plot3(n0disc_data(:,1),n0disc_data(:,3),n0disc_data(:,4),'.');
title('Discrepancy of cx^d+x over F_p excluding 0');
xlabel('p (prime)');
ylabel('d (degree)');
zlabel('discrepancy of value sequence');

figure;
plot(disc_data(:,1),disc_data(:,4),'.');
title('Discrepancy of cx^d+x over F_p');
xlabel('p (prime)');
ylabel('discrepancy of value sequence');

figure;
plot(n0disc_data(:,1),n0disc_data(:,4),'.');
title('Discrepancy of cx^d+x over F_p excluding 0');
xlabel('p (prime)');
ylabel('discrepancy of value sequence');