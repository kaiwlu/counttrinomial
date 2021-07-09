disc_data = readmatrix('discrepancy_data.csv');

figure;
plot3(disc_data(:,1),disc_data(:,3),disc_data(:,4),'.');
title('Discrepancy of cx^d+x over F_p');
xlabel('p (prime)');
ylabel('d (degree)');
zlabel('discrepancy of value sequence');