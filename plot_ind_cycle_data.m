cyc_data = readmatrix('ind_cyc_data.csv');

figure;
plot(cyc_data(:,2),cyc_data(:,6),'.');
title('Length of cycle');
xlabel('p (prime)');
ylabel('d (degree)');

figure;
plot(cyc_data(:,2),cyc_data(:,5),'.');
title('Roots of cycle');
xlabel('d (degree)');
ylabel('number of roots');