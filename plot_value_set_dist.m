val_data = readmatrix('a0_valueset_dist_n0.csv');

figure;
lhp = plot3(val_data(:,1),val_data(:,2),val_data(:,4), '.', 'DisplayName', 'Lower half');
lhp.Color = '#0072BD';
hold on;
uhp = plot3(val_data(:,1),val_data(:,2),val_data(:,5), '.', 'DisplayName', 'Upper half');
uhp.Color = '#EDB120';
title('Fraction of value set in each half');
xlabel('p (prime)');
ylabel('d (degree)');
zlabel('fraction of value set');
legend;
hold off;

figure;
uhp = plot(val_data(:,1),val_data(:,5), '.', 'DisplayName', 'Upper half');
uhp.Color = '#EDB120';
hold on;
lhp = plot(val_data(:,1),val_data(:,4), '.', 'DisplayName', 'Lower half');
lhp.Color = '#0072BD';
title('Fraction of value set in each half');
xlabel('p (prime)');
ylabel('fraction of value set');
legend;
hold off;

figure;
lhp = plot(val_data(:,2),val_data(:,4), '.', 'DisplayName', 'Lower half');
lhp.Color = '#0072BD';
hold on;
uhp = plot(val_data(:,2),val_data(:,5), '.', 'DisplayName', 'Upper half');
uhp.Color = '#EDB120';
title('Fraction of value set in each half');
xlabel('d (degree)');
ylabel('fraction of value set');
legend;
hold off;
