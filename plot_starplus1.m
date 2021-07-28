cyc_data = csvread('500_star_plus1.csv'); % Read entire file.

figure;
plot(gcd(cyc_data(:,1)-1, cyc_data(:,2)-1),cyc_data(:,4),'.');
title('number of cycles vs H_p(d)');
xlabel('H_p(d)');
ylabel('number of cycles');