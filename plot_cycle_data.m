% cyclen = readmatrix('avg_cyc_len.csv');
% cycmax = readmatrix('avg_cyc_max.csv');
% cycmin = readmatrix('avg_cyc_min.csv');
cycnum = readmatrix('avg_cyc_num.csv');

% for g = 2:10
%     figure;
%     [m,n]=size(cyclen);
%     halfp = false(m,1);
%     for k = 1:(g-1)
%         halfp = halfp | (((cyclen(:,1)*k+g-k)/g == cyclen(:,2)) & gcd(k,g-k) == 1);
%     end
%     plot(cyclen(halfp,1), cyclen(halfp,3), '.')
%     title(strcat('Average |Per(f)| with H_p(d)=', num2str(g)));
%     xlabel('p (prime)');
%     ylabel('Average |Per(f)|');
% end


% figure;
% halfp = (cyclen(:,1)+2)/3 == cyclen(:,2);
% plot(cyclen(halfp,1), cyclen(halfp,3))
% title('Average sum of cycle length of cx^{(p+2)/3}+x+a');
% xlabel('p (prime)');
% ylabel('d (degree)');
% zlabel('Average sum of cycle length');

% figure;
% halfp = (cycmax(:,1)+2)/3 == cycmax(:,2);
% plot(cycmax(halfp,1), cycmax(halfp,3))
% title('Average max of cycle length of cx^{(p+2)/3}+x+a');
% xlabel('p (prime)');
% ylabel('d (degree)');
% zlabel('Average max of cycle length');
% figure;
% 
% halfp = (cycmin(:,1)+2)/3 == cycmin(:,2);
% cycmin(halfp,1)
% cycmin(halfp,3)
% plot(cycmin(halfp,1), cycmin(halfp,3))
% title('Average min of cycle length of cx^{(p+2)/3}+x+a');
% xlabel('p (prime)');
% ylabel('d (degree)');
% zlabel('Average min of cycle length');
% figure;
% 
% halfp = (cycnum(:,1)+2)/3 == cycnum(:,2);
% plot(cycnum(halfp,1), cycnum(halfp,3))
% title('Average number of cycles of cx^{(p+2)/3}+x+a');
% xlabel('p (prime)');
% ylabel('d (degree)');
% zlabel('Average number of cycles');

% figure;
% halfp = (cyclen(:,1)+1)/2 == cyclen(:,2);
% plot(cyclen(halfp,1), cyclen(halfp,3))
% title('Average sum of cycle length of cx^{(p+1)/2}+x+a');
% xlabel('p (prime)');
% ylabel('d (degree)');
% zlabel('Average sum of cycle length');
% 
% figure;
% halfp = (cycmax(:,1)+1)/2 == cycmax(:,2);
% plot(cycmax(halfp,1), cycmax(halfp,3))
% title('Average max of cycle length of cx^{(p+1)/2}+x+a');
% xlabel('p (prime)');
% ylabel('d (degree)');
% zlabel('Average max of cycle length');
% figure;
% 
% halfp = (cycmin(:,1)+1)/2 == cycmin(:,2);
% cycmin(halfp,1)
% cycmin(halfp,3)
% plot(cycmin(halfp,1), cycmin(halfp,3))
% title('Average min of cycle length of cx^{(p+1)/2}+x+a');
% xlabel('p (prime)');
% ylabel('d (degree)');
% zlabel('Average min of cycle length');
% figure;
% 
% halfp = (cycnum(:,1)+1)/2 == cycnum(:,2);
% plot(cycnum(halfp,1), cycnum(halfp,3))
% title('Average number of cycles of cx^{(p+1)/2}+x+a');
% xlabel('p (prime)');
% ylabel('d (degree)');
% zlabel('Average number of cycles');

% figure;
% plot3(cyclen(:,1),cyclen(:,2),cyclen(:,3),'.');
% title('Average sum of cycle length of cx^d+x+a');
% xlabel('p (prime)');
% ylabel('d (degree)');
% zlabel('Average sum of cycle length');
% 
% figure;
% plot3(cycmax(:,1),cycmax(:,2),cycmax(:,3),'.');
% title('Average max of cycle length of cx^d+x+a');
% xlabel('p (prime)');
% ylabel('d (degree)');
% zlabel('Average max of cycle length');
% 
% figure;
% plot3(cycmin(:,1),cycmin(:,2),cycmin(:,3),'.');
% title('Average min of cycle length of cx^d+x+a');
% xlabel('p (prime)');
% ylabel('d (degree)');
% zlabel('Average min of cycle length');
% 
figure;
plot(gcd(cycnum(:,1)-1,cycnum(:,2)-1),cycnum(:,3),'.');
title('Average number of cycles of cx^d+x+a');
xlabel('H_p(d)');
ylabel('Average number of cycles');
