% Approximating in increments of 10^k where x is an decided below
format long
n = power(10, 9);
array_size = power(10, 5);

x = linspace(0, n);
y = power(1 + (1./x), x);

x2 = [0, n];
y2 = [exp(1), exp(1)];

figure;plot(x, y, x2, y2);
axis([0 n 2.7182814 2.7182822]);
xlabel('N','fontsize',16);
legend({'Approximation', 'Actual'},'fontsize',16,'Location','northwest');
title('Fig. 1a - Approximation vs Actual');

actual = exp(1);
absolute_error = abs(actual - y);
relative_error = absolute_error / actual;

figure;plot(x, absolute_error, x, relative_error);
axis([0 n 0 0.3*power(10, -6)]);
xlabel('N','fontsize',16);
legend({'Absolute Error', 'Relative Error'},'fontsize',16,'Location','northwest');
title('Fig. 1b - Approximation vs Actual Errors');

% Approximating in increments of 2^k
format long

% Fig 2a
k_max = 53;
k = 0:1:k_max;

x = power(2,k);
y = power(1 + (1./x), x);

x2 = 0;
y2 = [exp(1), exp(1)];
figure;plot(x, y);
xlabel('N','fontsize',16);
legend({'Approximation'},'fontsize',16,'Location','northeast');
title('Fig. 2a - Underflow error as k approaches 53');

% Fig 2b
k_max = 51;
k = 0:1:k_max;

x = power(2,k);
y = power(1 + (1./x), x);

x2 = [0, power(2, k_max-1)];
y2 = [exp(1), exp(1)];

figure;plot(x, y, x2, y2);
axis([0 power(10,15) 2.718281828459 2.7182818284591]);
xlabel('N','fontsize',16);
legend({'Approximation', 'Actual'},'fontsize',16,'Location','northwest');
title('Fig. 2b - Approximation vs Actual');

% Error plot

actual = exp(1);
absolute_error = abs(actual - y);
relative_error = absolute_error / actual;

figure;plot(x, absolute_error, x, relative_error);
axis([0 power(2, k_max-1) 0 .8*power(10, -14)]);
xlabel('N','fontsize',16);
legend({'Absolute Error', 'Relative Error'},'fontsize',16,'Location','northeast');
title('Fig. 2c - Approximation vs Actual Errors');



