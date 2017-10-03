format long
trials = 5; % Number of trials

N_min = 0; % Range of N
N_max = 10000;
span = 100; % Span between N
arr_size = (N_max / span) + 1 % Determine arr size needed

N_range = N_min:span:N_max
arr_size
N_range(1)

% For determining condition of A
A_condition = zeros(trials, 1);
A_mean_cond = zeros(arr_size, 1);

% QR Factorization

QR_times = zeros(trials,1); % Array of timing data for the QR factorization trials
QR_mean_time = zeros(arr_size, 1);

QR_errors = zeros(trials, 1);
QR_mean_error = zeros(arr_size, 1);

QR_condition = zeros(trials, 1);
QR_mean_cond = zeros(arr_size, 1);

% GE / LU factorization

GE_times = zeros(trials,1); % Array of timing data for the GE factorization trials
GE_mean_time = zeros(arr_size, 1);

GE_errors = zeros(trials, 1);
GE_mean_error = zeros(arr_size, 1);

GE_condition = zeros(trials, 1);
GE_mean_cond = zeros(arr_size, 1);

for j=2:arr_size % Skip N = 0
    
    j % Progress report
    
    N = N_range(j)
    
    for i=1:trials
        
        % Form a random matrix A and right-hand side b (normally distributed)
        A = randn(N,N);
        x = ones(N,1);
        
        b = A*x;
        
        A_condition(i) = cond(A);
        
        % Apply backslash and calculate time taken for QR factorization
        tic;
        [Q,R] = qr(A); % Compute the QR factorization of A        
        xhatQR = R\(Q.'*b);  % Solve the linear system (note that Q' is the transpose of Q)
        QR_times(i) = toc;
        
        QR_condition(i) = cond(Q) * cond(R);
        
        error = abs(xhatQR - x);
        QR_errors(i) = max(error);
        
        % Apply backslash and calculate time taken for GE factorization
        tic;
        xhatGE = A\b; % Solve the linear system using GE with partial pivoting
        GE_times(i) = toc;
        
        [L, U] = lu(A); % For calculating the condition of the LU matrix
        GE_condition(i) = cond(L) * cond(U);
        
        error = abs(xhatGE - x);
        GE_errors(i) = max(error);
        
    end
    
    A_mean_cond(j) = mean(A_condition);
    
    QR_mean_time(j) = mean(QR_times);
    QR_mean_error(j) = mean(QR_errors);
    QR_mean_cond(j) = mean(QR_condition);
    
    GE_mean_time(j) = mean(GE_times);
    GE_mean_error(j) = mean(GE_errors);
    GE_mean_cond(j) = mean(GE_condition);
    
end

semilogy(N_range, QR_mean_time, N_range, GE_mean_time);
title('Fig. 1a - Mean QR vs GE times');
legend({'QR', 'GE'},'fontsize',16,'Location','northwest');
xlabel('N','fontsize',16);
ylabel('Time (seconds)','fontsize',16);

figure;semilogy(N_range, QR_mean_error, '*', N_range, GE_mean_error, '*');
title('Fig. 1b - QR vs GE Actual - Mean Absolute Errors');
legend({'QR', 'GE'},'fontsize',16,'Location','northwest');
xlabel('N','fontsize',16);
ylabel('Error','fontsize',16);

figure;semilogy(N_range, QR_mean_error, N_range, GE_mean_error);
title('Fig. 1b - QR vs GE Actual - Mean Absolute Errors');
legend({'QR', 'GE'},'fontsize',16,'Location','northwest');
xlabel('N','fontsize',16);
ylabel('Error','fontsize',16);


figure;semilogy(N_range, QR_mean_cond, N_range, GE_mean_cond, N_range, A_mean_cond);
title('Fig. 1c - Condition of A vs QR vs LU');
legend({'QR', 'LU', 'A'},'fontsize',16,'Location','northwest');
xlabel('N','fontsize',16);
ylabel('Condition','fontsize',16);
