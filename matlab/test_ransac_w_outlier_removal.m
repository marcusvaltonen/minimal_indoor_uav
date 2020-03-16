addpath('helpers')
addpath('solvers')

% Standard deviation of noise
noise = 1e-4;

% Robust estimation using RANSAC
threshold = 1000 * noise;
nbr_iter = 200;

% Generate NN images on a the ground plane and m other non-planar
NN = 100;
m = 30;
N = NN + m;
[R1g, R2g, f_gt, H_gt, x1, x2] = generate_points_on_ground_plane(NN, m);

% Add noise (inhomogenous)
x1n = x1(1:2,:) + noise * randn(2,N);
x2n = x2(1:2,:) + noise * randn(2,N);

% Add outliers
p = 0.2;
idx_outliers = repmat(rand(1, N) < p, 2, 1);
x1n(idx_outliers) = 10 * randn(1,nnz(idx_outliers));
x2n(idx_outliers) = 10 * randn(1,nnz(idx_outliers));

% Make homogenous again
x1n = [x1n; ones(1,N)];
x2n = [x2n; ones(1,N)];

% Test our method
tic;
[H, f, inliners, history] = ransac_floor_fHf(x1n, x2n, R1g, R2g, ...
                                                   nbr_iter, threshold);
time_fHf = toc;

% Normalize
H_gt = H_gt / H_gt(end);
H = H / H(end);

% Display some information about the problem instance
fprintf(1, 'Noise standard deviation: %g\n', noise)
fprintf(1, 'Number of outliers (planar): %d\n', m+nnz(idx_outliers(1,1:NN)))
fprintf(1, 'Number of outliers    (all): %d\n\n', nnz(idx_outliers(1,:)))
fprintf(1, 'Estimated number of outliers %d\n', N - sum(inliners))
fprintf(1, 'Homography error: %g\n', norm(H - H_gt, 'fro') / norm(H_gt, 'fro'))
fprintf(1, 'Focal length error: %g\n', abs(f - f_gt) / f_gt)