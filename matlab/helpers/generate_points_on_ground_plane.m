function [R1g, R2g, f, H, x1, x2, R, t] = generate_points_on_ground_plane(N, m)
if nargin < 2
    m = 0;
end
if nargin < 1
    N = 100;
end

% Assure that the ground truth can be trusted by assuring that the
% reprojection error is not too large.
reproj = inf;

while norm(reproj, 'fro') >= 1e-14

    % Generate points on GT plane - y-axis is aligned with plane normal
    X = [randn(1,N); zeros(1,N); randn(1,N); ones(1,N)];

    % Generate rotation matrices (known overhead - IMU data)
    R1g = rotm('x',randn) * rotm('z',randn);
    R2g = rotm('x',randn) * rotm('z',randn);

    % Generate partially calibrated cameras
    f = 1 + 9 * rand;
    K = diag([f,f,1]);
    Kinv = diag([1 / f, 1 / f, 1]);

    Ry1 = rotm('y',randn);
    Ry2 = rotm('y',randn);

    R1 = R1g * Ry1;
    R2 = R2g * Ry2;

    % Relative rotation
    R = R2 * R1';

    t1 = randn(3,1);
    t2 = randn(3,1);

    % Relative translation
    t = -R * t1 + t2;
    %R = R1' * R2;

    P1 = K * [R1 t1];
    P2 = K * [R2 t2];

    % Ground truth homography
    H = K * [R2(:, [1 3]) t2] / [R1(:, [1 3]) t1] * Kinv;
    H = H / H(end);

    % Generate point correspondences (pinhole)
    x1 = pflat(P1 * X);
    x2 = pflat(P2 * X);

    % Check GT
    reproj = pflat(H * x1) - x2;
        
    % Add non-planar points
    if m > 0
        X = [X [randn(3,m); ones(1,m)]];
        x1 = pflat(P1 * X);
        x2 = pflat(P2 * X);
    end
end