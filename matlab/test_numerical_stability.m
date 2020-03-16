addpath('helpers')
addpath('solvers')

nbr_iter = 1e4;

error_H1 = zeros(nbr_iter,1);
error_f1 = zeros(nbr_iter,1);
error_H2 = zeros(nbr_iter,1);
error_f2 = zeros(nbr_iter,1);
error_H3 = zeros(nbr_iter,1);
error_f3 = zeros(nbr_iter,1);

N = 3;

for j = 1:nbr_iter
    % Generate points on ground plane (this generates f1 = f2, hence
    % compatible with all problem formulations)
    [R1g, R2g, f_gt, H_gt, p1, p2] = generate_points_on_ground_plane(N);

    % f1Hf2 (3 pt)
    data = get_floor_f1Hf2(p1, p2, R1g, R2g);
    
    % fHf (2.5 pt)
    [H2, f2] = get_floor_fHf(p1(1:2, :), p2(1:2, :), R1g, R2g);

    % Hf (2.5 pt)
    [H3, f3] = get_floor_Hf(p1(1:2, :), p2(1:2, :), R1g, R2g, f_gt);
    
    % Normalize
    H2 = H2 / H2(end);
    H3 = H3 / H3(end);
    H_gt = H_gt / H_gt(end);
    
    % Compute putative errors (f1Hf2)
    nbr_sols = length(data);
    tmp_err_f1 = zeros(1, nbr_sols); 
    tmp_err_f2 = zeros(1, nbr_sols); 
    tmp_err_H = zeros(1, nbr_sols); 
    for k = 1:nbr_sols
        tmp_err_f1(k) = abs(data(k).f1 -f_gt) / f_gt;
        tmp_err_f2(k) = abs(data(k).f2 -f_gt) / f_gt;
        Htmp = data(k).H;
        Htmp = Htmp / Htmp(end);
        tmp_err_H(k) = norm(Htmp - H_gt, 'fro') / norm(H_gt, 'fro');
    end 

    % Compute errors
    error_f1(j) = mean([min(tmp_err_f1), min(tmp_err_f2)]);
    error_H1(j) = min(tmp_err_H);
    error_H2(j) = norm(H2 - H_gt, 'fro') / norm(H_gt, 'fro');
    error_f2(j) = abs(f2 -f_gt) / f_gt;
    error_H3(j) = norm(H3 - H_gt, 'fro') / norm(H_gt, 'fro');
    error_f3(j) = abs(f3 -f_gt) / f_gt;
end

%% Histograms
M = -16:0.25:2;

figure(1)
[~,edges] = histcounts(log10(error_H1), M);
histogram(error_H1,10.^edges);
hold on
[~,edges] = histcounts(log10(error_H2), M);
histogram(error_H2,10.^edges);
[~,edges] = histcounts(log10(error_H3), M);
histogram(error_H3,10.^edges);
hold off
set(gca, 'xscale', 'log')
title('Homography error', 'interpreter', 'latex')
h = legend('$f_1Hf_2$','$fHf$','$Hf$');
set(h, 'fontsize', 14);
set(h, 'interpreter', 'latex');
set(gca,'TickLabelInterpreter', 'latex', 'fontsize', 14);

figure(2)
[~,edges] = histcounts(log10(error_f1), M);
histogram(error_f1,10.^edges);
hold on
[~,edges] = histcounts(log10(error_f2), M);
histogram(error_f2,10.^edges);
[~,edges] = histcounts(log10(error_f3), M);
histogram(error_f3,10.^edges);
hold off
set(gca, 'xscale', 'log')
title('Focal length error', 'interpreter', 'latex')
h = legend('$f_1Hf_2$','$fHf$','$Hf$');
set(h, 'fontsize', 14);
set(h, 'interpreter', 'latex');
set(gca,'TickLabelInterpreter', 'latex', 'fontsize', 14);