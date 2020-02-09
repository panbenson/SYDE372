% lab 2
% Benson Pan
% 3 Classifiers

% Case 1
% class A
N_A = 200;
cov_A = [8 0; 0 4];
m_A = [5 10]';
% class B
N_B = 200;
cov_B = [8 0; 0 4];
m_B = [10 15]';

% Case 2
% class C
N_C = 100;
cov_C = [8 4; 4 40];
m_C = [5 10]';
% class D
N_D = 200;
cov_D = [8 0; 0 8];
m_D = [15 10]';
% class E
N_E = 150;
cov_E = [10 -5; -5 20];
m_E = [10 5]';

% generate random points for the distributions
A = generate_clusters(N_A, m_A, cov_A);
B = generate_clusters(N_B, m_B, cov_B);
C = generate_clusters(N_C, m_C, cov_C);
D = generate_clusters(N_D, m_D, cov_D);
E = generate_clusters(N_E, m_E, cov_E);

% clean up old figures
close all
clf('reset');

% following the idea of creating a 2D plot using 0.1 and the 
% min and max of each dimension per plot
case_1_dim_max = max([max(A, [], 2) max(B, [], 2)], [], 2);
case_1_dim_min = min([min(A, [], 2) min(B, [], 2)], [], 2);
case_1_x1_range = floor(case_1_dim_min(1)):0.1:ceil(case_1_dim_max(1));
case_1_x2_range = floor(case_1_dim_min(2)):0.1:ceil(case_1_dim_max(2));

% the 2d plot according to the size of the data
case_1_plot_MED = zeros(size(case_1_x1_range, 2), size(case_1_x2_range, 2));
case_1_plot_MICD = zeros(size(case_1_x1_range, 2), size(case_1_x2_range, 2));
case_1_plot_MAP = zeros(size(case_1_x1_range, 2), size(case_1_x2_range, 2));
case_1_plot_NN = zeros(size(case_1_x1_range, 2), size(case_1_x2_range, 2));
case_1_plot_KNN = zeros(size(case_1_x1_range, 2), size(case_1_x2_range, 2));

for i = 1:size(case_1_x1_range, 2)
    for k = 1:size(case_1_x2_range, 2)
        x = [case_1_x1_range(i) case_1_x2_range(k)];
        A_dist = sqrt(sum((x' - m_A).^2));
        B_dist = sqrt(sum((x' - m_B).^2));
        A_MICD_dist = (x' - m_A)' * inv(cov_A) * (x' - m_A);
        B_MICD_dist = (x' - m_B)' * inv(cov_A) * (x' - m_B);
        P_A = N_A / (N_A + N_B);
        P_B = N_B / (N_A + N_B);
        MAP_LHS = exp((-1/2) * A_MICD_dist) / exp((-1/2) * B_MICD_dist);
        MAP_RHS = (det(cov_A)^(1/2) * P_B) / (det(cov_B)^(1/2) * P_A);
        NN_A = k_nearest(A, x, 1);
        NN_B = k_nearest(B, x, 1);
        KNN_A = k_nearest(A, x, 5);
        KNN_B = k_nearest(B, x, 5);
        
        % MED
        if A_dist < B_dist
            case_1_plot_MED(i, k) = 0;
        else
            case_1_plot_MED(i, k) = 1;
        end
        
        % MICD
        if A_MICD_dist < B_MICD_dist
            case_1_plot_MICD(i, k) = 0;
        else
            case_1_plot_MICD(i, k) = 1;
        end
        
        % MAP 
        if MAP_LHS > MAP_RHS
            case_1_plot_MAP(i, k) = 0;
        else
            case_1_plot_MAP(i, k) = 1;
        end
        
        % NN
        if NN_A < NN_B
            case_1_plot_NN(i, k) = 0;
        else
            case_1_plot_NN(i, k) = 1;
        end
        
        % KNN
        if KNN_A < KNN_B
            case_1_plot_KNN(i, k) = 0;
        else
            case_1_plot_KNN(i, k) = 1;
        end
    end
end

case_2_dim_max = max([max(C, [], 2) max(D, [], 2) max(E, [], 2)], [], 2);
case_2_dim_min = min([min(C, [], 2) min(D, [], 2) max(E, [], 2)], [], 2);
case_2_x1_range = floor(case_2_dim_min(1)):0.1:ceil(case_2_dim_max(1));
case_2_x2_range = floor(case_2_dim_min(2)):0.1:ceil(case_2_dim_max(2));

% the 2d plot according to the size of the data
case_2_plot_MED = zeros(size(case_2_x1_range, 2), size(case_2_x2_range, 2));
case_2_plot_MICD = zeros(size(case_2_x1_range, 2), size(case_2_x2_range, 2));
case_2_plot_MAP = zeros(size(case_2_x1_range, 2), size(case_2_x2_range, 2));
case_2_plot_NN = zeros(size(case_2_x1_range, 2), size(case_2_x2_range, 2));
case_2_plot_KNN = zeros(size(case_2_x1_range, 2), size(case_2_x2_range, 2));

for i = 1:size(case_2_x1_range, 2)
    for k = 1:size(case_2_x2_range, 2)
        x = [case_2_x1_range(i) case_2_x2_range(k)];
        C_dist = sqrt(sum((x' - m_C).^2));
        D_dist = sqrt(sum((x' - m_D).^2));
        E_dist = sqrt(sum((x' - m_E).^2));
        C_MICD_dist = (x' - m_C)' * inv(cov_C) * (x' - m_C);
        D_MICD_dist = (x' - m_D)' * inv(cov_D) * (x' - m_D);
        E_MICD_dist = (x' - m_E)' * inv(cov_E) * (x' - m_E);
        P_C = N_C / (N_C + N_D + N_E);
        P_D = N_D / (N_C + N_D + N_E);
        P_E = N_E / (N_C + N_D + N_E);
        MAP_CD_LHS = exp((-1/2) * C_MICD_dist) / exp((-1/2) * D_MICD_dist);
        MAP_CD_RHS = (det(cov_C)^(1/2) * P_D) / (det(cov_D)^(1/2) * P_C);
        MAP_CE_LHS = exp((-1/2) * C_MICD_dist) / exp((-1/2) * E_MICD_dist);
        MAP_CE_RHS = (det(cov_C)^(1/2) * P_E) / (det(cov_E)^(1/2) * P_C);
        MAP_DE_LHS = exp((-1/2) * D_MICD_dist) / exp((-1/2) * E_MICD_dist);
        MAP_DE_RHS = (det(cov_D)^(1/2) * P_E) / (det(cov_E)^(1/2) * P_D);
        NN_C = k_nearest(C, x, 1);
        NN_D = k_nearest(D, x, 1);
        NN_E = k_nearest(E, x, 1);
        KNN_C = k_nearest(C, x, 5);
        KNN_D = k_nearest(D, x, 5);
        KNN_E = k_nearest(E, x, 5);
        
        % MED
        if C_dist < D_dist && C_dist < E_dist
            case_2_plot_MED(i, k) = 0;
        elseif D_dist < C_dist && D_dist < E_dist
            case_2_plot_MED(i, k) = 1;
        else
            case_2_plot_MED(i, k) = 2;
        end
        
        % MICD
        if C_MICD_dist < D_MICD_dist && C_MICD_dist < E_MICD_dist
            case_2_plot_MICD(i, k) = 0;
        elseif D_MICD_dist < C_MICD_dist && D_MICD_dist < E_MICD_dist
            case_2_plot_MICD(i, k) = 1;
        else
            case_2_plot_MICD(i, k) = 2;
        end
        
        % MAP
        if MAP_CD_LHS > MAP_CD_RHS
            if MAP_CE_LHS > MAP_CE_RHS
                case_2_plot_MAP(i, k) = 0;
            else
                case_2_plot_MAP(i, k) = 2;
            end
        else
            if MAP_DE_LHS > MAP_DE_RHS
                case_2_plot_MAP(i, k) = 1;
            else
                case_2_plot_MAP(i, k) = 2;
            end
        end
        
        % NN
        if NN_C < NN_D && NN_C < NN_E
            case_2_plot_NN(i, k) = 0;
        elseif NN_D < NN_C && NN_D < NN_E
            case_2_plot_NN(i, k) = 1;
        else
            case_2_plot_NN(i, k) = 2;
        end
        
        %KNN
        if KNN_C < KNN_D && KNN_C < KNN_E
            case_2_plot_KNN(i, k) = 0;
        elseif KNN_D < KNN_C && KNN_D < KNN_E
            case_2_plot_KNN(i, k) = 1;
        else
            case_2_plot_KNN(i, k) = 2;
        end
    end
end

% plot case 1 (MED, MICD, MAP)
figure(1);
hold on;
scatter( A(1,:), A(2,:) );
plot_ellipse(m_A(1), m_A(2), 0, sqrt(cov_A(1, 1)), sqrt(cov_A(2, 2)));
scatter( B(1,:), B(2,:) );
plot_ellipse(m_B(1), m_B(2), 0, sqrt(cov_B(1, 1)), sqrt(cov_B(2, 2)));
[X, Y] = meshgrid(case_1_x1_range, case_1_x2_range);
contour(X, Y, case_1_plot_MED', 'k');
contour(X, Y, case_1_plot_MICD', 'r');
contour(X, Y, case_1_plot_MAP', 'g');
legend('class A', 'class A s.d.', 'class B', 'class B s.d.', 'MED', 'MICD', 'MAP');
hold off;

% plot case 1 (NN, KNN)
figure(2);
hold on;
scatter( A(1,:), A(2,:) );
plot_ellipse(m_A(1), m_A(2), 0, sqrt(cov_A(1, 1)), sqrt(cov_A(2, 2)));
scatter( B(1,:), B(2,:) );
plot_ellipse(m_B(1), m_B(2), 0, sqrt(cov_B(1, 1)), sqrt(cov_B(2, 2)));
[X, Y] = meshgrid(case_1_x1_range, case_1_x2_range);
contour(X, Y, case_1_plot_NN', 'k');
contour(X, Y, case_1_plot_KNN', 'r');
legend('class A', 'class A s.d.', 'class B', 'class B s.d.', 'NN', 'KNN');
hold off;

% plot case 2 (MED, MICD, MAP)
[V_C, D_C] = eig(cov_C);
[V_E, D_E] = eig(cov_E);
figure(3);
hold on;
scatter( C(1,:), C(2,:) );
plot_ellipse(m_C(1), m_C(2), atan(V_C(2, 1)/V_C(1, 1)), sqrt(cov_C(1, 1)), sqrt(cov_C(2, 2)));
scatter( D(1,:), D(2,:) );
plot_ellipse(m_D(1), m_D(2), 0, sqrt(cov_D(1, 1)), sqrt(cov_D(2, 2)));
scatter( E(1,:), E(2,:) );
plot_ellipse(m_E(1), m_E(2), atan(V_E(2, 1)/V_E(1, 1)), sqrt(cov_E(1, 1)), sqrt(cov_E(2, 2)));
[X, Y] = meshgrid(case_2_x1_range, case_2_x2_range);
contour(X, Y, case_2_plot_MED', 'k');
contour(X, Y, case_2_plot_MICD', 'r');
contour(X, Y, case_2_plot_MAP', 'g');
legend('class C', 'class C s.d.', 'class D', 'class D s.d.', 'class E', 'class E s.d.', 'MED', 'MICD', 'MAP');
hold off;

% plot case 2 (MED, MICD, MAP)
figure(4);
hold on;
scatter( C(1,:), C(2,:) );
plot_ellipse(m_C(1), m_C(2), atan(V_C(2, 1)/V_C(1, 1)), sqrt(cov_C(1, 1)), sqrt(cov_C(2, 2)));
scatter( D(1,:), D(2,:) );
plot_ellipse(m_D(1), m_D(2), 0, sqrt(cov_D(1, 1)), sqrt(cov_D(2, 2)));
scatter( E(1,:), E(2,:) );
plot_ellipse(m_E(1), m_E(2), atan(V_E(2, 1)/V_E(1, 1)), sqrt(cov_E(1, 1)), sqrt(cov_E(2, 2)));
[X, Y] = meshgrid(case_2_x1_range, case_2_x2_range);
contour(X, Y, case_2_plot_NN', 'k');
contour(X, Y, case_2_plot_KNN', 'r');
legend('class C', 'class C s.d.', 'class D', 'class D s.d.', 'class E', 'class E s.d.', 'NN', 'KNN');
hold off;