% lab 2
% Benson Pan
% 2 Generating Clusters

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
close all;
clf('reset');

% plot case 1
figure(1);
hold on;
scatter( A(1,:), A(2,:) );
plot_ellipse(m_A(1), m_A(2), 0, sqrt(cov_A(1, 1)), sqrt(cov_A(2, 2)));
scatter( B(1,:), B(2,:) );
plot_ellipse(m_B(1), m_B(2), 0, sqrt(cov_B(1, 1)), sqrt(cov_B(2, 2)));
legend('class A', 'class A s.d.', 'class B', 'class B s.d.');
hold off;

% plot case 2
[V_C, D_C] = eig(cov_C);
[V_E, D_E] = eig(cov_E);
figure(2);
hold on;
scatter( C(1,:), C(2,:) );
plot_ellipse(m_C(1), m_C(2), atan(V_C(2, 1)/V_C(1, 1)), sqrt(cov_C(1, 1)), sqrt(cov_C(2, 2)));
scatter( D(1,:), D(2,:) );
plot_ellipse(m_D(1), m_D(2), 0, sqrt(cov_D(1, 1)), sqrt(cov_D(2, 2)));
scatter( E(1,:), E(2,:) );
plot_ellipse(m_E(1), m_E(2), atan(V_E(2, 1)/V_E(1, 1)), sqrt(cov_E(1, 1)), sqrt(cov_E(2, 2)));
legend('class C', 'class C s.d.', 'class D', 'class D s.d.', 'class E', 'class E s.d.');
hold off;