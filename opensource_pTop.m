load example.mat
addpath('pTop')

r0 = r;
b0 = b;

len = size(nvr, 1);
coef = double(0);
nvr = double(nvr);
for i = 1 : len
   rr = double(r0(i, :).');
   bb = double(b0(i, :).');
   coef = coef + double(1) / len * coef_J_pure_pTop_sQPEP_merge_s(rr, bb, nvr(i, :).');
end
coef_J = coef;

coef_tq1 = coeftq1_pTop_sQPEP_merge_s_sQPEP(coef);
coef_tq2 = coeftq2_pTop_sQPEP_merge_s_sQPEP(coef);
coef_tq3 = coeftq3_pTop_sQPEP_merge_s_sQPEP(coef);
coefs_tq=[coef_tq1;
    coef_tq2;
    coef_tq3];
G = G_pTop_sQPEP_merge_s_sQPEP(coef);
pinvG = pinv(G);

coef_Jacob_qt1 = coef_Jacob_qt1_pTop_sQPEP_merge_s_sQPEP(coefs_tq, pinvG, coef_J);
coef_Jacob_qt2 = coef_Jacob_qt2_pTop_sQPEP_merge_s_sQPEP(coefs_tq, pinvG, coef_J);
coef_Jacob_qt3 = coef_Jacob_qt3_pTop_sQPEP_merge_s_sQPEP(coefs_tq, pinvG, coef_J);
coef_Jacob_qt4 = coef_Jacob_qt4_pTop_sQPEP_merge_s_sQPEP(coefs_tq, pinvG, coef_J);
coef_Jacob_qt = [coef_Jacob_qt1;
    coef_Jacob_qt2;
    coef_Jacob_qt3;
    coef_Jacob_qt4];

QQ = Q_sym_pTop_sQPEP_merge_s_sQPEP(coefs_tq, pinvG, coef);
WW = W_pTop_sQPEP_merge_s_sQPEP(coefs_tq, pinvG, coef);


[Q]=sQPEP_L2_Solver(WW,QQ, 123456);
cnt=0;
real_res_q=zeros(0,4);
for i = 1 : length(Q)
    tmp = Q{i};
    if (~isreal(tmp(1)))
        continue;
    end
    if (isnan(tmp(1)) || sum(Q{i})==0 )
        continue;
    end
    cnt = cnt + 1;
    real_s(cnt) = Q{i}.'*Q{i};
    real_res_q(cnt,:) = tmp/sqrt(real_s(cnt));
end


ret_q = real_res_q;
ret_s = real_s;


best_idx = 1;
best_cost = inf;
for i = 1 : length(ret_s)
    cur_s = ret_s(i);
    cur_q = ret_q(i,:);
    t1 = t1_pTop_sQPEP_merge_s_sQPEP(pinvG, coefs_tq, cur_q.'*sqrt(cur_s));
    t2 = t2_pTop_sQPEP_merge_s_sQPEP(pinvG, coefs_tq, cur_q.'*sqrt(cur_s));
    t3 = t3_pTop_sQPEP_merge_s_sQPEP(pinvG, coefs_tq, cur_q.'*sqrt(cur_s));
    cur_cost = J_func_pTop_fast((cur_q).'*(sqrt(cur_s)), ([t1;t2;t3]), (r0), (b0), (nvr));
    if (best_cost > cur_cost)
        best_cost = cur_cost;
        best_idx = i;
    end
    %best_cost
end

best_q = ret_q(best_idx,:);
best_s = ret_s(best_idx);

best_t1 = t1_pTop_sQPEP_merge_s_sQPEP(pinvG, coefs_tq, best_q.'*sqrt(best_s));
best_t2 = t2_pTop_sQPEP_merge_s_sQPEP(pinvG, coefs_tq, best_q.'*sqrt(best_s));
best_t3 = t3_pTop_sQPEP_merge_s_sQPEP(pinvG, coefs_tq, best_q.'*sqrt(best_s));
best_t = [best_t1; best_t2; best_t3];

best_s

function J = J_func_pTop_fast(q, t, r, b, nv)
    R = q2R(q); 
    t_matrix = repmat(t', size(r, 1), 1); 
    transformed_r = r * R' + t_matrix; 
    differences = transformed_r - b; 
    projections = sum(differences .* nv, 2); 
    J = sum(projections.^2); 
end
function J = J_func_pTop(q, t, r, b, nv)
len = size(r, 1);
R = q2R(q);
J = 0;
for i = 1 : len
   rr = r(i, :).';
   bb = b(i, :).';
   J = J + 1 / len * (nv(i, :) * (R * rr + t - bb))^2; 
end
end
function R = q2R(q)
q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
R = [
        q0^2 + q1^2 - q2^2 - q3^2,         2*q0*q3 + 2*q1*q2,         2*q1*q3 - 2*q0*q2;
                2*q1*q2 - 2*q0*q3, q0^2 - q1^2 + q2^2 - q3^2,         2*q0*q1 + 2*q2*q3;
                2*q0*q2 + 2*q1*q3,         2*q2*q3 - 2*q0*q1, q0^2 - q1^2 - q2^2 + q3^2];
end
