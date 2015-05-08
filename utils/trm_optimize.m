function optimized_trm = trm_optimize(trm, n_angles, iter_num, superpose)

n_angles = min(n_angles, size(trm.psi,1));

angle_indices = trmdistantangleindices(trm, n_angles);
f = @(x) trmobjfunc(trm, angle_indices, x, superpose);

initial_point = trm.psi(angle_indices,2:end-1);
initial_point = initial_point(:);
initial_point = reduceangles(initial_point);

options = optimoptions('fmincon');
options = optimoptions(options,'Display', 'iter');
options = optimoptions(options,'MaxIter', iter_num);
options = optimoptions(options,'MaxFunEvals', Inf);
options = optimoptions(options,'GradObj','on');

n = size(initial_point);
x = fmincon(f,initial_point,[],[],[],[],-pi*ones(n),pi*ones(n),[],...
    options);

optimized_trm = trm;
optimized_trm.psi(angle_indices, 2:end-1) = reshape(x, ...
    length(angle_indices), size(optimized_trm.psi, 2) - 2);
if superpose
    optimized_trm.U = trmsuperpos(optimized_trm);
end

end

