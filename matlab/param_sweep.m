param_list = {'param1','param2'};
p = parameters();
nominals = [1,1];
ratios = [.8, 1, 1.2];

for i = 1:length(param_list)
    for j = 1:length(ratios)
        p.(param_list{i}) = nominals(i) * ratios(i);
        x_star(i,j,:) = run_optimization(p);
    end
    dx_star = (x_star(i,end,:) - x_star(i,1,:)) ./ x_star(i,2,:);
    dx_star_dp(i,:) = dx_star / (ratios(end) - ratios(1));
end


