[v, f] = read_vtk('lh.white.vtk');
H = load('lh.white.H.txt');
write_property('lh.white.H.vtk', v, f, struct('H', H));
[v_inf, f_inf] = read_vtk('lh.inflated.vtk');
write_property('lh.inflated.H.vtk', v_inf, f_inf, struct('H', H));
 
connections = [];

for i = 1:size(f, 1)
    b = nchoosek(f(i, :), 2);
    connections = [connections; b; b(:, [2, 1])];
end

connections = unique(connections, 'rows');
connections = connections + 1;
edge_dists = sqrt(sum(power(connections(:, 1) - connections(:, 2), 2), 2));
avg = sum(edge_dists) / size(edge_dists, 1);
gaus = gaussian(edge_dists, 0, avg);
v_size = size(v, 1);

sp = sparse(connections(:, 1), connections(:, 2), gaus, v_size, v_size);
sum_rows = sum(sp, 2);
sum_rows(~sum_rows) = 1;
sp = sp ./ sum_rows;

for i=1:40
    H = sp * H;
    if i == 10
        write_property('lh.inflated.H.10.vtk', v_inf, f_inf, struct('H', H));
    elseif i == 20
        write_property('lh.inflated.H.20.vtk', v_inf, f_inf, struct('H', H));
    elseif i == 40
        write_property('lh.inflated.H.40.vtk', v_inf, f_inf, struct('H', H));
    end
end

for i=1:40
    v = sp * v;
    if i == 10
        write_vtk('lh.white.H.10.vtk', v, f);
    elseif i == 20
        write_vtk('lh.white.H.20.vtk', v, f);
    elseif i == 40
        write_vtk('lh.white.H.40.vtk', v, f);
    end
end

