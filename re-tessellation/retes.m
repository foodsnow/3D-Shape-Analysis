

for iso = [5]
    gen_iso(iso);
end

function gen_iso(iso_num)
    [v, f] = read_vtk("lh.sphere.vtk");
    f = f + 1;
    [v_iso, f_iso] = read_vtk(strcat('icosphere_', num2str(iso_num), '.vtk'));
    f_iso = f_iso + 1;

    % v = v_iso;
    % f = f_iso;

    MD = KDTreeSearcher(v);
    tr = triangulation(f, v);

    baries = zeros(size(v_iso, 1), 3);
    closest_triangles = zeros(size(v_iso, 1), 3);
    
    k = 1;
    Q = num2cell(v_iso, 2);
    numQ = length(Q);
    toRemove = false(1, numQ);
    while ~all(toRemove)
        for i = 1:numQ
            if toRemove(i)
                continue
            end
            q = Q{i};
            p = MD.knnsearch(q, 'k', k, 'distance', 'euclidean');
            faces = tr.vertexAttachments(p(k));
            for face_id = faces{:}
                face = f(face_id, :);
                q_proj = proj(q, face, v);
                br = tr.cartesianToBarycentric(face_id, q_proj);
                if all(br >= 0)
                    closest_triangles(i, :) = face;
                    baries(i, :) = br;
                    toRemove(i) = true;
                end

            end
        end
        disp(k);
        disp(sum(toRemove));
        k = k + 1;
    end


    [v_orig, f_orig] = read_vtk('lh.white.vtk');
    H = load('lh.white.H.txt');
    new_H = zeros(size(H));
    new_v_iso = zeros(size(v_iso, 1), 3);
    for i=1:size(v_iso, 1)
        A = v_orig(closest_triangles(i, 1), :);
        B = v_orig(closest_triangles(i, 2), :);
        C = v_orig(closest_triangles(i, 3), :);
        a = baries(i, 1);
        b = baries(i, 2);
        c = baries(i, 3);

        new_v_iso(i, :) = A * a + B * b + C * c;
        new_H(i) = H(closest_triangles(i, 1)) * a + H(closest_triangles(i, 2)) * b + H(closest_triangles(i, 3)) * c;
    end

    write_property(strcat('lh.white.iso', num2str(iso_num), '.H.vtk'), new_v_iso, f_iso-1, struct('H', new_H));
end

function Q_prime = proj(q, face, v)
    cords = v(face, :);
    p1 = cords(1, :);
    p2 = cords(2, :);
    p3 = cords(3, :);
    N = cross(p2 - p1, p3 - p1);
    
    A = p2;
    Q = q;
    proj_A = dot(A, N) * N / dot(N, N);
    proj_Q = dot(Q, N) * N / dot(N, N);

    dist_A_projA = norm(proj_A - A);
    dist_Q_projQ = norm(proj_Q - Q);
    scale_factor = dist_A_projA / dist_Q_projQ;
    Q_prime = A + scale_factor * (Q - proj_Q);

    % plot_line(cords(1, :), cords(2, :), cords(3, :));
    % hold on;
    % scatter3(cords(:, 1), cords(:, 2), cords(:, 3));
    % scatter3(Q(:, 1), Q(:, 2), Q(:, 3));
    % scatter3(Q_prime(:, 1), Q_prime(:, 2), Q_prime(:, 3));

    % Projection using orthogonal projection
    N = N / norm(N);
    d = dot(A, N);
    dist_Q = dot(Q, N);
    % Q_prime = Q - (dist_Q -d) * N;
end

function [normal, d] = plot_line(p1, p2, p3)
    normal = cross(p1 - p2, p1 - p3);
    d = p1(1)*normal(1) + p1(2)*normal(2) + p1(3)*normal(3);
    d = -d;
    x = -5:5; y = -5:5;
    [X,Y] = meshgrid(x,y);
    Z = (-d - (normal(1)*X) - (normal(2)*Y))/normal(3);
    mesh(X,Y,Z)
end