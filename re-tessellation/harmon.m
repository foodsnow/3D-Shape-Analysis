icos = [4, 5 ,6];
for ico=icos
    gen_harmon(ico);
end

function gen_harmon(ico)
    [v_sph, f_sph] = read_vtk("lh.sphere.vtk");
    [v, f] = read_vtk("lh.white.vtk");
    [v_ico, f_ico] = read_vtk(strcat("icosphere_", num2str(ico),".vtk"));
    H = load("lh.white.H.txt");
    
    L = 40;
    Y = zeros(size(v_sph, 1), (L+1)^2);
    for l=0:40
        Y(:, l^2 + 1: (l + 1)^2) = spharm_real(v_sph, l);
    end
    
    orig_system_cords = Y \ v;
    orig_system_features = Y \ H;
    
    Y_ico = zeros(size(v_ico, 1), (L+1)^2);
    for l=0:L
        Y_ico(:, l^2 + 1: (l + 1)^2) = spharm_real(v_ico, l);
    end
    
    v_ico = Y_ico * orig_system_cords;
    H_ico = Y_ico * orig_system_features;
    
    write_property(strcat("lh.white.har.ico", num2str(ico), ".H.vtk"), v_ico, f_ico, struct('H', H_ico));
end