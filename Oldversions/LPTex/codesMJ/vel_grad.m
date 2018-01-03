function det_p = vel_grad(euler_angle)
dummy6 =1;
for p=-2.5:0.01:2.5
    
    strain_external = [1 0 p;0 0 0;p 0 -1];
    g_mat = Eulertogmat(euler_angle);
    s_in= g_mat*strain_external*g_mat';
    for counter1=1:1:56 % Calculating m using bishop-hill table and picking max m as taylor factor 
            m(counter1) = round(sqrt(6)*( (bishop(counter1,1)*s_in(2,2)) - ((bishop(counter1,2)*s_in(1,1))) + ((bishop(counter1,4)*2*s_in(2,3))) + ((bishop(counter1,5)*2*s_in(3,1))) +((bishop(counter1,6)*2*s_in(1,2)))),8);
    end
    Taylor_factor(dummy6) = max(m);
    dummy6 = dummy6+1; 
end
[min_taylor, index_min_taylor] = min(Taylor_factor);
det_p=(-2.5) + 0.01*(index_min_taylor-1);
end