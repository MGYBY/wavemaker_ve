clear all; clc;

fr = 0.84;
swe_sol = load('out-5.4759.txt');
[sol_row, sol_col] = size(swe_sol);
reconst_vert_reso = 36;
reconst_points = reconst_vert_reso*sol_row;
reconst_v = zeros(reconst_points,3);

for i = 1:1:sol_row
    h = swe_sol(i, 2);
    u_ave = swe_sol(i, 3);
    x_coord = swe_sol(i, 1)*fr^2;
    for j = 1:1:reconst_vert_reso
        row_ind = (i-1)*reconst_vert_reso+j;
        reconst_v(row_ind, 1) = x_coord;
        reconst_v(row_ind, 2) = h/reconst_vert_reso*j;
        reconst_v(row_ind, 3) = 1.50*u_ave*(1.0-(1.0-reconst_v(row_ind, 2)/h)^2);
    end
end

outname = 'reconstructed_SWE_vel_Fr=0.84';
% fid = fopen(outname,'w');

csvwrite (outname, reconst_v)