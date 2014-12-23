function run_matlab_sims_mod( N, data_dir, myseed )
%run_matlab_sims Run simulations on generated data
%   N - number of replications
%   data_dir - UNIX directory where data lives

%rng(myseed);
rand('seed',myseed)
    for i = 1:N
       infile = fullfile(data_dir, sprintf('input_%06d.mat', i));
       outfile = fullfile(data_dir, sprintf('phmpl_ic_ms_output_%06d.mat', i));
       fprintf('Processing %s... -> %s\n', infile, outfile);
       load(infile);
       
       endpts = common_ymin_ymax(N, data_dir);
       
       ph_out = phmpl_ic_ms_mod(data.y, data.X, data.smooth, endpts);
       b = ph_out.b;
       bh = ph_out.bh;
       bhmod = ph_out.bhmod;
       the = ph_out.the;
       bch = ph_out.bch;
       bchmod = ph_out.bchmod;
       bs = ph_out.bs;
       bsmod = ph_out.bsmod;
       lam = ph_out.lam;
         
       save('-6', outfile, 'b', 'bh', 'bhmod', 'the', 'bch', 'bchmod', 'bs', 'bsmod', 'lam')
    end
    fprintf('Processed %d files.\n', N);
end

