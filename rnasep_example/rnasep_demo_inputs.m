function main(pdb, chain, aln_fp, fs_fp)
  addpath '../code_matlab/'
  addpath '../code_matlab/plot_functions'

  sleep_time = 0
  % input data
  %pdb = '3q1q'; 
  %chain = 'B';
  %aln_fp = 'input_data/RF00010.mat';
  %fs_fp = 'input_data/3Q1Q_RF00010_functional_site_COMMENTS.txt';

  % load pdb 
  %fp_pdb = sprintf('input_data/%s.mat', pdb);
  %loaded = load(fp_pdb);
  %structure = loaded.structure;

  %load alignment (MSF) and save it as .mat
  aligned = multialignread(aln_fp);
  aln_fp = 'input_data/RF00010.mat';
  save(aln_fp, 'aligned');

  % FETCH PDB FILE FROM PDB WEBSITE DIRECTLY 
  structure = getpdb(pdb);
  
  % calculate 4A adjacency matrix for the structure
  fprintf('\n');
  disp(sprintf('[1] Calculating matrix for %s chain %s', pdb, chain));
  atoms = struct2table(structure.Model(1).Atom);
  [adj_matrix, seq_table] = pdbtable2distmatrix(atoms, chain);
  disp('...done');
  fprintf('\n');

  % Profile sequence in PDB to alignment 
  fprintf('\n');
  disp(sprintf('[2] Profile %s chain %s to alignment %s', pdb, chain, aln_fp));
  format_flag = 'pdb';
  [aln, pdb_res, pdb_seq_prof] = profile(atoms, pdb,...
                                   chain, format_flag, aln_fp);
  profiled_msf = sprintf('data/profiled_3q1q.msf');
  multialignwrite(profiled_msf, aln, 'Format', 'MSF');
  disp(sprintf('profiled alignment saved to %s', profiled_msf));
  disp('...done');
  fprintf('\n');
  pause(sleep_time)

  % Trace the profiled alignment
  fprintf('\n');
  disp(sprintf('[3] Trace profiled alignment', pdb, chain, aln_fp));
  fprintf('\n');
  trace_result = et_wetc_wrapper(profiled_msf, 'rnase_p_example',...
                                 'realval', 'pdb_query');
  disp('...done');
  fprintf('\n');

  % trim matrix to account for traced residues only
  trace_res = pdb_res;
  matrix_res = seq_table.pdb_resSeq;
  try
    keep = find(ismember(matrix_res, trace_res));
  catch
    matrix_res = str2double(seq_table.pdb_resSeq);
    keep = find(ismember(matrix_res, trace_res));
  end
  trace_matrix = adj_matrix(keep,keep);

  % Calculate clustering z-score
  % Calculate expected clustering table
  % z-score code only needs half the matrix (off diagonal)
  lower_trig = trace_matrix;
  lower_trig = tril(trace_matrix, -1);
  exp_table = calc_expected_clust(lower_trig);
  % calculate clustering z-score
  clust_z = calc_clust_z_score(lower_trig, trace_result, exp_table);
  clust_z.z_score(end) = 0;
  if clust_z.m(1) == 1
     clust_z.z_score(1) = 0;
  end
  % get ET smoothness too
  smoothness = rank_smoothness(trace_result, lower_trig);

  % Calculate overlap z-scores
  % access active sites (pdb labels of active residues)
  %fs_fp = 'input_data/3Q1Q_RF00010_functional_site_COMMENTS.txt';
  fs = importdata(fs_fp);
  fs = fs.data;
  % find these residues in the trace residue list
  trace_fs = find(ismember(trace_res, fs));
  % calculate ET overlap
  olap_z = calc_overlap_z_score(trace_result.coverage, trace_fs);


  % calculate and print summary z-scores
  %[oz_max, oz_mean, oz_max35, oz_mean35] = zscore_stats(olap_z, true);
  %[cz_max, cz_mean, cz_max35, cz_mean35] = zscore_stats(clust_z, true);
  %out = [cz_max, cz_mean, cz_max35, cz_mean35];
  
  % print results as XSL
  writetable(clust_z, 'results/rnasep_3q1q_clustering.xls')
  writetable(olap_z, 'results/rnasep_3q1q_overlap.xls')

  % print figures
  plot_8cm(clust_z.cov_bin, clust_z.z_score,'clustering z-score',...
           'visualization/rnaseo_3q1q_clusteringz.png') 
  plot_8cm(olap_z.cov_bin, olap_z.z_score,'overlap z-score',...
           'visualization/rnaseo_3q1q_overlapz.png') 

function [adj_matrix, seq_table] = make_adj_matrix(pdb, chain) 

  % check if cif or pdb 
  
  % plot the matrix, save the plot
  h = imagesc(adj_matrix);
  fp = sprintf('../visualization/distance_matrix_%s_%s.png',...
               pdb, chain);
  saveas(h, fp);
  clear h;

function [aln, pdb_res, pdb_seq_prof] = profile(atoms, pdb, chain, format_flag, aln_fp)
  % load alignment
  loaded = load(aln_fp);
  aln = loaded.aligned;
  [~,pdb_res,pdb_seq_prof,~] = profile_pdbcif_chain_to_aln(aln,...
                                  atoms, chain, format_flag);
  % add profiled sequence to alignment
  nseq = length(aln);
  aln(nseq+1).Sequence = pdb_seq_prof;
  aln(nseq+1).Header = 'pdb_query';
