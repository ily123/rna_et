function main

  addpath('../code_matlab')
  addpath('../code_matlab/plot_functions')


  % load file paths
  fps = parse_file_paths('FILE_PATHS.TXT');
  
  % load trace scores, profiled to 2WDL pdb
  trace_result = load(fps.trace_result);
  trace_result = trace_result.trace_table;

  % load pdb resnum of the trace query sequence
  trace_resid = load(fps.pdb_resid);
  trace_resid = trace_resid.large_pdbres;

  % load and parse 16S functional sites
  fs = parse_functional_sites(fps.small_subunit_fs);

  % iterate over each site, and compute overlap with ET nts
  summary35_all = table();
  zscores = [];
  fs_categ = fieldnames(fs);
  fs_categ

  names =  {'Antibiotic sites', 'Bridges', 'GAC',...
            'Modified bases', 'PTC', 'R-proteins',...
            'SRL', 'Translation Factors', 'tRNA sites', 'All sites'};

  for i=1:length(fs_categ)
    disp(sprintf('Calculating ET overlap with functional site - %s', fs_categ{i}));
    resi = fs.(fs_categ{i});
    resi = unique(resi);
    % remaining functional sites
    rem_resi = setdiff(fs.all_sites, resi);
    % map pdb resi to trace result
    resi = convert_index(trace_resid, resi);
    rem_resi = convert_index(trace_resid, rem_resi);
    % calculate overlap
    zsc_table = olap_zscore_remove_all_other_pos(resi, rem_resi,...
                  trace_result.coverage);
    
    summary35 = get_bin_summary(zsc_table, 0.35);
    % add site name to summary row 
    t1=table(string(fs_categ{i}), 'VariableNames', {'site'});
    summary35 = [t1, summary35];

    fill = true;
    [zsc_max, zsc_mean, zsc_max35, zsc_mean35] = zscore_stats(zsc_table, fill);
    
    mean_rank = mean(trace_result.coverage(resi));
    [x,y,auc] = remove_non_poi_cons_and_other_pos(resi, rem_resi, trace_result.coverage);

    zrow = [zsc_max, zsc_mean, zsc_max35, zsc_mean35,auc,...
            mean_rank, length(resi),...
            length(trace_result.coverage) - length(rem_resi) - length(resi),...
            length(rem_resi)];

    zscores = [zscores; zrow];
    summary35_all = [summary35_all; summary35];
    % plot mini panel for overlap figure
    fig_fp = [fps.fig_small, '%s_overlap_lsu.png'];
    fig_fp = sprintf(fig_fp, fs_categ{i});
    plot_3cm(zsc_table.cov_bin, zsc_table.z_score, names{i}, fig_fp);
    % auc panel for supp info
    fig_fp_auc = [fps.fig_small, '%s_AUC_lsu.png'];
    fig_fp_auc = sprintf(fig_fp_auc, fs_categ{i});
    plot_3cm_auc(x,y,auc,names{i},fig_fp_auc);
  end

  % display as table
  out = array2table(zscores, 'VariableNames', {'max' 'mean' 'max35' 'mean35'...
   'auc' 'meanrank' 'resi' 'neg' 'rem_resi'});
  site_names = array2table(fs_categ, 'VariableNames', {'site'});
  out = [site_names, out]

  display('LEGEND:');
  display('z_max     peak of the z-score curve');
  display('z_mean    average score, calculated over the entire curve');
  display('z_max35   peak of the z-score curve within 0-35% coverage bins');
  display('z_mean35  average score, calculated over the 0-35% coverage bins [reported in main text]');
  display('auc       accuracy measured using ROC [reported in supplementary]');
  display('meanrank  average ET rank of nts defining functional site');
  display('resi      # nucleotides in functional site');
  display('rem_resi  # nucleotides in all other functional sites (excluded from analysis)');
  display('neg       # nucleotides not assinged to functional site (the negatives)');
  disp(' ');
  disp(sprintf('Individual z-score and AUC curves saved to visualization/'));
  
  % display summary table for ET=35% too
  summary35_all 
  % save summary table to results
  sfp_table = [fps.results, 'results_fs_overlap_summary.xls'];
  disp(sprintf('Writing overlap summary table to %s', sfp_table));
  writetable(out, sfp_table);
  % save summary table for cov=0.35
  sfp_table2 = [fps.results, 'results_fs_ET35_summary.xls'];
  disp(sprintf('Writing ET=35% site recovery summary table to %s', sfp_table));
  writetable(summary35_all, sfp_table2);

function summary = get_bin_summary(olap_table, bin)
  %find bin closest to 35% coverage
  inx = nnz(olap_table.cov_bin <= 0.35);
  m = olap_table.m(inx);
  cov_bin = olap_table.cov_bin(inx);
  L = olap_table.m(end);
  fs_total = olap_table.observed(end);
  fs_observed = olap_table.observed(inx);
  fs_expected = olap_table.expected(inx);
  zscore = olap_table.z_score(inx);
  % back calc the std dev
  std_dev = (fs_observed - fs_expected)/zscore;
  %fs_sdev = olap_table.
  % return one row table
  fs_missed = fs_total - fs_observed;
  summary = array2table([m cov_bin fs_observed fs_missed,...
                         fs_expected std_dev zscore],...
                         'VariableNames',...
                         {'m' 'cov_bin' 'fs_obs',...
                          'fs_miss' 'fs_exp' 'std' 'zscore'});

function trace_inx = convert_index(trace_resi, resi)
  trace_inx = find(ismember(trace_resi, resi));

function fs_struct = parse_functional_sites(fp)
  
  fs_struct = struct();
  % load xlsx file and populate struct
  fs_table = readtable(fp);
  fs_table(find(strcmp(fs_table.label, 'let')),:) = [];
  fs_table(find(strcmp(fs_table.label, 'ben')),:) = [];
  fs_table(find(strcmp(fs_table.label, 'mg')),:) = [];
  labels = unique(fs_table.label);
  for i=1:length(labels)
    cat_label = labels{i};
    % find all resi that belong to category
    resi_indx = find(strcmp(fs_table.label, cat_label));
    %length(resi_indx)
    cat_resi = fs_table.resi(resi_indx);
    fs_struct.(cat_label) = cat_resi;
  end
  
  fs_struct.all_sites = unique(fs_table.resi);


function zsc_table = olap_zscore_remove_all_other_pos(poi, other_pos, scores)
  labels = zeros(length(scores),1);
  labels(poi) = 1;
  % remove non-poi sites from the negative sets (except overlapping res)
  scores(other_pos) = [];
  labels(other_pos) = [];
  % recover indices of the active site residues
  as_index = find(labels);
  zsc_table = calc_overlap_z_score(scores, as_index);

function [x,y, auc] = remove_non_poi_cons_and_other_pos(poi, other_pos, scores)
  scores = 1-scores;
  max_sc = max(scores);

  % find universally conserved residues
  cons = scores==max_sc;

  % set labels
  labels = zeros(length(scores),1);
  labels(poi) = 1;

  % remove non-poi sites from the negative sets (except overlapping res)
  aa = [scores, cons, labels];
  aa(other_pos, :) = [];
  % remove entirely conserved residues from experiment, unless poi
  aa_table = array2table(aa);
  rows = (aa_table.aa2 == 1 & aa_table.aa3 == 0);
  removei = find(rows);
  aa(removei, :) = [];
  scores = aa(:,1);
  labels = aa(:,3);
  [x,y,~,auc] = perfcurve(labels, scores, 1);
