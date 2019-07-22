function [aln, pdb_res, pdb_seq_prof] = profilex(atoms, pdb, chain, format_flag, aln_fp)
  % load alignment
  loaded = load(aln_fp);
  aln = loaded.aligned;
  [~,pdb_res,pdb_seq_prof,~] = profile_pdbcif_chain_to_aln(aln,...
                                  atoms, chain, format_flag);
  % add profiled sequence to alignment
  nseq = length(aln);
  aln(nseq+1).Sequence = pdb_seq_prof;
  aln(nseq+1).Header = 'pdb_query';