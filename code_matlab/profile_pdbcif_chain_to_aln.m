

function [keep_track, keep_track2, pdb_seq_profiled, atoms_chain] = main(alignment, atoms, chain_id, format_flag)
%   profile_pdbcif_chain_to_aln
%   Profiles PDB sequence of a given chain to alignment.
%   Input  is (1) PDB strcture (from pdbread)
%             (2) chain you want to trace
%             (3) alignment path, MUST ME MSF with '.' for gaps
%             (4) 'cif' flag for .cif file, and 'pdb' flag for .pdb structure
%
%   Output is (1) cell with pdb sequence and pdb res nums (not all of these align)
%             (2) array of pdb res nums that WERE ALIGNED and traced
%             (3) trace struct with ET scores (corresponding pdb resmuns are in (2))
%             (4) atoms? not sure why 

  % we only want chain 'A'
  chain_index = strcmp(atoms.chainID, chain_id);
  atoms_chain = atoms(chain_index, :);
  % find and remove INSERTS (their iCode is not blank)
  if format_flag == 'pdb'
    rows = (strcmp(atoms_chain.iCode, ''));
    atoms_chain = atoms_chain(rows, :);
  elseif format_flag == 'cif'
    rows = (strcmp(atoms_chain.iCode, '?'));
    atoms_chain = atoms_chain(rows, :);
  end
  
  % remove any residues with negative index numbers
  atoms_chain(atoms_chain.resSeq < 1, :) = [];

  resSeq = atoms_chain.resSeq;
  resName = atoms_chain.resName;

  % [resSeq, atoms_chain.resName]

  unique_resSeq = unique(resSeq);
  keep_track = {};
  for i=1:length(unique_resSeq)
    % get the number by which the residue goes in the PDB
    residue_pdb_number = unique_resSeq(i);
    % find the indecies of the residue in resSeq array
    residue_indecies = find(resSeq==residue_pdb_number);
    % using the indecies access the residue identity
    residue_names = resName(residue_indecies);
    % all of these should be the same one name (A, or U, or G, or C)
    residue_name = unique(residue_names);
    if length(residue_name) > 1
      error('[ERROR 1] more than one nucleic acid assigned to residue num');
    end
    residue_name = residue_names{1};
    keep_track_row = {residue_pdb_number, residue_name};
    keep_track = [keep_track; keep_track_row];
  end
  kt = keep_track';
  sequence = strjoin(kt(2,:), '');

  % now take the sequence, and profile it to the alignment

  % pairwise align pdb sequence to every sequence in alignment
  bit_scores = [];
  for i=1:length(alignment)
    sequence_from_aln = alignment(i).Sequence;
    sequence_from_aln_no_gaps = strrep(sequence_from_aln, '.','');
    [score, alignment2seq] = nwalign(sequence, sequence_from_aln_no_gaps, 'ALPHABET', 'NT',...
        'GapOpen', 8, 'ExtendGap', 0.5);
    alignment2seq;
    score;
    bit_scores = [bit_scores; score];
  end

  % choose sequence from alignment that aligns to pdb sequence with highest bit score
  highest_bit_index = find(bit_scores == max(bit_scores));

  if length(highest_bit_index) > 1
    highest_bit_index = highest_bit_index(1);
    disp('[WARNING] multiple higest bit');
  end

  % withdraw this sequence from the alignment
  best_bit_sequence_ori = alignment(highest_bit_index).Sequence;
  % align it to the pdb (again)
  best_bit_sequence = strrep(best_bit_sequence_ori, '.', '');
  [score, alignment2seq] = nwalign(sequence, best_bit_sequence, 'ALPHABET', 'NT');

  % retrieve sequences from the alignment
  pdb_seq_aln = alignment2seq(1,:);
  best_bit_aln = alignment2seq(3,:);
  % iterate over the residue number array (keep track), and insert
  % gaps to reflect the aligned PDB sequence (need this for later)
  len_aln = length(pdb_seq_aln);
  keep_track2 = ones(len_aln, 1);
  gaps_in_aln_pdb = find(pdb_seq_aln=='-');
  keep_track2(gaps_in_aln_pdb) = 0;
  blahx = find(keep_track2==1);
  res_num = cell2mat(keep_track(:,1));
  keep_track2(blahx) = res_num;

  % find gaps in best_bit_aln, remove these columns from aligned PDB seq
  remove_gaps = find(best_bit_aln=='-');
  pdb_seq_aln(remove_gaps) = [];
  keep_track2(remove_gaps) = [];
  keep_track2 = keep_track2(keep_track2>0);

  %AAXXXAA--XX
  %--XXX--BBXX
  % becomes
  %XXX--XX
  %XXXBBXX
  
  % now iterate over the best bit sequence in the alignemtn and
  % introduce alignment gaps into the profiled PDB sequence
  pdb_seq_profiled = '';
  non_gap_count = 0;
  for i=1:length(best_bit_sequence_ori)
    letter_in_ori = best_bit_sequence_ori(i);
    if letter_in_ori == '.' % line source of errs depending on gap format
      pdb_seq_profiled = [pdb_seq_profiled, '-'];
    else
      non_gap_count = non_gap_count + 1;
      pdb_seq_profiled = [pdb_seq_profiled, pdb_seq_aln(non_gap_count)];
    end
  end

  % best_bit_sequence_ori
  pdb_seq_profiled = strrep(pdb_seq_profiled, '-', '.');

  %[pdb_seq_profiled', best_bit_sequence_ori']
