
function [adjmat, seq_table] = pdbtable2distmatrix(atoms, chain)
% pdbtable2distmatrix  Convert atoms table of the chain into distance matrix 
%                      Default distance is 4A
%
%                      Output is matrix A
%                      where A(i,j) = 1 if residue i and j are within
%                      4A of each other (and 0 otherwise)
%
%                      The secound output is a sequence table to 
%                      keep track of residues

  % we only want chain 'A'
  chain_index = strcmp(atoms.chainID, chain);
  chain_atoms = atoms(chain_index, :);
  % remove inserted atoms     
  %inserts = find(cellfun(@(x) ~isequal('', x), chain_atoms.iCode));
  %chain_atoms(inserts, :) = [];
  chain_len = length(unique(chain_atoms.resSeq));
  res = unique(chain_atoms.resSeq);
  adjmat = zeros(chain_len, chain_len);

  coords = [chain_atoms.X, chain_atoms.Y, chain_atoms.Z];
  d = pdist(coords);
  d = squareform(d);
  d(d<=4) = 1;
  d(d>4) = 0;
  [ii, jj] = find(d);
  
  hardindex = [];
  sequence = {};
  pdb_resSeq = [];
  resseq_last = 'none';
  hrdinx = 0;
  hardinx2 = [];
  for inx = 1:height(chain_atoms)
    resname = chain_atoms.resName(inx);  
    resseq = chain_atoms.resSeq(inx);
    if resseq ~= resseq_last
      hrdinx = hrdinx + 1;
      hardindex = [hardindex; hrdinx];
      sequence = [sequence; resname];
      pdb_resSeq = [pdb_resSeq; resseq];
    end
    hardinx2 = [hardinx2, hrdinx];
    resseq_last = resseq;
    
  end

  seq_table = table(hardindex, pdb_resSeq, sequence);

  cluster_ii = hardinx2(ii)';
  cluster_jj = hardinx2(jj)';

  [ii, jj, cluster_ii, cluster_jj];

  ind = sub2ind(size(adjmat), cluster_ii, cluster_jj);

  adjmat(ind) = 1;
