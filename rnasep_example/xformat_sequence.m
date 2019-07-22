function matrix_aln = format_sequence(aln)
    matrix_aln = [];
    for i=1:length(aln)
        int_seq = nt2int(strrep(aln(i).Sequence, '.', '-'));
        matrix_aln = [matrix_aln; int_seq];
    end