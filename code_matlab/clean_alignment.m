function [clean_aln, remove] = clean_alignment(aln, ref_len)
%remove non-canonical sequences and fragments/dupes

	% remove sequences that are too long or too short
	remove = [];
	int_aln = [];%zeros(length(aln(1).Sequence),1);

	for i=1:length(aln)
		% get full sequence
		full_seq = aln(i).Sequence;
		full_seq = strrep(full_seq, '.', '-');
		aln(i).Sequence = full_seq;

		% find length of full sequence
		seq = strrep(full_seq, '-', '');
		if (length(seq) > ref_len*1.15)
			remove = [remove; i];
		elseif (length(seq) < ref_len*0.3)
			remove = [remove; i];
		end
		%keep track of gaps
		int_aln = [int_aln; nt2int(full_seq)];
	end

	clean_aln = aln;
	clean_aln(remove) = [];
	int_aln(remove,:) = [];

	int_aln(int_aln == 16) = 0;

	all_gap = find(sum(int_aln,1) == 0);

	for i=1:length(clean_aln)
		seq = clean_aln(i).Sequence;
		seq(all_gap) = [];
		clean_aln(i).Sequence = seq;
	end




	%int_aln(:,1)
	%int_aln(:,2)



	% remove all-gap columns
	%remove_columns = [];
	%counts = zeros(length(clean_aln),1);
	%for i=1:length(clean_aln)
	%seq = clean_aln(i).Sequence;
	%	seq = strrep(seq, '.', '-');
	%end
	%clean_aln = aln;
	%clean_aln(remove) = [];