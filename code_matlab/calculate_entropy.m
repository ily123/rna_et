function entropy = calculate_entropy(int_matrix)
    entropy = [];
    [h,~] = size(int_matrix);
    int_matrix(:,int_matrix(h,:)==16) = [];
    [~,w] = size(int_matrix);
    for i=1:w
        column = int_matrix(:,i); % there is probably a premade library for this
        counts = zeros(1,max(column));
        for j=1:length(column) %get proportion of each element
            base = column(j);
            counts(base) = counts(base) + 1;
        end
        counts_ = counts;
        counts_(counts_==0) = [];
        counts_/sum(counts_);
        ent = -1*sum(counts_.*log(counts_/sum(counts_)));
        entropy = [entropy; ent];
    end