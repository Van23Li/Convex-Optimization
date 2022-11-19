function [X_c,y_c] = get_comparable_data(y,X_observed, y_observed)
%GET_COMPARABLE_DATA Build matrix of y_j and X_j, where (i,j)in E
comparable_y = max(y_observed(y_observed<y));
idx = randi([1,length(y_observed(y_observed==comparable_y))]);
...
X_c = X_observed(y_observed==comparable_y,:);
X_c = X_c(idx,:);
y_c = y_observed(y_observed==comparable_y);
y_c = y_c(idx);
end

