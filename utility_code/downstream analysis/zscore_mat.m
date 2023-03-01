function Mz=zscore_mat(M)

Mz=(M-nanmean(M(:)))/std(M(:));