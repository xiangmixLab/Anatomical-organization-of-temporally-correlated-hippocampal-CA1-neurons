function C=interleave_matrix(A,B)

for ii = numel(A):-1:1
    C(2*ii-[0 1]) = [B(ii),A(ii)];
end