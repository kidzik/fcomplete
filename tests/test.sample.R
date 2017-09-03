res = fcomplete:::fc.sample(ftrue)
res$test.mask

A = ftrue
A[!is.na(res$test)] = (A - res$test)[!is.na(res$test)]
A[!is.na(res$train)] = (A - res$train)[!is.na(res$train)]
