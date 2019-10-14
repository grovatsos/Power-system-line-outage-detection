function D = computeKL(post, pre, k)

    D = 1/2*(trace(pre\post) - k + log(det(pre/post)));

end