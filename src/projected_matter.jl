
function w_ell_tullio(c,T)
    return @tullio w[i,j,k] := c[j,k,l] * T[i,j,k,l]
end