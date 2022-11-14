# Eliminates the variable y in the system
# Ax + By = u
# Cx + Dy = v
# whenever D is invertible. 
function schur_complement(A, B, C, D, u, v)
    D_fact = lu(D)

    A_other = -(D_fact \ C)
    b_other = D_fact \ v

    A_new = A + B*A_other
    b_new = u - B*b_other

    return A_new, b_new, A_other, b_other, D_fact
end

# Eliminates all variables with indices not found in idx_schur
function schur_complement_full(M, b, idx_schur)
    m, n = size(M)

    @assert m == n 
    @assert m == length(b)
    @assert m == length(idx_schur)
    
    n_idx_schur = .~idx_schur

    A = M[idx_schur, idx_schur]
    B = M[idx_schur, n_idx_schur]
    C = M[n_idx_schur, idx_schur]
    D = M[n_idx_schur, n_idx_schur]

    u, v = b[idx_schur], b[n_idx_schur]

    return schur_complement(A, B, C, D, u, v)
end

function schur_complement(M, b, idx_schur)
    A_new, b_new, _, _, D_fact = schur_complement_full(M, b, idx_schur)
    return A_new, b_new, D_fact
end