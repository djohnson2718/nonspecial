def transition(H0gens, f, coefs_to_check):
    """
    This function will find the coefficients of the expression of f in the basis given H0gens. 
    coefs_to_check is a list of pairs (point, exponent) of terms in the expansion to use to 
    determine the coefficients.
    
    Note: The user is responsible for making sure that the parameters are compatible!!
    """
    H0gens = [ry(g) for g in H0gens]
    f = ry(f)
    
    M = Matrix(Uring, len(coefs_to_check),len(H0gens),)
    b = vector(Uring,len(coefs_to_check))

    for i, (point, expon) in enumerate(coefs_to_check):
        b[i] = expand_at(f,point)[expon]
        for j,gen in enumerate(H0gens):
            M[i,j] = expand_at(gen,point)[expon]
            
    result = M.augment(b).echelon_form()
    for row in result.rows():
        print row
        
    return result