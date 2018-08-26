lms = []
#the leading monomials

#This only workds for genus 2
for S in Subsets(indexes,2):
    #Ss = list(S)
    #Ss.sort()
    for j in indexes:
        if j in S:
            continue
        for i in S:
            if i == S[0]:
                k = S[1]
            else:
                k = S[0]
            #lms.append((v("f",i,Ss)*v("h",Ss,j,Ss), fhS, (Ss,i,j)))
            lms.append((ry("f_i_ik*H_j_ik",i=i,k=k,j=j), fH, (i, k, j)))
            
            #lms.append((v("h",i,Ss)*v("h",Ss,j,Ss), hhS, (Ss,i,j)))
            lms.append((ry("h_i_ik*H_j_ik",i=i,k=k,j=j), hH, (i, k, j))) 
            
        for jp in indexes:
            if jp in S or jp == j:
                continue
            #lms.append((v("h",Ss,j)*v("h",Ss,jp), hShS, (Ss,j,jp)))
            lms.append((ry("H_j_kl*H_p_kl", j=j,p=jp,k=S[0],l=S[1]), HH, (j,jp,S[0],S[1])))

    
    #lms.append( (v("f",i,S)*v("f",j,S), make_rel , ("ff", i, j,S)))
    lms.append((ry("f_i_ij*f_j_ij", i=S[0], j=S[1]), ff, (S[0],S[1])))
    #lms.append( (v("h",i,S)*v("h",j,S), make_rel, ("hh", i, j,S)))
    lms.append((ry("h_i_ij*h_j_ij", i=S[0], j=S[1]), hh, (S[0],S[1])))

    for i in S:
        if i == S[0]:
            k = S[1]
        else:
            k = S[0]
            
        #lms.append( (v("h",i,S)^2, make_rel , ("h2", i, None, S)))  
        lms.append((ry("h_i_ik^2", i=i,k=k), h2, (i, k)))

        #lms.append( (v("f",i,S)*v("h",j,S), make_rel,("fh", i, j,S)))
        lms.append((ry("f_i_ik*h_k_ik",i=i,k=k), fh, (i,k)))
        
        
   

def doS(r1,r2):
    """
    Implements the Buchberger algorithm to find relations between the coefficients.
    """
    L = lcm(r1[0].change_ring(QQ), r2[0].change_ring(QQ))
    factor1 = L.quo_rem(r1[0])[0]
    factor2 = L.quo_rem(r2[0])[0]
    S = factor1 * r1[1] - factor2 * r2[1]
    newS = reduce_lm(S)

    print "Relations Found:"
    #print newS
    
    for c in newS.coefficients():
        yield c
        
def reduce_lm(fff):
    """
    Outsource this code from doS.
    """
    fff1 = fff
    done = False
    
    while not done:
        new_fff = 0
        for m in fff1.monomials():
            can_be_reduced = False
            for lm in lms:
                if lm[0].divides(m):
                    can_be_reduced = True
                    #newS += make_rel(lm[1], lm[2], lm[3])[1] *(m/lm[0]).numerator()*S1.monomial_coefficient(m)
                    new_fff += (lm[1](*lm[2]))[1] * (m/lm[0]).numerator()*fff1.monomial_coefficient(m)
                    break
            if not can_be_reduced:
                new_fff += m*fff1.monomial_coefficient(m)


        if fff1 == new_fff:
            done = True
        else:
            #print "need to go again"
            fff1 = new_fff
    return new_fff

def eliminate_vars_buch(r1,r2):
    for rrr in doS(r1,r2):
        add_rel_to_dict(rrr,U_rel_dict,True)