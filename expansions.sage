Tring.<t> = LaurentSeriesRing(Uring)

def expand_at(expr, point):
    #print "entering expand_at", expr, point
        
    if expr in Uring:
        return Tring(expr)
    #else it should be in Rring
    
    expand_dict = dict()
    

    for v in expr.variables():
        vsplit = str(v).split("_")
        vind = vsplit[1]
        
        S = vsplit[2][0], vsplit[2][1]
        
        if S[0] == vind:
            not_vind = S[1]
        elif S[1] == vind:
            not_vind = S[0]
        else:
            not_vind = None
        
        if vsplit[0] == "f":
            if vind is None:
                raise Exception("vind not in S!")
            
            if vind == point:
                #expand_dict[v] = ryt("t^0  + ga_ii_ij*t^2 + del_ii_ij*t^3 + Pi_ii_ij*t^4",i=vind,j=not_vind)/t^2 + O(t^3) #Annnoying workaround to deal with bug in LaurentSeries parser.
                expand_dict[v] = ryt(t, -2, "1", "0", "GA_ii_ij", "DEL_ii_ij", "PI_ii_ij", i = vind, j = not_vind)
            elif point in S:
                #expand_dict[v] = ryt("al_ij*t^-1 + ga_ij + del_ij*t + pi_ij*t^2 + O(t^3)",i=vind,j=not_vind)
                expand_dict[v] = ryt(t, -1, "al_ij", "ga_ij", "del_ij", "Pi_ij",i=vind,j=not_vind) 
            else: #point is not in S
                #expand_dict[v] = ryt("b_ik_ij + O(t)", i=vind, j=not_vind, k=point)
                expand_dict[v] = ryt(t,0, "b_ik_ij", i=vind, j=not_vind, k=point)
                
        elif vsplit[0] == "h":
            if vind is None:
                raise Exception("vind not in S!")
                
            if vind == point:
                #expand_dict[v] = ryt("t^-3  + ep_ii_ij + theta_ii_ij*t + O(t^2)",i=vind,j=not_vind)
                expand_dict[v] = ryt(t,-3,"1","0","0","EP_ii_ij","THETA_ii_ij",i=vind,j=not_vind)
            elif point in S:
                #expand_dict[v] = ryt("be_ij*t^-1 + ep_ij + theta_ij*t + O(t^2)",i=vind,j=not_vind)
                expand_dict[v] = ryt(t,-1,"be_ij","ep_ij","theta_ij",i=vind,j=not_vind)
            else: #point is not in S
                #expand_dict[v] = ryt("e_ik_ij + O(t)", i=vind, j=not_vind, k=point)
                expand_dict[v] = ryt(t, 0, "e_ik_ij",i=vind, j=not_vind, k=point)
            
        elif vsplit[0] == "H":
            if vind == point:
                #expand_dict[v] = ryt("t^-1 + x_k_ij", k = vind, i = S[0], j=S[1])
                expand_dict[v] = ryt(t,-1, "1", "x_k_ij",k = vind, i = S[0], j=S[1])
            elif point in S:
                #expand_dict[v] = ryt("a_ij_jk*t^-1 + ac_ij_jk + alin_ij_jk*t + asq_ij_jk*t^2 + O(t^3)", i = point, j = vind, k = not_vind)
                expand_dict[v] = ryt(t,-1,"a_ij_kl", "ac_ij_kl", "alin_ij_kl", "asq_ij_kl", i = point, j = vind, k = S[0], l = S[1])
            else:
                #expand_dict[v] = ryt("c_ij_kl",i=vind,j=point, k = S[0], l =S[1])
                expand_dict[v] = ryt(t,0,"c_ij_kl",i=vind,j=point, k = S[0], l =S[1])
                
    return Tring(expr.subs(expand_dict))
                
        
def find_rels_by_expansion(eqn, point_index,):
    if type(eqn) == tuple:
        if not len(eqn) == 2:
            raise Exception("BAD input!!!!")
        eqn = eqn[0] - eqn[1]
    
    expansion = expand_at(eqn, point_index)

    
    print "Some relations"
    for value in expansion.coefficients():
        add_rel_to_dict(value, U_rel_dict)      