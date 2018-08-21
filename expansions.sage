def expand_at(expr, point):
    if expr in Uring:
        return expr
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
                expand_dict[v] = ryt("t^-2  + ga_ii_ij + del_ii_ij*t + Pi_ii_ij*t^2 + O(t^3)",i=vind,j=not_vind)
            elif point in S:
                expand_dict[v] = ryt("al_ij*t^-1 + ga_ij + del_ij*t + pi_ij*t^2 + O(t^3)",i=vind,j=not_vind)
            else: #point is not in S
                expand_dict[v] = ryt("b_ik_ij + O(t)", i=vind, j=not_vind, k=point)
                
        elif vsplit[0] == "h":
            if vind is None:
                raise Exception("vind not in S!")
                
            if vind == point:
                expand_dict[v] = ryt("t^-3  + ep_ii_ij + theta_ii_ij*t + O(t^2)",i=vind,j=not_vind)
            elif point in S:
                expand_dict[v] = ryt("be_ij*t^-1 + ep_ij + theta_ij*t + O(t^2)",i=vind,j=not_vind)
            else: #point is not in S
                expand_dict[v] = ryt("e_ik_ij + O(t)", i=vind, j=not_vind, k=point)
            
        elif vsplit[0] == "H":
            if vind == point:
                expand_dict[v] = ryt("t^-1 + x_k_ij", k = vind, i = S[0], j=S[1])
            elif point in S:
                expand_dict[v] = ryt("a_ij_jk*t^-1 + ac_ij_jk + alin_ij_jk*t + asq_ij_jk*t^2 + O(t^3)", i = point, j = vind, k = not_vind)
            else:
                expand_dict[v] = ryt("c_ij_kl",i=vind,j=point, k = S[0], l =S[1])
                
    return expr.subs(expand_dict)
                
        
        
        
        


expand_at_dict = dict()
def make_expand_dict():
    for j in indexes:
        expand_at_dict[j] = dict()
        for T in Subsets(indexes,2):
            Ts = list(T)
            Ts.sort()

            for i in indexes:
                if i in Ts:
                    f = v("f",i,Ts)
                    h = v("h",i,Ts)
                    if j in Ts:
                        if i != j:
                            expand_at_dict[j][f] = v("al",i,j)/t + v("ga",i,j) + v("del",i,j)*t + v("pi",i,j)*t^2 + O(t^3)
                            expand_at_dict[j][h] = v("be",i,j)/t + v("ep",i,j) + v("theta",i,j)*t + O(t^2)
                        else:
                            expand_at_dict[j][f] = 1/t^2 + v("ga",i,j,Ts) + v("del",i,j,Ts)*t + v("pi",i,j,Ts)*t^2 + O(t^3)
                            expand_at_dict[j][h] = 1/t^3 + v("ep",i,j,Ts) + v("theta",i,j,Ts)*t + O(t^2)
                    else: #j is not in Ts
                        expand_at_dict[j][f] = v("b",i,j,Ts) + O(t)
                        expand_at_dict[j][h] = v("e",i,j,Ts) + O(t)
                else: #i is not in Ts
                    hS = v("h",Ts,i)
                    if j in Ts:
                        expand_at_dict[j][hS] = v("a",j,i,Ts)/t + v("ac",j,i,Ts) + v("alin",j,i,Ts)*t + v("asq", j,i,Ts)*t^2 + O(t^3) #check the index order one more time!!!


                    else: #j not in Ts
                        if i != j:
                            expand_at_dict[j][hS] = v("c",i,j,Ts) + v("cc",i,j,Ts)*t + O(t^2)
                        else:
                            expand_at_dict[j][hS] = 1/t + v("x",i,Ts) + O(t)
                        
make_expand_dict()
                
def expand_at(func, point):
    return func.subs(expand_at_dict[point])


def find_rels_by_expansion(point_index, eqn):
    expansion = expand_at(eqn[0]-eqn[1], point_index)
    
    print "Some relations"
    for value in expansion.coefficients():
        #print value
        add_rel_to_dict(value)