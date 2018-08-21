used_names = []
def check_name(name):
    if name in used_names:
        raise Exception(name + " is already used!!")
    used_names.append(name)

def one_index_names(name):
    check_name(name)
    return [name +"_"+ i  for i in indexes]
    
def two_index_names(name):
    check_name(name)
    return [name +"_"+ i + j for i in indexes for j in indexes if i != j]

def two_index_names_repeat(name):
    check_name(name)
    names = []
    for i in indexes:
        for j in indexes:
            if i != j:
                names.append(name + "_" + i + j)
            else:
                for k in indexes:
                    if k != i:
                        T = [i,k]
                        T.sort()
                        names.append(name + "_" + i + i + "_" + "".join(T))
    return names
    #return [name +"_"+ i + j for i in indexes for j in indexes]

def three_index_names(name):
    check_name(name)
    return [name +"_"+ i + j + k for i in indexes for j in indexes for k in indexes if i != j and i != k and j !=k]

def jjpS_names(name):
    check_name(name)
    names = []
    for S in Subsets(indexes, g):
        Ss = list(S)
        Ss.sort()
        for jj in Subsets(indexes,2):
            jjs = list(jj)
            jjs.sort()
            j = jjs[0]
            jp = jjs[1]
            
            if j not in S and jp not in S:
                names.append(name + "_" + j + jp + "_" + "".join(Ss))
                names.append(name + "_" + jp + j + "_" + "".join(Ss))

    return names

def ijS_names(name):
    check_name(name)
    names = []
    for S in Subsets(indexes, g):
        Ss = list(S)
        Ss.sort()
        for i in Ss:
            for j in indexes:
                if j not in Ss:
                    names.append(name + "_" + i + j + "_" + "".join(Ss))

    return names

def iS_names(name):
    check_name(name)
    names = []
    for S in Subsets(indexes, g):
        Ss = list(S)
        Ss.sort()
        for i in Ss:
            names.append(name + "_" + i + "_" + "".join(Ss))
    return names

def inotinS_names(name):
    check_name(name)
    names = []
    for S in Subsets(indexes, g):
        Ss = list(S)
        Ss.sort()
        for i in indexes:
            if i not in Ss:
                names.append(name + "_" + i + "_" + "".join(Ss))
    return names

def ijlS_names(name):
    check_name(name)
    names = []
    for S in Subsets(indexes, g):
        Ss = list(S)
        Ss.sort()
        for i in S:
            for l in S:
                for j in indexes:
                    if not j in S:
                        names.append(name + "_"+ i + j +l + "_" + "".join(Ss))
    return names
                    

        
    
Uring = PolynomialRing(QQ, two_index_names("al") + two_index_names_repeat("ga") +  two_index_names("A") + two_index_names("D") + two_index_names("t") + two_index_names("v") + two_index_names("R") + two_index_names_repeat("del") + two_index_names("B") + two_index_names("be") + two_index_names_repeat("ep") + two_index_names("psi")  + two_index_names("u") + two_index_names("g") + iS_names("pi") + two_index_names("k") + iS_names("S") 
                      + jjpS_names("c") + jjpS_names("cc") + ijS_names("a") + ijS_names("ac") + ijS_names("alin") + ijS_names("asq") + jjpS_names("hshsc") + ijS_names("b") +ijlS_names("d") + ijS_names("fhsc") + ijS_names("e") + ijS_names("r") + iS_names("s") +ijS_names("hhsc")
                      + two_index_names_repeat("theta") + two_index_names_repeat("Pi") + inotinS_names("x") + ["aaa","bbb","ccc","ddd","eee"])

Rring = PolynomialRing(Uring, iS_names("f") + iS_names("h")+ iS_names("H"))
    
#Udict = Uring.gens_dict()
Rdict = Rring.gens_dict()

U_rel_dict = dict()
#e.g. a_12 = b_12 + c_12 is represented by {"a": ("12", "b_12 + c_12")}

def change_indexes(matchobj, old_to_new):
    new_indexes = ""
    for i in matchobj.group(2):
        new_indexes += old_to_new[i]
    if matchobj.group(3) is None:
        #return "v('" + matchobj.group(1) + "', '" +  "','".join(new_indexes) + "')"
        return matchobj.group(1) + "_" + "".join(new_indexes)

    #else
    new_indexes2 = []
    for i in matchobj.group(3):
        new_indexes2.append(old_to_new[i])
        
    new_indexes2.sort()

    #return "v('" + matchobj.group(1) + "', " + str(new_indexes) + "', '" + ", '".join(new_indexes2) + "')"
    return matchobj.group(1) + "_" + new_indexes + "_" + "".join(new_indexes2)
    #return "v('" + matchobj.group(1) + "', '" +  "','".join(new_indexes) + "', " + str(new_indexes2) + ")"
                

def apply_known_relations(c):
    """
    c should be an element of Uring.
    """
    subs_dict = dict()
    
    for v in c.variables():
        vsplit = str(v).split("_")
        vname = vsplit[0]
        if vname in U_rel_dict.keys():
            stored_indexes, stored_rel = U_rel_dict[vsplit[0]]
            v_indexes = "".join(vsplit[1:])
        
            stored_to_new = dict(zip(stored_indexes,v_indexes))
            
            
            subs_dict[v] = sage_eval("Uring('{0}')".format(re.sub(r"([a-zA-Z]+)_(\d+)(?:_(\d+))?", lambda m: change_indexes(m,stored_to_new), stored_rel)),globals())

    return c.subs(subs_dict)
        
        
import re
def ry(expr_str, **vindexes):
    
    pos_ind = "[\d" + "".join(vindexes.keys()) + "]"
    vindexes.update({str(d):str(d) for d in indexes})
    
    expr = sage_eval("Rring('{0}')".format(re.sub(r"([a-zA-Z]+)_("+ pos_ind + r"+)(?:_("+ pos_ind + "+))?", lambda m: change_indexes(m,vindexes), expr_str)), globals())
    
    result = expr.map_coefficients(apply_known_relations)
    
    if result in Uring:
        return Uring(result)
    #else it in in Rring
    return result


#the relations
def ff(i,j):
    return ry("f_i_ij*f_j_ij",i=i,j=j), ry("al_ji*h_i_ij + al_ij*h_j_ij + ga_ji*f_i_ij + ga_ij*f_j_ij + A_ij",i=i,j=j)

def fh(i,j):
    return ry("f_i_ij*h_j_ij",i=i,j=j), ry("D_ij*f_j_ij^2 + t_ji*h_i_ij + v_ij*h_j_ij + R_ji*f_i_ij + del_ij*f_j_ij + B_ij",i=i,j=j)

def hh(i,j):
    return ry("h_i_ij*h_j_ij",i=i,j=j), ry("be_ji*f_i_ij^2 + be_ij*f_j_ij^2 + ep_ji*h_i_ij + ep_ij*h_j_ij + psi_ji*f_i_ij + psi_ij*f_j_ij + u_ij",i=i,j=j)

def h2(i,j):
    #Note we are normalizing by q=r=u=0.
    return ry("h_i_ij^2",i=i,j=j), ry("f_i_ij^3 + g_ij*h_j_ij + pi_i_ij*f_i_ij + k_ij*f_j_ij + s_i_ij",i=i,j=j)