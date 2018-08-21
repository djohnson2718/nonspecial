vars_to_save = ["al", "be", "ga", "pi", "ep"]
def add_rel_to_dict(rrr, rel_dict):
    do_monos =[]
    for m in rrr.monomials():
        print m, m.degree()
        if m.degree() == 1:
            do_monos.append(m)

    doe_monos = []
    for m in do_monos:
        if str(m).split("_")[0] not in vars_to_save:
            doe_monos.append(m)

    if len(doe_monos) == 0:
        if len(do_monos) >0:
            raise Exception("Must eliminate reserved var!!!!!")
        print "This relation cannot be solved for a single variable:", rrr
    else:
        m = doe_monos[0]
        
        msplit = str(m).split("_")
        
        m_name = msplit[0]

        m_ind = "".join(msplit[1:]) 

        
        if str(m) not in rel_dict.keys():
        
            print "Adding relations of the form:",m, "=",m-1/rrr.monomial_coefficient(m)*rrr
        
            rel_dict[m_name] = (m_ind, str(m-1/rrr.monomial_coefficient(m)*rrr))
            
        else:
            print "Variable {0} is already solved, taking no action, you should check this!!!".format(m), rrr

        reduce_dict(rel_dict)

            
def reduce_dict(rel_dict):
    to_update = dict()
    for k, val in rel_dict.items():
        ry_val1 = str(ry(val[1]))
        if ry_val1 != val[1]:
            to_update[k] = (val[0],ry_val1)
            
    rel_dict.update(to_update)