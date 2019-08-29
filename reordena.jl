function reordena(x,CDC,DESL_PR)

# Function to reorder de vectors in accordance with boundary conditions

n_el_total = length(CDC[:,1])
n_nos =2*n_el_total
n_des_pr = length(DESL_PR[:,1])
cdc_val=zeros(6*n_el_total,1)
De=zeros(n_nos,3)
T=zeros(n_el_total,7)

# Generation of vector of boundady conditions
for el = 1 : n_el_total
    for no = 1 : 3
        for gdl = 1 : 2
            cdc_val[6*(el-1)+2*(no-1)+gdl,1] = CDC[el,4*no+2*gdl-3]
        end
    end
end

for de=1 : n_des_pr
    n_no=DESL_PR[de,1]
    n_gdl=DESL_PR[de,2]
    n_el=DESL_PR[de,3]
    n_no_loc=DESL_PR[de,4]
    ind_x=(2*n_no-2)+n_gdl
    ind_cd=(6*n_el+2*n_no_loc-8)+n_gdl
    troca = cdc_val[ind_cd]
    cdc_val[ind_cd] = x[ind_x]
    valor_x=x[ind_x]
    x[ind_x] = troca
    if(DESL_PR[de,5]!=0)
        n_el=DESL_PR[de,5]
        n_no_loc=DESL_PR[de,6]
        ind_cd=(6*n_el+2*n_no_loc-8)+n_gdl
        cdc_val[ind_cd]=valor_x
    end
end

# Set the values of boundary displacements [De], internal nodes
# displacements [Di] & boundary tractions [T]

ndof=2*n_nos
desl = x[1:ndof,1]
trac = cdc_val

for no = 1 : n_nos
    De[no,1] = no
    De[no,2] = x[2*no-1]
    De[no,3] = x[2*no]
end

for el = 1 : n_el_total
    T[el,1] = el
    for no = 1 : 3
        for gdl = 1 : 2
            T[el,2*(no-1)+gdl+1] = cdc_val[6*(el-1)+2*(no-1)+gdl]
        end
    end
end
return desl,trac,De,T
end
