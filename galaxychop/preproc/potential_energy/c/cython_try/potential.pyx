# Desde acá, reescribo el Octree de C en Cython (!!!)

# Bruno:
# 1er intetno de transcribir el código de C a un .pyx de Cython.
# Dios y la Virgen sabrán si está bien y qué puede llegar a fallar...

cimport cython

# Para variables globales:
cdef float Thetamax = 0.45
cdef float GCONS = 6.67300e-20  # Constante de Gravitación [km³ / kg / seg²]
cdef float Msol = 1.9891e30  # Masa del Sol [Kg]
cdef float Kpc = 3.08568025e16  # Kiloparsec . Kilometro
cdef int KERNEL_LENGTH = 10000
cdef int SOFT = 3

# Para structs:
cdef struct NODE:
    float[3] s;  # center of mass
    long partind;
    float[3] center,
    float len;  # center and sidelength of treecubes
    float mass;  # mass
    float oc;  # variable for opening criter
    struct NODE *next,*sibling,*father,*suns[8];


################################################################
################################################################

cdef void force_setkernel(float *knlrad, float *knlpot):
    cdef int i;
    cdef float u;
    cdef float u2;
    cdef float u3;
    cdef float u4;
    cdef float u5;

    for i in range(KERNEL_LENGTH):
        u = <float>i/<float>KERNEL_LENGTH;

        knlrad[i] = u;

        if u <= 0.5:
            u2 = u*u;
            u4 = u2*u2;
            u5 = u4*u;
            knlpot[i]  = 16.0/3.0*u2;
            knlpot[i] -= 48.0/5.0*u4;
            knlpot[i] += 32.0/5.0*u5;
            knlpot[i] -= 14.0/5.0;
        else:
            u2 = u*u;
            u3 = u2*u;
            u4 = u2*u2;
            u5 = u4*u;
            knlpot[i]  = 1.0/15.0/u;
            knlpot[i] += 32.0/3.0*u2;
            knlpot[i] -= 16.0*u3;
            knlpot[i] += 48.0/5.0*u4;
            knlpot[i] -= 32.0/15.0*u5;
            knlpot[i] -= 16.0/5.0;

################################################################
################################################################

cdef void add_particle_props_to_node(struct NODE *no, float *pos, float mp):

    cdef int i;
    cdef float delta;

    for i in range(3):
        delta = pos[i] - no.center[i]
        no.s[i] += mp * delta

    no.mass += mp;

    return

################################################################
################################################################

cdef void force_setupnonrecursive(struct NODE *no, struct NODE **last)

    cdef int i;
    cdef struct NODE *nn;
    
    if *last:
        *last.next = no;

    *last = no;
    
    for i in range(8):
        nn = no.suns[i];
        if nn != NULL:
            force_setupnonrecursive(nn, last);


################################################################
# packs the particles of group 'gr' into into BH-trees
################################################################

cdef float force_treeevaluate_potential(const struct NODE *no, const float pos [], float *knlrad, float *knlpot)

    cdef float r2,dx,dy,dz,r,u,h,ff;
    cdef float wp;
    cdef float h_inv;
    cdef float local_pot;
    cdef int ii;

    h = 2.8*<float>SOFT;
    h_inv = 1.0 / h;

    local_pot = 0.;

    while no:
        dx = no.s[0] - pos[0];     
        dy = no.s[1] - pos[1];     
        dz = no.s[2] - pos[2];

        r2 = dx*dx + dy*dy + dz*dz;

        if no.partind >= 0:  # single particle
            r = <float>sqrt(r2);  
            
            u = r * h_inv;

            if u >= 1:
                local_pot -= no.mass/r;
            else:
                ff = u*KERNEL_LENGTH;
                ii = <int>ff;
                ff -= knlrad[ii]*KERNEL_LENGTH;
                wp = knlpot[ii] + (knlpot[ii+1] - knlpot[ii])*ff;
                
                local_pot += no.mass*h_inv*wp;
            
            no = no.sibling;
        else:
            if r2 < no.oc:
                no = no.next;  # open cell
            else:
                r = <float>sqrt(r2);  
                u = r*h_inv;
            
                if u >= 1:  # ordinary quadrupol moment
                    local_pot += -no.mass/r;
                else:  # softened monopole moment
                    ff = u*KERNEL_LENGTH;
                    ii = <int>ff;
                    ff -= knlrad[ii]*KERNEL_LENGTH;
                    wp = knlpot[ii] + (knlpot[ii+1]-knlpot[ii])*ff;

                    local_pot += no.mass*h_inv*wp;
                    
                no = no.sibling;

    return local_pot;


################################################################
################################################################


cdef int tree_potential(
    const int npart, const float *mp, const float *x, const float *y, const float *z, float *Ep
    ):
    cdef int i, j, subp,subi,p,subnode,fak;
    cdef float[3] icoord
    cdef float imass
    cdef float[3] pcoord
    cdef float pmass
    cdef float length;
    cdef float dx,dy,dz;
    cdef float[3] xmin
    cdef float[3] xmax;
    cdef struct NODE *nodes, *last, *nfree,*th,*nn,*ff;
    cdef int MaxNodes, numnodestotal;
    cdef float *knlrad, *knlpot; 
    cdef knlrad = <float*> malloc((KERNEL_LENGTH+1) * sizeof(float));
    cdef knlpot = <float*> malloc((KERNEL_LENGTH+1) * sizeof(float));

    for i in range(npart):
        Ep[i] = 0.0;

    MaxNodes = 2 * npart + 200;
    nodes = <struct NODE*> malloc(MaxNodes*sizeof(struct NODE));
    assert nodes != NULL;

    force_setkernel(knlrad, knlpot);

    nfree = nodes;
    numnodestotal = 0;

    # find enclosing rectangle
    xmin[0] = xmax[0] = x[0];
    xmin[1] = xmax[1] = y[0];
    xmin[2] = xmax[2] = z[0];

    for i in range(npart):
        xmin[0] = x[i] if x[i] < xmin[0] else xmin[0]
        xmin[1] = y[i] if y[i] < xmin[1] else xmin[1]
        xmin[2] = z[i] if z[i] < xmin[2] else xmin[2]
        xmax[0] = x[i] if x[i] > xmax[0] else xmax[0]
        xmax[1] = y[i] if y[i] > xmax[1] else xmax[1]
        xmax[2] = z[i] if z[i] > xmax[2] else xmax[2]
    
    length = xmax[0]-xmin[0]
    for j in range(3):  # determine maxmimum extension
        if (xmax[j]-xmin[j]) > xmax[0]-xmin[0]:
            length = xmax[j]-xmin[j];

    length *= 1.01
    assert length > 0

    # insert first particle in root node
    for j in range(3):
        nfree.center[j] = (xmax[j]+xmin[j])/2;

    nfree.len = length;
    
    # inicializa variables
    nfree.father = 0;
    for j in range(8):
        nfree.suns[j] = 0;

    nfree.partind = 0;
    nfree.mass = mp[0];
    pcoord[0] = x[0];
    pcoord[1] = y[0];
    pcoord[2] = z[0];
    for j in range(3):
        nfree.s[j] = nfree.mass*(pcoord[j] - nfree.center[j]);
    
    # inicializa la variable que apunta al hermano
    nfree.sibling = 0;

    numnodestotal += 1
    nfree += 1
    
    if (numnodestotal >= MaxNodes):
        print(stderr,"Maximum number %d of tree-nodes reached. file: %s line: %d\n",
                numnodestotal,__FILE__,__LINE__);
        free(nodes);
        free(knlrad);
        free(knlpot);

    # insert all other particles
    for i in range(npart):
        th        = nodes;
        icoord[0] = x[i];
        icoord[1] = y[i];
        icoord[2] = z[i];
        imass     = mp[i];
    
        while True:
            add_particle_props_to_node(th,icoord,imass);

            if th.partind >= 0:
                break;
            
            subnode = 0
            fak = 1
            for j in range(3):
                if icoord[j] > th.center[j]:
                    subnode += fak;
                fak <<= 1

            nn = th.suns[subnode];
            if nn != NULL:
                th = nn;
            else:
                break;
        
        if th.partind >= 0:  # cell is occcupied with one particle
            while True:
                p         = th.partind;
                pmass     = mp[p];
                pcoord[0] = x[p];
                pcoord[1] = y[p];
                pcoord[2] = z[p];

                subp = 0
                fak = 1
                for j in range(3):
                    if pcoord[j] > th.center[j]:
                        subp += fak;
                    fak <<= 1

                nfree.father = th;
                
                for j in range(8):
                    nfree.suns[j] = 0;

                nfree.sibling = 0;
                
                nfree.len = th.len/2;
            
                for j in range(3):
                    nfree.center[j] = th.center[j];

                for j in range(3):
                    if (pcoord[j] > nfree.center[j]):
                        nfree.center[j] += nfree.len/2;
                    else:
                        nfree.center[j] -= nfree.len/2;

                nfree.partind = p;
                nfree.mass    = pmass;
                for j in range(3):
                    nfree.s[j]  = (pcoord[j] - nfree.center[j]);
                    nfree.s[j] *= pmass;

                th.partind = -1;
                th.suns[subp] = nfree;
            
                numnodestotal += 1 
                nfree += 1
                if numnodestotal >= MaxNodes:
                    print(stderr,"Maximum number %d of tree-nodes reached. i=%d, npart=%d. file: %s line: %d\n",
                            numnodestotal,i,npart,__FILE__,__LINE__);
                    free(nodes);
                    free(knlrad);
                    free(knlpot);
                    return False; 

                subi = 0
                fak = 1
                for j in range(3):
                    if icoord[j] > th.center[j]:
                        subi += fak;
                    fak <<= 1

                if subi == subp:  # the new particle lies in the same sub-cube
                    th = nfree-1;
                    add_particle_props_to_node(th,icoord,imass);
                else:
                    break;

        subi = 0
        fak = 1
        for j in range(3):
            if icoord[j] > th.center[j]:
                subi += fak;
            fak <<= 1
        
        nfree.father = th;
        
        for j in range(8):
            nfree.suns[j] = 0;

        nfree.sibling = 0;

        nfree.len = th.len/2;

        for j in range(3):
            nfree.center[j] = th.center[j];

        for j in range(3):
            if icoord[j] > nfree.center[j]:
                nfree.center[j] += nfree.len/2;
            else:
                nfree.center[j] -= nfree.len/2;

        nfree.partind = i;
        nfree.mass = imass;
        for j in range(3):
            nfree.s[j]  = (icoord[j] - nfree.center[j]);
            nfree.s[j] *= nfree.mass;

        th.suns[subi] = nfree;
        
        numnodestotal += 1
        nfree += 1

        if numnodestotal >= MaxNodes:
            print(stderr,"Maximum number %d of tree-nodes reached. i=%d, npart=%d. file: %s line: %d\n",
                    numnodestotal,i,npart,__FILE__,__LINE__);
            free(nodes);
            free(knlrad);
            free(knlpot);
            return False; 

    # now finish-up center-of-mass and quadrupol computation
    th = nodes
    for i in range(numnodestotal):
        for j in range(3):
            th.s[j] /= th.mass;
        
        if th.partind < 0:  # cell contains more than one particle
            dx = th.s[0];
            dy = th.s[1];
            dz = th.s[2];
            
            th.oc  = <float>sqrt(dx*dx + dy*dy + dz*dz);
            th.oc += th.len/(Thetamax); 
            th.oc *= th.oc;  # used in cell-opening criterion

        th.s[0] += th.center[0];
        th.s[1] += th.center[1];
        th.s[2] += th.center[2];
    
        # (j = 7, ; j >= 0; j--)
        for j in reversed(range(8)):  # preparations for non-recursive walk
            if th.suns[j]:
                th.suns[j].sibling = nn;
                nn = th.suns[j];
            
            if nn == 0:
                break

        th += 1

    last = 0;
    force_setupnonrecursive(nodes, &last);  # set up non-recursive walk
    last.next = 0;
    
    th = nodes
    for i in range(numnodestotal):
        if (~th.sibling):
            ff = th;
            nn = ff.sibling;

            while ~nn:
                ff = ff.father;
                if ~ff:
                    break;
                nn = ff.sibling;
        
            th.sibling = nn;
        th += 1

    for i in range(npart):
        pcoord[0] = x[i];
        pcoord[1] = y[i];
        pcoord[2] = z[i];

        Ep[i]  = force_treeevaluate_potential(nodes, pcoord, knlrad, knlpot);

        Ep[i] *= GCONS*Msol/Kpc;

        assert Ep[i] < 0.0;

    free(nodes);
    free(knlrad);
    free(knlpot);

    return True;
