// Modifico los include para poder wrappear Python
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <omp.h>
#include <Python.h>
#include <numpy/numpyconfig.h>
#include <numpy/arrayobject.h>

////#define PERIODIC // Las posiciones estan bien suponemos
#define Thetamax 0.45
#define KERNEL_LENGTH 10000
#define SOFT 3
#define GCONS 6.67300E-20    /* Constante de Gravitación [km³ / kg / seg²]     */
#define Msol 1.9891E30       /* Masa del Sol [Kg]                              */
#define Kpc 3.08568025E16    /* Kiloparsec -> Kilometro                        */  

struct NODE{
  float s[3];                     /* center of mass */
  long  partind;
  float center[3],len;            /* center and sidelength of treecubes */
  float mass;                     /* mass*/
  float oc;                       /* variable for opening criter*/
  struct NODE *next,*sibling,*father,*suns[8];
};

/************************************************************************/
/************************************************************************/

static void force_setkernel(float *knlrad, float *knlpot)
{
  int i;
  float u;
  float u2;
  float u3;
  float u4;
  float u5;

  for(i = 0; i <= KERNEL_LENGTH; i++){
    u = (float)i/(float)KERNEL_LENGTH;

    knlrad[i] = u;

    if(u <= 0.5){
      u2 = u*u;
      u4 = u2*u2;
      u5 = u4*u;
      knlpot[i]  = 16.0/3.0*u2;
      knlpot[i] -= 48.0/5.0*u4;
      knlpot[i] += 32.0/5.0*u5;
      knlpot[i] -= 14.0/5.0;
    }else{
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
    }
  }
}

/************************************************************************/
/************************************************************************/

static void add_particle_props_to_node(struct NODE *no, float *pos, float mp)
{
  int i;
  float delta;

  for( i = 0 ; i < 3 ; i++)
  {
    delta = (pos[i] - no->center[i]);
    no->s[i] += mp * delta;
  }

  no->mass += mp;

  return;
}

/************************************************************************/
/************************************************************************/

static void force_setupnonrecursive(struct NODE *no, struct NODE **last)
{
  int i;
  struct NODE *nn;
  
  if(*last)
    (*last)->next = no;

  *last = no;
  
  for(i = 0; i < 8; i++){
    nn = no->suns[i];
    if(nn != NULL)
      force_setupnonrecursive(nn, last);
  }
}

/****************************************/
/* packs the particles of group 'gr' into into BH-trees */
/****************************************/

static float force_treeevaluate_potential(const struct NODE *no, const float pos [], float *knlrad, float *knlpot)
{
  float r2,dx,dy,dz,r,u,h,ff;
  float wp;
  float h_inv;
  float local_pot;
  int ii;

  h = 2.8*(float)SOFT;
  h_inv = 1.0 / h;

  local_pot = 0.;

  while(no){
    dx = no->s[0] - pos[0];     
    dy = no->s[1] - pos[1];     
    dz = no->s[2] - pos[2];

    r2 = dx*dx + dy*dy + dz*dz;

    if(no->partind >= 0){   /* single particle */
      r = (float)sqrt(r2);  
     
      u = r * h_inv;

      if(u >= 1){
        local_pot -= no->mass/r;
      }else{
        ff = u*KERNEL_LENGTH;
        ii = (int)ff;
        ff -= knlrad[ii]*KERNEL_LENGTH;
        wp = knlpot[ii] + (knlpot[ii+1] - knlpot[ii])*ff;
        
        local_pot += no->mass*h_inv*wp;
      }
      no = no->sibling;
    }else{
      if(r2 < no->oc){
        no = no->next;  /* open cell */
      }else{
        r = (float)sqrt(r2);  
        u = r*h_inv;
    
        if(u >= 1){  /* ordinary quadrupol moment */
          local_pot += -no->mass/r;
        }else{    /* softened monopole moment */
          ff = u*KERNEL_LENGTH;
          ii = (int)ff;
          ff -= knlrad[ii]*KERNEL_LENGTH;
          wp = knlpot[ii] + (knlpot[ii+1]-knlpot[ii])*ff;

          local_pot += no->mass*h_inv*wp;
        }
        no = no->sibling;
      }
    }
  }

  return local_pot;
}

/****************************************************************/
static int tree_potential(const int npart, const float *mp, const float *x, const float *y, const float *z, float *Ep)
{
  int    i, j, subp,subi,p,subnode,fak;
  float  icoord[3], imass;
  float  pcoord[3], pmass;
  float  length;
  float  dx,dy,dz;
  float xmin[3],xmax[3];
  struct NODE *nodes, *last, *nfree,*th,*nn,*ff;
  int  MaxNodes, numnodestotal;
  float *knlrad, *knlpot; 
  knlrad = (float *) malloc((KERNEL_LENGTH+1) * sizeof(float));
  knlpot = (float *) malloc((KERNEL_LENGTH+1) * sizeof(float));

  for(i = 0; i < npart; i++)
    Ep[i] = 0.0;

  MaxNodes = 2 * npart + 200;
  nodes = (struct NODE *) malloc(MaxNodes*sizeof(struct NODE));
  assert(nodes!=NULL);

  force_setkernel(knlrad, knlpot);

  nfree = nodes;
  numnodestotal = 0;

  /* find enclosing rectangle */
  xmin[0] = xmax[0] = x[0];
  xmin[1] = xmax[1] = y[0];
  xmin[2] = xmax[2] = z[0];

  for(i = 1; i < npart; i++)
  {
    xmin[0] = x[i] < xmin[0] ? x[i] : xmin[0];
    xmin[1] = y[i] < xmin[1] ? y[i] : xmin[1];
    xmin[2] = z[i] < xmin[2] ? z[i] : xmin[2];
    xmax[0] = x[i] > xmax[0] ? x[i] : xmax[0];
    xmax[1] = y[i] > xmax[1] ? y[i] : xmax[1];
    xmax[2] = z[i] > xmax[2] ? z[i] : xmax[2];
  }
  
  for(j = 1 , length = xmax[0]-xmin[0] ; j < 3 ; j++)  /* determine maxmimum extension */
    if((xmax[j]-xmin[j]) > length)
      length = xmax[j]-xmin[j];
  length *= 1.01;
  assert(length > 0);

  /* insert first particle in root node */
  for(j = 0 ; j < 3 ; j++)
    nfree->center[j] = (xmax[j]+xmin[j])/2;
  nfree->len = length;
  
  /*inicializa variables*/
  nfree->father = 0;
  for(j = 0 ; j < 8 ; j++)
    nfree->suns[j] = 0;

  nfree->partind = 0;
  nfree->mass = mp[0];
  pcoord[0] = x[0];
  pcoord[1] = y[0];
  pcoord[2] = z[0];
  for(j = 0 ; j < 3 ; j++)
    nfree->s[j] = nfree->mass*(pcoord[j] - nfree->center[j]);
  
  /*inicializa la variable que apunta al hermano*/
  nfree->sibling = 0;

  numnodestotal++; nfree++;
  
  if(numnodestotal >= MaxNodes)
  {
    fprintf(stderr,"Maximum number %d of tree-nodes reached. file: %s line: %d\n",
            numnodestotal,__FILE__,__LINE__);
    free(nodes);
    free(knlrad);
    free(knlpot);
  }

  /* insert all other particles */
  for(i = 1 ; i < npart; i++)
  {
    th        = nodes;
    icoord[0] = x[i];
    icoord[1] = y[i];
    icoord[2] = z[i];
    imass     = mp[i];
 
    while(1)
    {
      add_particle_props_to_node(th,icoord,imass);

      if(th->partind >= 0)
        break;
    
      for(j = 0 , subnode = 0 , fak = 1 ; j < 3 ; j++ , fak <<= 1)
        if(icoord[j] > th->center[j])
          subnode += fak;

      nn = th->suns[subnode];
      if(nn != NULL)
        th = nn;
      else
        break;
    }
      
    if(th->partind >= 0)  /* cell is occcupied with one particle */
    {
      while(1)
      {
        p         = th->partind;
        pmass     = mp[p];
        pcoord[0] = x[p];
        pcoord[1] = y[p];
        pcoord[2] = z[p];

        for(j = 0 , subp = 0 , fak = 1 ; j < 3 ; j++ , fak <<= 1)
          if(pcoord[j] > th->center[j])
            subp += fak;

        nfree->father = th;
        
        for(j = 0 ; j < 8 ; j++)
          nfree->suns[j] = 0;

        nfree->sibling = 0;
        
        nfree->len = th->len/2;
    
        for(j = 0 ; j < 3 ; j++)
          nfree->center[j] = th->center[j];

        for(j = 0 ; j < 3 ; j++)
          if(pcoord[j] > nfree->center[j])
            nfree->center[j] += nfree->len/2;
          else
            nfree->center[j] -= nfree->len/2;

        nfree->partind = p;
        nfree->mass    = pmass;
        for(j = 0 ; j < 3 ; j++)
        {
          nfree->s[j]  = (pcoord[j] - nfree->center[j]);
          nfree->s[j] *= pmass;
        }
        th->partind = -1;
        th->suns[subp] = nfree;
     
        numnodestotal++; nfree++;
        if(numnodestotal >= MaxNodes)
        {
          fprintf(stderr,"Maximum number %d of tree-nodes reached. i=%d, npart=%d. file: %s line: %d\n",
                  numnodestotal,i,npart,__FILE__,__LINE__);
          free(nodes);
          free(knlrad);
          free(knlpot);
          return 0; 
        }

        for(j = 0 , subi = 0 , fak = 1 ; j < 3 ; j++ , fak <<= 1)
          if(icoord[j] > th->center[j])
            subi += fak;

        if(subi == subp)   /* the new particle lies in the same sub-cube */
        {
          th = nfree-1;
          add_particle_props_to_node(th,icoord,imass);
        }
        else
          break;
      }
    }
      
    for(j = 0 , subi = 0 , fak = 1 ; j < 3 ; j++ , fak <<= 1)
      if(icoord[j] > th->center[j])
        subi += fak;
      
    nfree->father = th;
      
    for(j = 0 ; j < 8 ; j++)
      nfree->suns[j] = 0;

    nfree->sibling = 0;

    nfree->len = th->len/2;

    for(j = 0 ; j < 3 ; j++)
      nfree->center[j] = th->center[j];

    for(j = 0 ; j < 3 ; j++)
      if(icoord[j] > nfree->center[j])
        nfree->center[j] += nfree->len/2;
      else
        nfree->center[j] -= nfree->len/2;

    nfree->partind = i;
    nfree->mass = imass;
    for(j = 0 ; j < 3 ; j++)
    {
      nfree->s[j]  = (icoord[j] - nfree->center[j]);

      nfree->s[j] *= nfree->mass;
    }
    th->suns[subi] = nfree;
      
    numnodestotal++; nfree++;

    if(numnodestotal >= MaxNodes)
    {
       fprintf(stderr,"Maximum number %d of tree-nodes reached. i=%d, npart=%d. file: %s line: %d\n",
              numnodestotal,i,npart,__FILE__,__LINE__);
       free(nodes);
       free(knlrad);
       free(knlpot);
       return 0; 
    }
  }

  /* now finish-up center-of-mass and quadrupol computation */
  
  for(i = 0, th = nodes; i < numnodestotal; i++, th++)
  {
    for(j = 0; j < 3; j++)
      th->s[j] /= th->mass;
      
    if(th->partind < 0)   /* cell contains more than one particle */
    {
      dx = th->s[0];
      dy = th->s[1];
      dz = th->s[2];
    
      th->oc  = (float)sqrt(dx*dx + dy*dy + dz*dz);
      th->oc += th->len/(Thetamax); 
      th->oc *= th->oc;     /* used in cell-opening criterion */
    }

    th->s[0] += th->center[0];
    th->s[1] += th->center[1];
    th->s[2] += th->center[2];
 
    for(j = 7, nn = 0; j >= 0; j--)    /* preparations for non-recursive walk */
    {
      if(th->suns[j])
      {
        th->suns[j]->sibling = nn;
        nn = th->suns[j];
      }
    }
  }

  last = 0;
  force_setupnonrecursive(nodes, &last);    /* set up non-recursive walk */
  last->next = 0;
  
  for(i = 0, th = nodes; i < numnodestotal; i++, th++)
  {
    if(!(th->sibling))
    {
      ff = th;
      nn = ff->sibling;

      while(!nn)
      {
        ff = ff->father;
        if(!ff)
        break;
        nn = ff->sibling;
      }
  
      th->sibling = nn;
    }
  }

  for(i = 0; i < npart; i++)
  {
    pcoord[0] = x[i];
    pcoord[1] = y[i];
    pcoord[2] = z[i];

    Ep[i]  = force_treeevaluate_potential(nodes, pcoord, knlrad, knlpot);

    Ep[i] *= GCONS*Msol/Kpc;

    assert(Ep[i]<0.0);
  }

  free(nodes);
  free(knlrad);
  free(knlpot);

  return 1;
}

// ###############################################
// Mucho problemas + parallelization omp jode la build -> lo commentareo
// ###############################################

/* static void force_brute(const int npart, const float *mp, const float *x, const float *y, const float *z, float *Ep)
{
  int  i, j;
  float pos[3], dx[3], dis;
  const int NTHREADS = omp_get_max_threads();


  for(i = 0; i < npart; i++)
    Ep[i] = 0.0;

  #pragma omp parallel for schedule(dynamic) num_threads(NTHREADS) \
  default(none) private(i,j,pos,dx,dis) \
  shared(npart,x,y,z,mp,Ep)

  for(i = 0; i < npart; i++)
  {
    pos[0] = x[i];
    pos[1] = y[i];
    pos[2] = z[i];

    // Calcula energia potencial
    for(j = 0; j < npart; j++)
    {

      if(i==j) continue;

      dx[0] = (pos[0] - x[j]);
      dx[1] = (pos[1] - y[j]);
      dx[2] = (pos[2] - z[j]);

      dis  = dx[0]*dx[0];
      dis += dx[1]*dx[1];
      dis += dx[2]*dx[2];
      dis  = sqrt(dis);

      Ep[i] += (mp[j])/dis;
    }

    Ep[i] *= (GCONS*Msol/Kpc);
    Ep[i] *= (-1.);                   // Cambio de signo para que Ep sea negativa
  }

  return;
} */

extern void calculate_potential(const int npart, const float *mp, const float *x, const float *y, const float *z, float *Ep)
{
  /* Si el halo tiene menos de 1000 particulas calcula la energia
  de forma directa (N^2). Si tiene mas de 1000 particulas usa 
  un octree. */
  // Mentira, no hagá nada para < 1000 parts
  if(npart > 1000)
  {
    if(!tree_potential(npart, mp, x, y, z, Ep))
    {
      fprintf(stderr,"Exit %d\n", npart);
      //force_brute(npart, mp, x, y, z, Ep);
    }
  
  }else{

    //force_brute(npart, mp, x, y, z, Ep);
  }

  return;
}


/* ----------------------------
Wrap para que lea Python, pseudo-copy
---------------------------- */

/* 1er intento de la func. Fallido...*/
// Func para pasar de np.arrays a C...
static float* PyToCArray(PyArrayObject *pyArr) {
    void* c_array_data = PyArray_DATA(pyArr);
    return (float*)c_array_data;
}


// En Python todo es un objeto, así que acá
// genero el método que es la función que me genera el 
// output deseado (la energía potencial)
static void pypot(const PyArrayObject *mp,
                         const PyArrayObject *x, const PyArrayObject *y,
                         const PyArrayObject *z, PyArrayObject *Ep){

    int npart = PyArray_NDIM(mp);

    float* c_mp = PyToCArray(mp);
    float* c_x = PyToCArray(x);
    float* c_y = PyToCArray(y);
    float* c_z = PyToCArray(z);
    float* c_Ep = PyToCArray(Ep);

    calculate_potential(npart, c_mp, c_x, c_y, c_z, c_Ep);

}

// 2do intento (a las ñapis):
/*
static PyObject*
pypot(PyObject* self, const PyObject* mp, const PyObject* x,
     const PyObject* y, const PyObject* z, PyObject* Ep)
{
    int npart = PyArray_NDIM(*mp);
    float* c_mp = PyArray_FromAny(*mp);
    float* c_x = PyArray_FromAny(*x);
    float* c_y = PyArray_FromAny(*y);
    float* c_z = PyArray_FromAny(*z);
    float* c_Ep = PyArray_FromAny(*Ep);
    
    PyObject *Ep = calculate_potential(npart, c_mp, c_x, c_y, c_z, c_Ep);
    
    return Ep;
}*/

// Defino los métodos
/* 1st try
static PyMethodDef ModuleMethods[] = {
  {"calculate_potential",
  method_calculate_potential,
  METH_VARARGS,
  "Python interface for Octree potential calculation in C"},
  {NULL, NULL, 0, NULL}
};
*/
// Genero el módulo
/* 1st try:
static struct PyModuleDef potential_module = {
    PyModuleDef_HEAD_INIT,
    "potential_module",
    "Interface for the Octree potential C function",
    -1,
    ModuleMethods};
*/
/*
static struct PyModuleDef potential_module = {
    PyModuleDef_HEAD_INIT,
    "potential_module",
    "Interface for the Octree potential C function",
    -1,
    ModuleMethods};



// Lo instancio
PyMODINIT_FUNC PyInit_calculate_potential_c(void) {
    import_array(); // Initialise Numpy
    if (PyErr_Occurred()) {
        return NULL;
    }
    return PyModule_Create(&module);
} */

// 2nd try (cambio "fib" por "pot"):
PyDOC_STRVAR(pot_doc, "computes the potential");
PyMethodDef pot_method = {
    "pot",                /* The name as a C string. */
    (PyCFunction) pypot,  /* The C function to invoke. */
    METH_O,               /* Flags telling Python how to invoke ``pyfib`` */
    pot_doc,              /* The docstring as a C string. */
};

PyMethodDef methods[] = {
    {"pot", (PyCFunction) pypot, METH_O, fpot_doc},
    {NULL},
};

PyDoc_STRVAR(pot_module_doc, "provides a potential function");
PyModuleDef pot_module = {
    PyModuleDef_HEAD_INIT,
    "pot",
    pot_module_doc,
    -1,
    methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit_pot(void)
{
    return PyModule_Create(&pot_module);
}