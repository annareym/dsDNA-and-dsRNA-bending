# dsDNA and dsRNA axial bending constraint

The dsDNA and dsRNA axial bending constraint can be used as any other [PLUMED Collective Variable](https://www.plumed.org/doc-v2.6/user-doc/html/colvarintro.html). The constraint is implemented within [PLUMED](https://www.plumed.org/) free energy library environment and can be used to perform all-atom [steered MD](https://en.wikipedia.org/wiki/Molecular_dynamics#Steered_molecular_dynamics_(SMD)) simulations of controlled bending for dsDNA and dsRNA molecules. 


## Author

Anna Reymer, [Reymer Lab](https://cmb.gu.se/english/about_us/staff?userId=xreyan), Department of Chemistry & Molecular Biology, University of Gothenburg, Sweden.


## Code dependencies (For installation instruction see below)

* [Eigen 3.3.7](http://eigen.tuxfamily.org/)
* [Autodiff 0.6.5](https://github.com/autodiff/autodiff)
* [Plumed 2.7.3](https://www.plumed.org/)


## How to use

The bending constraint works analogously to any collective variable (colvar) implemented in PLUMED. BEND colvar can monitor or control the value of the axial bend between a chosen set of base pairs, corresponding to the top _t_ and bottom _b_ in a dsDNA/dsRNA fragment. As an input to the colvar, namely plumed_bend.dat file, provide total 60 atoms' numbers: 8\*2\*3 atoms' numbers from DNA/RNA nucleotides that will be restrained during a simulation, numbers correspond to 8 – 4 restrained base levels at the top and 4 at the bottom, \* 2  DNA/RNA strands, \* 3 atoms per base (C1',N1/N9 and C6/C8 depending whether purines or pyrimidines, correspondingly), also provide additional 6\*2 auxiliary atoms evenly spaced  between the top and bottom atoms for each strand. The auxiliary atoms are important to avoid MD trajectory re-imaging issues during postprocessing. Additionally, in the input file provide a force constant (_KAPPA_) and desired value (_AT_) of the axial bend (in degrees). 
To push a system into the desired conformation, an energy penalty will be added to the potential energy functional: E<sub>Deform</sub>=0.5\*k\*(x<sub>0</sub>-x)<sup>2</sup>. If you want to just monitor the total twist or total stretch values while running MD, provide only 3\*2\*N atom numbers.

#### Example of plumed.dat input file for BEND:
```
bnd: BEND ATOMS=44,42,45,74,72,75,106,104,107,139,137,140,191,242,293,344,395,446,393,391,394,423,421,424,455,453,456,487,485,488,579,577,580,612,610,613,644,642,645,676,674,677,728,779,830,881,932,983,930,928,931,963,961,964,993,991,994,1025,1023,1026
bnd_r: RESTRAINT ARG=bnd KAPPA=0.3 AT=25
PRINT STRIDE=100 ARG=bnd,bnd_r.bias FILE=Bend_DNA_25deg.txt
```

Provided atom numbers should respect the following order. We start for Watson DNA/RNA strand (5'->3' backbone direction): first 3 atom numbers represent the base level _i_, which corresponds to the top of the restrained nucleic acids region, second 3 atom numbers represent the base level _i+1_, next 3 atom numbers – the base level _i+2_, then the base level _i+3_, followed by 6 auxiliary atom numbers, followed by the triads of atom numbers representing the base levels _j-3_, _j-2_, _j-1_, and _j_ on Watson strand. The _j_ base level is the bottom b.p. level in the restrained nucleic acid region. Then we switch to the atom numbers representing Crick DNA/RNA strand (3'->5' backbone direction) and provide the triads of atom numbers representing the base levels _j_, _j-1_, _j-2_, _j-3_, followed by 6 auxiliary atom numbers, finalised by triads of atom numbers corresponding to the base levels _i+3_, _i+2_, _i+1_, _i_. 


## Ho to install

The provided instructions (should be adjusted for any specific system) may not work for all systems. 
1. Download and install the Eigen and Autodiff software according to the respective instructions.
2. Download PLUMED. In PLUMED configure file edit on lines 4620-22 replace "-std=c++11" to "-std=c++17".
3. Set corresponding environmental variables (below are the paths that work for a correct system, change accordingly):

export LDFLAGS="-L/usr/local/opt/llvm/lib"
export CXXFLAGS="-std=c++17 -O3"
export CPPFLAGS="-I/usr/local/opt/llvm/include -I/Users/anna/eigen-3.3.7 -I/Users/anna/autodiff-0.6.5"

4. Follow instructions for PLUMED configuration and installation.


## Acknowledgements
This work was supported by Swedish Foundation for Strategic Research SSF Grant [ITM17-0431](https://strategiska.se/en/research/ongoing-research/instrument-technique-and-method-development-2017/project/9774/) to Dr. Anna Reymer.


## License

GNU LGPL.
