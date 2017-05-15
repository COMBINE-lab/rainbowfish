# rainbowfish

## A succinct data structure to store and query colored DBG.

**Note** : Rainbowfish builds atop the BOSS representation of the dBG, as implemented in [VARI](https://github.com/cosmo-team/cosmo/tree/VARI).  Currently, installing and using Rainbowfish requires first properly installing VARI and then
checking out the Rainbowfish repository as a subdirectory of that project.

**Installation Guide**

Dependencies:
1. [VARI](https://github.com/cosmo-team/cosmo/tree/VARI) and all its dependencies:
	1. KMC2
	2. sdsl-lite
	3. stxxl
	4. tclap
2. [sparsepp](https://github.com/greg7mdp/sparsepp.git)

How to install:
1. go through all the steps of installing VARI
2. git clone sparsepp into 3rd_party_inst/include
4. `make` (in cosmo directory)
5. git clone rainbowfish *inside the main cosmo directory*.
6. cd into the rainbowfish sub-directory and `make`.

How to run:
1. Run all the first four steps in VARI, building KMC2 files and sorting kmers. 
After that building color matrix using command "cosmo-build"

To compare rb results with VARI you need to run "pack-color" command from VARI too.

2. pack colors using the color file created by cosmo and save them in hard disc:

```
> rb-pack-color kmc2_list_1000.colors 1000 <Address to bitvector dir> <1pass/2pass>
```

3. validating our query of (color & edge) vs VARI's:

```
> rb-validate kmc2_list_1000.dbg kmc2_list_1000.colors.sd_vector <Address to bitvector dir> <validation-type> > rb-validate 
```

validation-type should be one of the following words : 
* compare (to compare results with cosmo), 
* query (to go over all pairs of edge/color sequentially), 
* cosmo-query (to go over all pairs of edge/color on cosmo data structure),
* random-query (to go over a random pair of edge/color)

4. find bubbles:

```
> rb-find-bubble kmc2_list_1000.dbg kmc2_list_1000.colors.sd_vector \<Address to bitvector dir\> > rb-bubbles
```
