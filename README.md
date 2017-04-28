# rainbowfish

## A succinct data structure to store and query colored DBG.

**Installation Guide**

Dependencies:
1. VARI ( https://github.com/cosmo-team/cosmo/tree/VARI ) and all its dependencies:
	1. KMC
	2. sdsl-lite
	3. stxxl
	4. tclap
2. sparsepp (https://github.com/greg7mdp/sparsepp.git)
3. ~~rankselect (https://github.com/efficient/rankselect.git)~~

How to install:
1. go through all the steps of installing VARI up until last step (building VARI)
2. git clone sparsepp and rankselect into 3rd_party_src
3. move sparsepp and rankselect hpp files into 3rd_party_inst/include
~~4. put all the files in rainbowfish repository into home directory of cosmo (including Makefile)~~
4. `make` (in cosmo directory)
5. git clone rainbowfish *inside the main cosmo directory*.
6. cd into the rainbowfish sub-directory and `make`.

How to run on NEWTON:
First of all create a directory in which you want to store all bitvectors (label, rank, and equivalence table bitvectors; compressed or uncompressed). 

The current queries we can run on newton are as follows (here we run it for 1000 colors):
1. pack colors using the color file created by cosmo and save them in hard disc:
LD_LIBRARY_PATH=/home/rob/cosmo/3rd_party_inst/lib gdb --args ~/projects/cosmo/rb-pack-color kmc2_list_1000.colors 1000 \<Address to bitvector dir\>
2. validating our query of (color & edge) vs VARI's:
LD_LIBRARY_PATH=/home/rob/cosmo/3rd_party_inst/lib ~/projects/cosmo/rb-validate kmc2_list_1000.dbg kmc2_list_1000.colors.sd_vector \<Address to bitvector dir\> > rb-validate
3. find bubbles:
LD_LIBRARY_PATH=/home/rob/cosmo/3rd_party_inst/lib ~/projects/cosmo/rb-find-bubble kmc2_list_1000.dbg kmc2_list_1000.colors.sd_vector \<Address to bitvector dir\> > rb-bubbles

