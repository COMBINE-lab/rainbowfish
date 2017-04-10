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
3. rankselect (https://github.com/efficient/rankselect.git)

How to install:
1. go through all the steps of installing VARI up until last step (building VARI)
2. git clone sparsepp and rankselect into 3rd_party_src
3. move sparsepp and rankselect hpp files into 3rd_party_inst/include
4. put all the files in rainbowfish repository into home directory of cosmo (including Makefile)
5. make (in cosmo directory)

