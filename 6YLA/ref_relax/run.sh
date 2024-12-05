ROSETTA_BIN=/home/zengyun1/App/rosetta.binary.ubuntu.release-371/main/source/bin
ROSETTA_BIN=/home/zengyun1/App/rosetta.binary.ubuntu.release-371/main/source/bin 
ROSETTA_DATABASE=/home/zengyun1/App/rosetta.binary.ubuntu.release-371/main/database

mkdir ref_relax

$ROSETTA_BIN/relax.default.linuxgccrelease -s origin.pdb  -nstruct 10 @relax.xml -out:path:all ref_relax
