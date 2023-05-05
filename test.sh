conda activate openduck

basedir=$(pwd)
testdir=$basedir/test

echo """OpenDUck testing

"""


## Protein Chunking
echo """Testing chunking

"""
cd ${testdir}/do_chunk
for i in {1..2}; do 
	cd test${i}
	echo "TEST ${i}"
	echo "openduck chunk -y test${i}.yaml"
	openduck chunk -y test${i}.yaml
	cd ..


## AMBER prearation
echo """Testing amber-prepare

"""
cd ${testdir}/amber
for i in {1..3}; do 
	cd test${i}
	echo "TEST ${i}"
	echo "openduck amber-prepare -y test${i}.yaml"
	openduck amber-prepare -y test${i}.yaml
	cd ..


## OPENMM
cd ${testdir}/openmm-prepare
echo """Testing openmm-prepare

"""
cd test1
echo "openduck openmm-prepare -y test1.yaml"
openduck openmm-prepare -y test1.yaml



echo """Testing openmm-full-protocol

"""
cd ${testdir}/openmm-full-protocol
for i in {1..2}; do 
	cd test${i}
	echo "TEST ${i}"
	echo "openduck openmm-full-protocol -y test${i}.yaml"
	openduck openmm-full-protocol -y test${i}.yaml
	cd ..



echo """Testing openmm-from-amber

"""
cd ${testdir}/openmm-from-amber
cd test1
echo "openduck openmm-from-amber -y test$1.yaml"
openduck openmm-from-amber -y test$1.yaml


echo """Testing report

"""
cd ${testdir}
openduck report -f openmm -p 'openmm-from-amber/test1' --plot
openduck report -f openmm -p 'openmm-full-protocol/test1' -d single

