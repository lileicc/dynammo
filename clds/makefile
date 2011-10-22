CLDSFILE=backward.m \
 forward.m \
 learn_clds.m \
 logdet.m \
 MLE_clds.m \
 ../condentropy.m \
 demo.m \
 makefile \
 readme.txt \
 sample_clds.m \
 dtw.m \
 dtw_cluster.m \
 test_clds_synthetic.m \
 test_clds_mocap16.m \
 test_clds_mocap35.m \
 ../colorGray.m \
 ../license.txt

MATLABPATH="./:../"

REVNUM:=$(shell svn info |grep Revision: |cut -c11-)

VERSION:=$(shell echo "scale=1;(${REVNUM} - ${REVNUM} % 100) / 1000" | bc)


demo:
	export MATLABPATH=${MATLABPATH}; matlab -r demo

zip: tar

tar: ${CLDSFILE}
	mkdir -p clds
	cp ${CLDSFILE} clds/
	zip -r clds-r${REVNUM}.zip clds
	rm -r -f clds/

sync:
	svn up
	svn commit -m "commit"
