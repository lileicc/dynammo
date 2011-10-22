OTHERFILE=EMSVD.m \
 readme \
 makefile \
 test_lds_fast.m \
 ../license.txt

DYNAMMOFILE=${OTHERFILE} \
 backward.m \
 compress_dynammo.m \
 compress_pca.m \
 cummin.m \
 decompress_dynammo.m \
 demo_dynammo.m \
 estimate_missing.m \
 forward.m \
 learn_lds.m \
 learn_lds_dynammo.m \
 linear_interp.m \
 logdet.m \
 MLE_lds.m  

DYNAMMOPFILE=${DYNAMMOFILE} \
 demo_dynammop.m \
 forward_fly.m \
 learn_lds_dynammop.m
 
REVNUM:=$(shell svn info |grep Revision: |cut -c11-)

VERSION:=$(shell echo "scale=1;(${REVNUM} - ${REVNUM} % 100) / 1000" | bc)


demo:
	matlab -r demo_dynammo
#	matlab -r demo_dynammop

zip: tar

tar: ${DYNAMMOFILE} ${DYNAMMOPFILE} 
	mkdir -p dynammo
	cp ${DYNAMMOFILE} dynammo/
	mkdir -p dynammop
	cp ${DYNAMMOPFILE} dynammop/
	zip dynammo-r${REVNUM}.zip dynammo
	zip dynammop-r${REVNUM}.zip dynammop
	rm -r -f dynammo
	rm -r -f dynammop

sync: 
	svn up
	svn commit -m "commit" 
	svn up
