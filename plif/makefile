OTHERFILE=readme \
 makefile \
 test_fingerprint_classify_mocap.m \
 ../condentropy.m

DYNAMMOFILE=../dynammo/backward.m \
 ../dynammo/compress_dynammo.m \
 ../dynammo/compress_pca.m \
 ../dynammo/cummin.m \
 ../dynammo/decompress_dynammo.m \
 ../dynammo/demo_dynammo.m \
 ../dynammo/estimate_missing.m \
 ../dynammo/forward.m \
 ../dynammo/learn_lds.m \
 ../dynammo/learn_lds_dynammo.m \
 ../dynammo/linear_interp.m \
 ../dynammo/logdet.m \
 ../dynammo/MLE_lds.m  

PLIFFILE=demo.m \
 fingerprint.m \
 fingerprint_classify.m \
 fingerprint_compress.m \
 fingerprint_compress_dynamic.m

MATLABPATH:=.:../:../dynammo/

REVNUM=`svn info |grep Revision: |cut -c11-`

VERSION=`echo (${REVNUM} - ${REVNUM} % 100) / 1000 | bc`

demo:
	export MATLABPATH=${MATLABPATH}; matlab -r 'demo'

zip: tar

tar: ${DYNAMMOFILE} ${PLIFFILE} ${OTHERFILE}
	mkdir plif
	cp ${DYNAMMOFILE} plif/
	cp ${PLIFFILE} plif/
	cp ${OTHERFILE} plif/
	zip -r plif${VERSION}-r${REVNUM}.zip plif
	rm -r -f plif

sync: 
	svn up
	svn commit -m "commit"
	svn up 
