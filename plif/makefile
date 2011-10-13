OTHERFILE=readme \
 makefile \
 test_fingerprint_classify_mocap.m 

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

PLIFFILE=condentropy.m \
 fingerprint.m \
 fingerprint_classify.m \
 fingerprint_compress.m \
 fingerprint_compress_dynamic.m
 

demo:
	# no demo now

zip: tar

tar: ${DYNAMMOFILE} ${PLIFFILE} ${OTHERFILE}
	mkdir plif
	cp ${DYNAMMOFILE} plif/
	cp ${PLIFFILE} plif/
	cp ${OTHERFILE} plif/
	zip -r plif.zip plif
	rm -r -f plif

sync: 
	svn up
	svn commit -m "commit" 
