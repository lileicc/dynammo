OTHERFILE=readme \
 makefile \
 test_fingerprint_classify_mocap.m 

DYNAMMOFILE=backward.m \
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

PLIFFILE=condentropy.m \
 fingerprint.m \
 fingerprint_classify.m \
 fingerprint_compress.m \
 fingerprint_compress_dynamic.m
 

demo:
	# no demo now

zip: tar

tar: ${DYNAMMOFILE} ${PLIFFILE}
	zip plif.zip ${DYNAMMOFILE} ${PLIFFILE}

sync: 
	svn up
	svn commit -m "commit" 
