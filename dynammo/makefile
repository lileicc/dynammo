OTHERFILE=EMSVD.m \
 readme \
 makefile \
 test_lds_fast.m

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
 

demo:
	matlab -r demo_dynammo
	matlab -r demo_dynammop

zip: tar

tar: ${DYNAMMOFILE} ${DYNAMMOPFILE}
	zip dynammo.zip ${DYNAMMOFILE}
	zip dynammop.zip ${DYNAMMOPFILE}

sync: 
	svn up
	svn commit -m "commit" 
