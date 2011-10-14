OTHERFILE= readme \
 makefile

DYNAMMOFILE=${OTHERFILE} \
 ../dynammo/backward.m \
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
 ../dynammo/MLE_lds.m \
 ../dynammo/forward_fly.m \
 ../dynammo/learn_lds_dynammop.m


BOLEROFILE=${DYNAMMOFILE} \
 constructBodySkeleton.m \
 draw_skel_highlight.m \
 play_mocap_skel_highlight.m \
 draw_skel_length.m \
 play_mocap_skel_length.m \
 draw_skel.m \
 play_mocap_skel.m \
 draw_skel_print.m \
 get_bones.m \
 run_c3d_occlusion_bone.m \
 learn_lds_dynammop_bone.m \
 skeleton.mat \
 learn_lds_dynammop_bone_newton_direct.m \
 test_learn_bone_simulated.m \
 learn_lds_dynammop_bone_newton.m \
 test_run_c3d_occlusion_bone.m \
 learn_lds_dynammop_soft_bone.m

 
REVNUM:=$(shell svn info |grep Revision: |cut -c11-)

VERSION:=$(shell echo "scale=1;(${REVNUM} - ${REVNUM} % 100) / 1000" | bc)


demo:
#	matlab -r demo_dynammo
#	matlab -r demo_dynammop
	echo "No demo ready yet, sorry"

zip: tar

tar: ${BOLEROFILE} 
	mkdir -p bolero
	cp ${BOLEROFILE} bolero/	
	zip -r bolero-r${REVNUM}.zip bolero
	rm -r -f bolero/

sync: 
	svn up
	svn commit -m "commit" 
	svn up
