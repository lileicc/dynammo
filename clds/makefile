CLDSFILE=backward.m \
 forward.m \
 learn_clds.m \
 logdet.m \
 MLE_clds.m \
 condentropy.m \
 demo.m \
 makefile \
 readme.txt \
 sample_clds.m \
 dtw.m \
 dtw_cluster.m \
 test_clds_synthetic.m \
 test_clds_mocap16.m \
 test_clds_mocap35.m \
 colorGray.m 

demo:
	matlab -r demo

zip: tar

tar: ${CLDSFILE}
	zip clds.zip ${CLDSFILE} 

sync:
	svn up
	svn commit -m "commit"
