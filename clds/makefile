CLDSFILE=backward.m \
 forward.m \
 learn_clds.m \
 logdet.m \
 MLE_clds.m \
 condentropy.m \
 demo.m \
 makefile \
 readme.txt

demo:
	matlab -r demo

zip: tar

tar: ${CLDSFILE}
	zip clds.zip ${CLDSFILE} 

sync:
	svn up
	svn commit -m "commit"
