
DIRS = ann_knn kdtree_knn cqt_knn sqt_knn cqt_render_tree sqt_render_tree

all:
	for dir in $(DIRS); do cd $$dir; make; cd ..; done

clean:
	for dir in $(DIRS); do cd $$dir; make clean; cd ..; done



