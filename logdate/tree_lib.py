import sys
#from dendropy import Tree
import logging
from treeswift import *

logger = logging.getLogger("tree_lib")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False

def prune_node(T,node):
    if node is not T.root:
        p = node.parent
        p.remove_child(node)
        if p.num_children() == 1:
            v = p.child_nodes()[0]
            p.remove_child(v)
            if p is T.root:
                T.root = v
            #    p.remove_child(v)
            else:
                u = p.parent
                l = p.edge_length + v.edge_length
                u.remove_child(p)
                u.add_child(v)
                v.edge_length = l



def prune_tree(T,RS):
# prune the taxa in the removing set RS from tree T
    L = list(T.traverse_leaves())
    for leaf in L:
        if leaf.label in RS:
            prune_node(T,leaf)

def get_taxa(tree_file,scheme='newick'):
    #a_tree = Tree.get_from_path(tree_file,scheme,preserve_underscores=True)
    a_tree = read_tree(tree_file,schema=scheme)
    return [leaf.label for leaf in a_tree.traverse_leaves()]

def report_taxa(tree_file,scheme='newick',listing=True,counting=True):
    #a_tree = Tree()
    #a_tree.read_from_path(tree_file,scheme)
    a_tree = read_tree(tree_file,schema=scheme)
    if listing:
        for leaf in a_tree.traverse_leaves():
            logger.info(leaf.label)
        if counting:
            logger.info('Taxa #: ' + str(len(list(a_tree.traverse_leaves()))))

'''
def tree_as_newick(a_tree,outfile=None,append=False):
# dendropy's method to write newick seems to have problem ...
	if outfile:
		outstream = open(outfile,'a') if append else open(outfile,'w')
	else:
		outstream = sys.stdout

	__write_newick(a_tree.root,outstream)

	outstream.write(";\n")
	if outfile:
		outstream.close()
'''

def tree_as_newick(a_tree,outfile=None,append=False):
    if outfile:
        outstream = open(outfile,'a') if append else open(outfile,'w')
    else:
        outstream = sys.stdout

    outstream.write(a_tree.newick() + "\n")

    if outfile:
        outstream.close()


'''
def __write_newick(node,outstream):
	if node.is_leaf():
		outstream.write(str(node.label))
	else:
		outstream.write('(')
		is_first_child = True
		for child in node.child_nodes():
			if is_first_child:
				is_first_child = False
			else:
				outstream.write(',')
			__write_newick(child,outstream)
		outstream.write(')')
	if not node.is_leaf() and node.label is not None:
			outstream.write(str(node.label))

	if not node.edge_length is None:
		outstream.write(":"+str(node.edge_length))
'''
