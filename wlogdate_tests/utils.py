from fastroot.Tree_extend import MPR_Tree, OGR_Tree
from treeswift import *
from logdate.logD_lib import logDate_with_random_init,f_wLogDate,f_wLogDate_changeVar,setup_constraint_changeVar,get_init_from_dated_tree,logIt,scale_tree
from logdate.fixed_init_lib import preprocess


EPSILON_t = 1e-4
EPS = 0.02

def test_logIt(inTreeFile,startTreeFile,smplTimeFile,refFile):
    tree = read_tree_newick(inTreeFile)
    startTree = read_tree_newick(startTreeFile)
    smplTimes = {}
    with open(smplTimeFile,'r') as fin:
        for line in fin:
            sp,t = line.strip().split()
            smplTimes[sp] = float(t)

    t0 = compute_divergence_time(startTree,smplTimes)
    preprocess(tree,smplTimes)
    f_obj = f_wLogDate_changeVar
    x = get_init_from_dated_tree(tree,startTree)
    
    cons_eq, b = setup_constraint_changeVar(tree, smplTimes)
    x0 = x + [t0]
    z0 = [x_i * b_i for (x_i, b_i) in zip(x0[:-2], b)] + [x0[-2], x0[-1]]
    _, f, z = logIt(tree, f_obj, cons_eq, b, x0=z0,verbose=False)
    x_opt = [z_i / b_i for (z_i, b_i) in zip(z[:-2], b)] + [z[-2], z[-1]]
    
    with open(refFile,'r') as fin:
        D = {}
        for line in fin:
            lb,true,wld = line.strip().split()
            D[lb] = (float(true),float(wld))
    scale_tree(tree,x_opt)
    failures = 0
    max_violate = 0
    for node in tree.traverse_preorder():
        if not node.is_root():
            e = node.edge_length 
            t,r = D[node.label]
            violate = (abs(t-e)-abs(t-r))/abs(t-r)
            max_violate = max(max_violate,violate)
            failures += violate < EPS

    return failures, max_violate

def compute_divergence_time(tree,sampling_time,bw_time=False,as_date=False):
# compute the divergence of the tree
# must have at least one sampling time. Assumming the tree branches have been
# converted to time unit and are consistent with the given sampling_time
    calibrated = []
    for node in tree.traverse_postorder():
        node.time,node.mutation_rate = None,None
        lb = node.label
        if lb in sampling_time:
            node.time = sampling_time[lb]
            calibrated.append(node)

    stk = []
    # push to stk all the uncalibrated nodes that are linked to (i.e. is parent or child of) any node in the calibrated list
    for node in calibrated:
        p = node.parent
        if p is not None and p.time is None:
            stk.append(p)
        if not node.is_leaf():
            stk += [ c for c in node.child_nodes() if c.time is None ]
    
    # compute divergence time of the remaining nodes
    while stk:
        node = stk.pop()
        lb = node.label
        p = node.parent
        t = None
        if p is not None:
            if p.time is not None:
                t = p.time + node.edge_length
            else:
                stk.append(p)    
        for c in node.child_nodes():
            if c.time is not None:
                t1 = c.time - c.edge_length
                t = t1 if t is None else t
                if abs(t-t1) > EPSILON_t:
                    print(t,t1,abs(t-t1))
                #assert abs(t-t1) < EPSILON_t, "Inconsistent divergence time computed for node " + lb
            else:
                stk.append(c)
        node.time = t
    return tree.root.time    
