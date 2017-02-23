from collections import deque


BIOLOGICAL_PROCESS = 'GO:0008150'
MOLECULAR_FUNCTION = 'GO:0003674'
CELLULAR_COMPONENT = 'GO:0005575'
FUNC_DICT = {
    'cc': CELLULAR_COMPONENT,
    'mf': MOLECULAR_FUNCTION,
    'bp': BIOLOGICAL_PROCESS}
EXP_CODES = set(['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'TAS', 'IC'])


def get_ontology(filename='go.obo'):
    # Reading Ontology from OBO Formatted file
    ont = dict()
    obj = None
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line == '[Term]':
                if obj is not None:
                    ont[obj['id']] = obj
                obj = dict()
                obj['is_a'] = list()
                obj['is_obsolete'] = False
                continue
            elif line == '[Typedef]':
                obj = None
            else:
                if obj is None:
                    continue
                l = line.split(": ")
                if l[0] == 'id':
                    obj['id'] = l[1]
                elif l[0] == 'is_a':
                    obj['is_a'].append(l[1].split(' ! ')[0])
                elif l[0] == 'name':
                    obj['name'] = l[1]

                elif l[0] == 'is_obsolete' and l[1] == 'true':
                    obj['is_obsolete'] = True
    if obj is not None:
        ont[obj['id']] = obj
    for node_id in ont.keys():
        if ont[node_id]['is_obsolete']:
            del ont[node_id]
    for node_id, val in ont.iteritems():
        if 'children' not in val:
            val['children'] = set()
        for n_id in val['is_a']:
            if n_id in ont:
                if 'children' not in ont[n_id]:
                    ont[n_id]['children'] = set()
                ont[n_id]['children'].add(node_id)

    return ont


def get_anchestors(ont, node_id):
    anchestors = set()
    q = deque()
    q.append(node_id)
    while(len(q) > 0):
        n_id = q.popleft()
        anchestors.add(n_id)
        for parent_id in ont[n_id]['is_a']:
            if parent_id in ont:
                q.append(parent_id)
    return anchestors


def get_parents(ont, node_id):
    parents = set()
    for parent_id in ont[node_id]['is_a']:
        if parent_id in ont:
            parents.add(parent_id)
    return parents


def get_subset(ont, node_id):
    subset = set()
    q = deque()
    q.append(node_id)
    while len(q) > 0:
        n_id = q.popleft()
        subset.add(n_id)
        for ch_id in ont[n_id]['children']:
            q.append(ch_id)
    return subset


def get_specifics(ont, node_set):
    res_set = node_set.copy()
    for node_id in node_set:
        anchestors = get_anchestors(ont, node_id)
        anchestors.discard(node_id)
        for anch_id in anchestors:
            res_set.discard(anch_id)
    return res_set
