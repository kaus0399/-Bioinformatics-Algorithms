def kmer_compsition(text, k):
    kmers = {}
    for i in range(len(text) - k + 1):
        kmer = text[i: i + k]
        kmers[kmer] = True
    # we have a dictionary that doesn't have duplicates
    kmer_list = []
    for kmer in kmers:
        kmer_list.append(kmer)
    # for x in kmer_list:
    #     print(x)    
    return kmer_list


# print(kmer_compsition("enter text here", "enter number"))


def spell_assembly(kmers):
    k = len(kmers[0])
    base = kmers[0]
    for kmer in kmers[1:]:
        base = base + kmer[k - 1]
    return base

# kmers = ['ACCGA',
#  'CCGAA',
#  'CGAAG',
#  'GAAGC',
#  'AAGCT'
# ]
# print(spell_assembly(kmers))


def create_overlap_graph(kmers):
    # create matrix
    a = []
    k = len(kmers[0])
    for i in range(k):
        a.append([])
        for j in range(k):
            a[i].append(0)
    # shorter version. Try using this, what goes wrong?
    #a = [[0] * len(kmers)] * len(kmers)
    for i, p in enumerate(kmers):
        for j, q in enumerate(kmers):
            if p[1:] == q[:k - 1]:
                a[i][j] = 1
    print(a)
    for i, p in enumerate(kmers):
        for j, q in enumerate(kmers):
            if a[i][j] == 1:
                print(p, '->', q)

# kmers = [
#     'ATGCG',
#     'GCATG',
#     'CATGC',
#     'AGGCA',
#     'GGCAT',
# ]


# create_overlap_graph(kmers)                


def build_de_bruijn_graph(k, text):
    kmers = []
    for i in range(len(text) - k + 1):
        kmer = text[i: i + k]
        kmers.append(kmer)
    matrix = {}
    for kmer in kmers:
        node = kmer[:k - 1]
        if not node in matrix:
            matrix[node] = []
        matrix[node].append(kmer[1:])
    for p in matrix:
        s = p + ' -> '
        for i, q in enumerate(matrix[p]):
            s += q
            if i != len(matrix[p]) - 1:
                s += ','
        print(s)


# print(build_de_bruijn_graph(4, 'AAGATTCTCTAC'))    


def build_de_bruijn_graph_from_kmers(kmers):
    matrix = {}
    k = len(kmers[0])
    for kmer in kmers:
        node = kmer[:k - 1]
        if not node in matrix:
            matrix[node] = []
        matrix[node].append(kmer[1:])
    for p in matrix:
        s = p + ' -> '
        for i, q in enumerate(matrix[p]):
            s += q
            if i != len(matrix[p]) - 1:
                s += ','
        print(s)

kmers = [
    'GAGG',
    'CAGG',
    'GGGG',
    'GGGA',
    'CAGG',
    'AGGG',
    'GGAG',
]
build_de_bruijn_graph_from_kmers(kmers)
        


def find_cycle(graph, v):
    cycle = []
    while len(graph[v]) != 0:
        u = graph[v][0]
        cycle.append(u)
        graph[v].pop(0)
        v = u
    return cycle

def append_cycle(graph, cycle, node):
    path = []
    v = node
    while len(graph[v]) != 0:
        u = graph[v][0]
        path.append(u)
        graph[v].pop(0)
        if u in cycle:
            i = cycle.index(node)
            for p in path:
                cycle.insert(i, p)
                i += 1
            break
        v = u
    print(cycle)
    return cycle

def eulerian_cycle(graph):
    n = len(graph)    
    cycle = find_cycle(graph, 0)
    print(graph)
    print(cycle)
    while True:
        for i, node in enumerate(cycle):
            if len(graph[node]) != 0:
                cycle = append_cycle(graph, cycle, node)
                break
        if len(cycle) == n + 1:
            break
    print(cycle)



# g = [0] * 10
# g[0] = [3]
# g[1] = [0]
# g[2] = [1,6]
# g[3] = [2]
# g[4] = [2]
# g[5] = [4]
# g[6] = [5,8]
# g[7] = [9]
# g[8] = [7]
# g[9] = [6]
# print(g)

# print(eulerian_cycle(g))    