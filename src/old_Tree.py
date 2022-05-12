class Tree():
    dist = None
    n_partitions = 0
    n_cuts = 0
    he = []
    regions = []
    def __init__(self, vertex, neighbors, edges, h, middle_edge):
        self.vertex = vertex
        self.neighbors = neighbors
        self.h = h
        self.edges = {}
        if type(edges) is list:
            self.create_arch(edges)
        else:
            self.edges = edges
        self.middle_edge = middle_edge
        if self.middle_edge is None:
            middle_edge, best = None, 0
            for edge in self.edges.keys():
                middle = min(len(self.edges[edge][min(edge)]),len(self.edges[edge][max(edge)]))
                if middle > best:
                    middle_edge, best = edge, middle
            self.middle_edge = middle_edge

    def heterogeneity(self,vertex):
        return Tree.dist[vertex].loc[vertex].sum().sum() / (len(vertex) * 2)

    def h_after_cut(self,edges,e):
        return self.heterogeneity(edges[e][e[0]]), self.heterogeneity(edges[e][e[1]])

    def create_arch(self,T):
        for e in T:
            e = (min(e),max(e))
            self.edges[e] = {e[0]:[e[0]],e[1]:[e[1]]}
            self.collect(e,e[0],e[1])
            self.collect(e,e[1],e[0])
            self.spread(self.edges[e][e[0]],e[1],e[0])
            self.spread(self.edges[e][e[1]],e[0],e[1])

    def collect(self,e,To,From):
        for v in self.neighbors[To]:
            if v != From:
                new_e = (min(To,v),max(To,v))
                if new_e in self.edges:
                    self.edges[e][To] += self.edges[new_e][v]

    def spread(self,val,To,From):
        for v in self.neighbors[To]:
            if v != From:
                new_e = (min(To,v),max(To,v))
                if new_e in self.edges:
                    self.edges[new_e][To] += val
                    self.spread(val,v,To)

    def find_regions(self, limit_vertex, z_min, z_max, w, sc):
        Tree.regions += [self.vertex]
        Tree.he += [self.h]
        if len(self.vertex) > limit_vertex:
            split_sizes = range(z_min,z_max+1)
            if w < len(split_sizes):
                split_sizes = np.random.choice(split_sizes,w,False)
            for z in split_sizes:
                new_trees = self.tree_partitioning(z,limit_vertex,sc)
                Tree.n_partitions += 1
                for tree in new_trees:
                    tree.find_regions(limit_vertex, z_min, z_max, w, sc)

    def tree_partitioning(self,z,limit_vertex,sc):
        trees_vertex = [self.vertex.copy()]
        trees_edges = [self.edges.copy()]
        trees_middle_edge = [self.middle_edge]
        neighbors = self.neighbors.copy()
        trees_h = [self.h]
        for k in range(1,z):
            best_ho = 0
            best_e = None
            best_tree = 0
            best_h1 = 0
            best_h2 = 0
            for i in range(k):
                if len(trees_vertex[i]) > limit_vertex:
                    e,h1,h2 = self.find_best_cut(trees_edges[i],trees_middle_edge[i],neighbors,sc)
                    Tree.n_cuts += 1
                    new_ho = trees_h[i] - h1 - h2
                    if new_ho > best_ho:
                        best_ho = new_ho
                        best_e = e
                        best_tree = i
                        best_h1, best_h2 = h1, h2
            if best_e is None:
                new_trees = []
                for i in range(len(trees_vertex)):
                    tree = Tree(trees_vertex[i], neighbors, trees_edges[i], trees_h[i], trees_middle_edge[i])
                    new_trees += [tree]
                return new_trees
            new_trees_vertex = []
            new_trees_edges = []
            new_trees_h = []
            new_trees_middle_edge = []
            for i in range(k):
                if i != best_tree:
                    new_trees_vertex += [trees_vertex[i]]
                    new_trees_edges += [trees_edges[i]]
                    new_trees_h += [trees_h[i]]
                    new_trees_middle_edge += [trees_middle_edge[i]]
                else:
                    vertex1, vertex2 = trees_edges[i][best_e][best_e[0]], trees_edges[i][best_e][best_e[1]]
                    edges1, edges2 = {}, {}
                    neighbors[best_e[0]] = [v for v in neighbors[best_e[0]] if v!=best_e[1]]
                    neighbors[best_e[1]] = [v for v in neighbors[best_e[1]] if v!=best_e[0]]
                    del trees_edges[i][best_e]
                    middle_edge1, middle_edge2 = None, None
                    best_middle1, best_middle2 = 0, 0
                    for edge in trees_edges[i].keys():
                        if edge[0] in vertex1:
                            edges1[edge] = {edge[0]:[v for v in trees_edges[i][edge][edge[0]] if v in vertex1],
                                            edge[1]:[v for v in trees_edges[i][edge][edge[1]] if v in vertex1]}
                            middle = min(len(edges1[edge][edge[0]]),len(edges1[edge][edge[1]]))
                            if middle > best_middle1:
                                best_middle1 = middle
                                middle_edge1 = edge
                        else:
                            edges2[edge] = {edge[0]: [v for v in trees_edges[i][edge][edge[0]] if v in vertex2],
                                            edge[1]: [v for v in trees_edges[i][edge][edge[1]] if v in vertex2]}
                            middle = min(len(edges2[edge][edge[0]]), len(edges2[edge][edge[1]]))
                            if middle > best_middle2:
                                best_middle2 = middle
                                middle_edge2 = edge
                    new_trees_vertex += [vertex1,vertex2]
                    new_trees_edges += [edges1,edges2]
                    new_trees_h += [best_h1,best_h2]
                    new_trees_middle_edge += [middle_edge1,middle_edge2]
            trees_vertex = new_trees_vertex
            trees_edges = new_trees_edges
            trees_h = new_trees_h
            trees_middle_edge = new_trees_middle_edge
        new_trees = []
        for i in range(z):
            tree = Tree(trees_vertex[i],neighbors,trees_edges[i],trees_h[i],trees_middle_edge[i])
            new_trees += [tree]
        return new_trees

    def find_best_cut(self, edges, middle_edge, neighbors, sc):
        best_h1, best_h2 = self.h_after_cut(edges,middle_edge)
        best_cut = middle_edge
        best_h = best_h1 + best_h2
        calculated = [middle_edge]
        i_since_upgrade = 0
        L = {}
        to_expand = middle_edge
        while i_since_upgrade < sc and len(calculated) < len(edges):
            for v1 in to_expand:
                for v2 in neighbors[v1]:
                    e = (min(v1,v2),max(v1,v2))
                    if e not in calculated:
                        calculated += [e]
                        h1, h2 = self.h_after_cut(edges, e)
                        h = h1 + h2
                        L[e] = max(h1,h2)
                        if h < best_h:
                            best_h, best_h1, best_h2, best_cut = h, h1, h2, e
                            i_since_upgrade = 0
            to_expand = min(L,key=L.get)
            del L[to_expand]
            i_since_upgrade += 1
        return best_cut, best_h1, best_h2

    def regions_to_onehot(self):
        R = pd.DataFrame(0,index=self.vertex,columns=range(len(Tree.regions)))
        for i in range(len(Tree.regions)):
            for el in Tree.regions[i]:
                R.loc[el,i] = 1
        return R

    def clear(self):
        Tree.regions = []
        Tree.n_cuts = 0
        Tree.n_partitions = 0
        Tree.he = []