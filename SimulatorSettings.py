__author__ = 'sjha1'

from  Tkinter import *
from constants import *
import random as rd
from Node import *
import snap
import networkx as nx
import numpy.linalg
import matplotlib.pyplot as plt
# from shortestpath import *
import shortest_path_sim as sps
import fastdetCommunity as fd
import os
from collections import defaultdict
import numpy
from scipy.cluster import hierarchy
from scipy.spatial import distance

# Nodes_List = []
COLOR_LIST = ["Red", "Green", "Blue","Pink","LavenderBlush","HotPink","DeepPink","DarkViolet",
"Purple"]
# COLOR_LIST = ["Red", "Green", "Blue","Pink"]
COUNTER = 1


class simulatorWidget:
    def __init__(self, master, pre_canvas, textWidget):
        top = self.top = Toplevel(master)

        self.initValue_Checkbox = 1
        self.textWidget = textWidget
        self.top.title("Simulator Settings")
        # self.top.geometry('%dx%d+%d+%d' % (300, 420, 10, 500))
        # self.top.resizable(False, False)
        self.top.resizable(TRUE, TRUE)
        self.top.configure(background=BACKGROUND)
        self.Nodes_List = []

        self.pre_canvas = pre_canvas
        self.canvas = pre_canvas.mainCanvas

        self.generator_frame = LabelFrame(top, text="Network Generator", background=BACKGROUND)
        self.generator_frame.pack(padx=10, pady=10, expand=TRUE, fill=BOTH)

        self.sna_frame = LabelFrame(top, text="Network Analysis", background=BACKGROUND)
        self.sna_frame.pack(padx=10, pady=10, expand=TRUE, fill=BOTH)

        # self.det_frame = LabelFrame(self.sna_frame, text="Community Detection", background=BACKGROUND)
        # self.det_frame.pack(padx=10, pady=10, expand=TRUE, fill=BOTH)

        self.lb_nodes = Label(self.generator_frame, text="Number of Nodes:", background=BACKGROUND)
        self.lb_nodes.grid(row=0, column=0, sticky=E, pady=10, padx=(10,0) )

        self.txt_nodes = Entry(self.generator_frame, font="Helvetica 14")
        self.txt_nodes.grid(row=0, column=1, sticky=E+W, pady=10, padx=(0,10))

        self.lbl_Network = Label(self.generator_frame, text="Network:", background=BACKGROUND)
        self.lbl_Network.grid(row=1, column=0, sticky=E, pady=10)

        self.var_network = StringVar(top)
        self.var_network.set("Select Network")

        self.network_option = OptionMenu(self.generator_frame, self.var_network, "GenStar", "GenRndGnm", "GenForestFire"
                                         , "GenFull")
        self.network_option.grid(row=1, column=1, sticky=E+W, pady=10, padx=(0,10))

        self.drawNodes = Button(self.generator_frame, text="Generate Network", command=self.draw, background=BACKGROUND)
        self.drawNodes.grid(row=2, column=1, sticky=E+W, pady=10, padx=(0,10))

        self.resetButton = Button(self.generator_frame, text="Reset", command=self.reset, background=BACKGROUND)
        self.resetButton.grid(row=2, column=0, padx=10, pady=10)

        # self.lbl_commDetect = Label(self.sna_frame, text="Community Detection", background=BACKGROUND)
        # self.lbl_commDetect.grid(row=3, column=0, sticky=E, pady=10)

        self.txt_commDetect = Entry(self.sna_frame, font="Helvetica 14")
        self.txt_commDetect.grid(row=3, column=0, sticky=W, padx=10)

        self.btn_commDetect = Button(self.sna_frame, text="Detect Community", background=BACKGROUND,
                                     command=self.commDetection)
        self.btn_commDetect.grid(row=3, column=1, sticky=E+W, padx=10)

        self.btn_fastcommDetect = Button(self.sna_frame, text="Fast Detect Community", background=BACKGROUND,
                                         command=self.fastcommDetection)
        self.btn_fastcommDetect.grid(row=4, column=1, sticky=E+W, padx=10)

        self.shortDButton = Button(self.sna_frame, text="Shortest Path", background=BACKGROUND,
                                   command=self.find_short_path)
        self.shortDButton.grid(row=7, column=1, sticky=E, padx=10)

        self.eigenValueButton = Button(self.sna_frame, text="EigenValue Histogram", background=BACKGROUND,
                                       command=self.draw_eigen_hist)
        self.eigenValueButton.grid(row=5, column=0, sticky=E+W, padx=10, pady=10)

        self.rcmButton = Button(self.sna_frame, text="rcm", background=BACKGROUND,
                                command=self.rcm_matrix)
        self.rcmButton.grid(row=5, column=1, sticky=E+W, padx=10, pady=10)

        self.btwnessButton = Button(self.sna_frame, text="Betweenness Centrality", background=BACKGROUND,
                                    command=self.get_betweenness_centrality)
        self.btwnessButton.grid(row=6, column=0, sticky=E+W, padx=10)

        self.blockModelButton = Button(self.sna_frame, text="Block Model", background=BACKGROUND,
                                       command=self.block_model)
        self.blockModelButton.grid(row=6, column=1, sticky=E+W, padx=10)

        self.degreeButton = Button(self.sna_frame, text="Degree Centrality", background=BACKGROUND,
                                   command=self.get_degree_centrality)
        self.degreeButton.grid(row=7, column=0, sticky=E+W, padx=10)

        self.degree_hist_Button = Button(self.sna_frame, text="Show Degree Histogram", background=BACKGROUND,
                                         command=self.degree_histogram)
        self.degree_hist_Button.grid(row=7, column=1, sticky=E+W, padx=10)

        self.closenessButton = Button(self.sna_frame, text="Closeness Centrality", background=BACKGROUND,
                                      command=self.get_closeness_centrality)
        self.closenessButton.grid(row=8, column=0, sticky=E+W, padx=10)

        self.propertiesButton = Button(self.sna_frame, text="Network Properties", background=BACKGROUND,
                                       command=self.get_graph_properties)
        self.propertiesButton.grid(row=8, column=1, sticky=E+W, padx=10)

        self.checkEdges = Checkbutton(self.sna_frame, text="Enable Edges", background=BACKGROUND
                                      , command=self.show_edges)  # variable = self.varEdges
        self.checkEdges.grid(row=9, column=1, sticky=E, padx=10)
        self.checkEdges.select()

    def show_edges(self):
        if self.initValue_Checkbox:
            # print self.initValue_Checkbox
            try:
                for Node in self.Nodes_List:
                    Node.show_edges_toggle(self.canvas, self.initValue_Checkbox)
            except:
                pass
            self.initValue_Checkbox = 0
        else:
            # print self.initValue_Checkbox
            try:
                for Node in self.Nodes_List:
                    Node.show_edges_toggle(self.canvas, self.initValue_Checkbox)
            except:
                pass
            self.initValue_Checkbox = 1

    def read_nodes(self):
        self.Nodes = self.txt_nodes.get()
        if self.Nodes <> "":
            self.lst_Nodes = self.Nodes.split(",")
        else:
            print "Please enter nodes"

    def create_nodes(self, Graph):
        color = "white"  # rd.choice(COLOR_LIST)
        for node in Graph.Nodes():
            polygon = rd.randint(1, len(self.pre_canvas.polygon_dict))
            direction = self.pre_canvas.polygon_dict[polygon]

            if direction[0][0] == direction[0][1]:
                direction[0][0] = direction[0][0] - 1
            if (direction[1][0] == direction[1][1]):
                direction[1][0] = direction[1][0] - 1

            tempVar = False
            i = 1
            while not tempVar:
                x = rd.randrange(direction[0][0], direction[0][1])
                y = rd.randrange(direction[1][0], direction[1][1])
                tempVar = self.point_inside_polygon(x, y, direction[2])
                i = i + 1
                if i > 25:
                    break

            # print "x ",x," y ",y,
            p = Node(node.GetId(), [x, y], color)
            p.draw(self.canvas)
            self.canvas.tag_bind(p.itemNo, '<ButtonPress-1>', self.__showAttriInfo)
            self.Nodes_List.append(p)

        for node in Graph.Nodes():
            follower = []
            for EI in Graph.Edges():
                if EI.GetSrcNId() == node.GetId():
                    if EI.GetSrcNId() <> EI.GetDstNId():
                        follower.append(self.Nodes_List[EI.GetDstNId()])
            nodeid = node.GetId()
            self.Nodes_List[nodeid].followers = follower
            self.Nodes_List[nodeid].draw_edges(self.canvas)
            # print nodeid,len(follower)

    def draw(self):
        self.read_nodes()
        self.Network = self.var_network.get()
        self.nPoints = sum(int(i) for i in self.lst_Nodes)

        if self.Network == "GenStar":
            writeCalculations(self.textWidget, "GenStar is the network with points :" + str(self.nPoints), False)
            # print "GenStar is the network with points ",self.nPoints
            self.graph = snap.GenStar(snap.PNGraph, self.nPoints, True)

        if self.Network == "GenRndGnm":
            # print "GenRndGnm is the network with points ",self.nPoints
            writeCalculations(self.textWidget, "GenRndGnm is the network with points :" + str(self.nPoints), False)
            self.graph = snap.GenRndGnm(snap.PNGraph, self.nPoints, self.nPoints)

        if self.Network == "GenForestFire":
            # print "GenForestFire is the network with points ",self.nPoints
            writeCalculations(self.textWidget, "GenForestFire is the network with points :" + str(self.nPoints), False)
            self.graph = snap.GenForestFire(self.nPoints, 0.5, 0.5)

        if self.Network == "GenFull":
            writeCalculations(self.textWidget, "GenFull is the network with points :" + str(self.nPoints), False)
            # print "GenFull is the network with points ",self.nPoints
            self.graph = snap.GenFull(snap.PNGraph, self.nPoints)

        if self.Network == "GenCircle":
            writeCalculations(self.textWidget, "GenCircle is the network with points :" + str(self.nPoints), False)
            # print "GenCircle is the network with points ",self.nPoints
            self.graph = snap.GenCircle(snap.PNGraph, self.nPoints, 10, 10)

        self.create_nodes(self.graph)

        # for node in Graph.Nodes():
        #     follower = []
        #     for EI in Graph.Edges():
        #         if EI.GetSrcNId() == node.GetId():
        #             if EI.GetSrcNId() <> EI.GetDstNId():
        #                 follower.append(EI.GetDstNId())
        #     print node.GetId(),follower
        #
        #
        # for node in self.Nodes_List:
        #     str1 = ""
        #     for foll in node.followers:
        #         str1 = str1 + "," + str(foll.id)
        #     print node.id,str1

    def fastcommDetection(self):
        self.reset_nodes_edges()
        for node in self.Nodes_List:
            id = node.itemNo
            self.canvas.itemconfig(id,fill="white")
        nodeInfo,nodeInfo_dict=fd.getNodeInfoTable(self.graph)
        communities=fd.fastDetCommunity(nodeInfo,nodeInfo_dict,self.graph)
        for c_n, c in communities.items():
            snap.SaveEdgeList(c, 'fastdetComm//fastCommunity_'+str(c_n)+'.txt', 'Save as tab-separated list of edges')

            nodelist = [node.GetId() for node in c.Nodes()]
            col = rd.choice(COLOR_LIST)
            for i in nodelist:
                self.canvas.itemconfig(self.Nodes_List[i].itemNo,fill=col)

    def commDetection(self):
        # print "Community Detection Function"
        self.reset_nodes_edges()
        for node in self.Nodes_List:
            id = node.itemNo
            self.canvas.itemconfig(id, fill="white")

        conn_degree = 4  # Change after Lu code
        threshold = self.txt_commDetect.get()
        try:
            threshold = int(threshold)
            if threshold > 0 and type(threshold) == int:
                self.detectCommunity(threshold, conn_degree)
                # g = {}
                # for node in self.Nodes_List:
                #     if len(node.followers) > 0:
                #         g[node.id]=[]
                #         #print node.id,node.followers
                #         for foll in node.followers:
                #             g[node.id].append(foll.id)
                #
                # #print g
                # graph = Graph(g,self.Nodes_List,self.canvas,threshold)
                # graph.find_community()
                # graph.change_color()
                # else:
                #    print "Enter positive threshold"
        except Exception as e:
            print e  # "Enter positive threshold"

    def point_inside_polygon(self, x, y, poly):

        n = len(poly) / 2
        inside = False

        p1x = poly[0]
        p1y = poly[1]
        # print p1x,p1y
        for i in range(0, n + 1, 1):
            p2x = poly[(i % n) * 2]
            p2y = poly[(i % n) * 2 + 1]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y

        return inside

    def find_short_path(self):
        # dict = Generate_Dictionary_Simulation(self.Nodes_List)
        # Generate_Graph(dict)
        self.sp_Graph = sps.Graph()
        for node in self.Nodes_List:
            self.sp_Graph.add_node(str(node.id))

        for fromNode in self.Nodes_List:
            for toNode in fromNode.followers:
                self.sp_Graph.add_edge(str(fromNode.id), str(toNode.id), 1)

        self.select_nodes = []

    def reset(self):
        for node in self.Nodes_List:
            self.canvas.delete(node.itemNo)
            for edges in node.lineItemNo:
                self.canvas.delete(edges[0])
        self.Nodes_List = []

    def __showAttriInfo(self, event):
        """
        Show attribute information of clicked unit
        """
        widget_id = event.widget.find_closest(event.x, event.y)

        # print "x ",self.canvas.gettags(widget_id)[0]," y ",self.canvas.gettags(widget_id)[1]
        if self.canvas.gettags(widget_id)[0] == "id":
            id = self.canvas.gettags(widget_id)[1]
            print "Node id :", id, [i.id for i in self.Nodes_List[int(id)].followers]
            self.canvas.itemconfig(self.Nodes_List[int(id)].itemNo, fill="red")
            self.select_nodes.append(str(id))
            if len(self.select_nodes) > 1 and len(self.select_nodes) == 2:

                if (self.select_nodes[0] == self.select_nodes[1]):
                    idd = int(self.select_nodes[0])
                    print self.Nodes_List[idd].lineItemNo
                    for edge in self.Nodes_List[idd].lineItemNo:
                        self.canvas.itemconfig(edge[0], state=NORMAL)

                else:
                    try:
                        steps, path = sps.shortest_path(self.sp_Graph, self.select_nodes[0], self.select_nodes[1])
                        print "Number of steps: ", steps, "Path :", "-->".join(path)
                        print "Find path between ", self.select_nodes[0], self.select_nodes[1]

                        for fromNode in range(len(path) - 1):
                            for toNode in self.Nodes_List[int(path[fromNode])].lineItemNo:
                                if toNode[1] == int(path[fromNode + 1]):
                                    self.canvas.itemconfig(toNode[0], state=NORMAL)

                        self.canvas.itemconfig(self.Nodes_List[int(path[0])].itemNo, fill="orange")
                        self.canvas.itemconfig(self.Nodes_List[int(path[-1])].itemNo, fill="orange")


                    except Exception as e:
                        print "Error: No path exits to", e
            elif len(self.select_nodes) > 2:
                del self.select_nodes[0]
                del self.select_nodes[0]
                self.reset_nodes_edges()
                id = self.select_nodes[0]
                self.canvas.itemconfig(self.Nodes_List[int(id)].itemNo, fill="red")

    def reset_nodes_edges(self):
        for eachNode in self.Nodes_List:
            self.canvas.itemconfig(eachNode.itemNo, fill="white")
            for eachEdge in eachNode.lineItemNo:
                self.canvas.itemconfig(eachEdge[0], state=HIDDEN)
        filepath='fastdetComm' 
        filelist=os.listdir(filepath)
        for f in filelist:
            ifile=os.path.join(filepath,f)
            if os.path.isfile(ifile):
                os.remove(ifile)

    def detectCommunity(self, comm_size, conn_degree):
        communities = {}  # Declaring communities dictionary
        k = 0  # some value k = 0
        g = self.graph  # Getting a random graph (Network) g from initalNet()

        while not g.Empty():  # checking till graph g is not empty
            comm = []  # Declaring the community variable comm list
            n1 = self.getMaxDegree(g)  # Getting max Degree from getMaxDegree(g) into n1
            comm += [n1]
            while True:
                max_id, max_degree = self.getMaxDegreetoComm(comm, g)
                if max_degree != 0:
                    if len(comm) < comm_size or max_degree > conn_degree:
                        comm += [max_id]
                    else:
                        break
                else:
                    break
            if g.GetEdges() == 0:
                break
            else:
                pass
            if g.GetNodes() < comm_size + len(comm):
                for n in g.GetNodes():
                    comm += [n.GetId()]
            else:
                pass
            commk_name = 'community' + str(k)
            g, commk = self.getCurrentandRestNet(comm, g)
            communities[commk_name] = commk
            k += 1

        # print "node list length",len(self.Nodes_List)
        f = open("results.txt", 'w')
        for c_n, c in communities.items():
            nodelist = [node.GetId() for node in c.Nodes()]
            if len(nodelist) == comm_size:
                col = rd.choice(COLOR_LIST)
                for i in nodelist:
                    self.canvas.itemconfig(self.Nodes_List[i].itemNo, fill=col)

            string = ",".join(
                    str(n) + "(" + str(self.Nodes_List[n].position[0]) + "," + str(self.Nodes_List[n].position[1]) + ")"
                    for n in nodelist)
            # print string
            string = str(c_n)[:9] + ":" + str(c_n)[
                                          9:] + " " + string + "\n"  # " : " +",".join(str(n) for n in nodelist) + "\n"
            f.write(string)
            # print 'community: ' + c_n + ' ', nodelist

        f.close()

    def getMaxDegree(self, g):
        n = 0
        idi = 0

        for i in g.Nodes():
            if i.GetDeg() > n or i.GetDeg() == n:
                n = i.GetDeg()
                idi = i.GetId()
            else:
                pass

        return idi

    def getMaxDegreetoComm(self, comm, g):
        n_dgree = {}
        for n in comm:
            n_g = g.GetNI(n)
            for nth in range(0, n_g.GetDeg()):
                key = n_g.GetNbrNId(nth)
                # print str(key)+'key'
                if key not in comm:
                    if n_dgree.has_key(key):
                        n_dgree[key] += 1
                    else:
                        n_dgree[key] = 1
                else:
                    pass
        max_id = 0
        max_degree = 0
        for k, v in n_dgree.items():
            if v > max_degree:
                max_id = k
                max_degree = v
            else:
                pass
        return max_id, max_degree

    def getCurrentandRestNet(self, n_list, g):
        lis_g = []
        NIdV1 = snap.TIntV()
        for i in g.Nodes():
            lis_g += [i.GetId()]
        lis_rest = [i for i in lis_g if i not in n_list]
        for i in lis_rest:
            NIdV1.Add(i)
        SubG1 = snap.GetSubGraph(g, NIdV1)
        NIdV2 = snap.TIntV()
        for i in n_list:
            # print i
            NIdV2.Add(i)
        SubG2 = snap.GetSubGraph(g, NIdV2)
        return SubG1, SubG2

    def draw_eigen_hist(self):
        networkx_g = nx.Graph()
        g = self.graph
        for EI in g.Edges():
            networkx_g.add_edge(EI.GetSrcNId(), EI.GetDstNId())
            # print "edge is added", EI.GetSrcNId(), EI.GetDstNId()
        L = nx.normalized_laplacian_matrix(networkx_g)
        e = numpy.linalg.eigvals(L.A)
        print_seperator(self.textWidget, "EigenValue Histogram", False)
        # print("Largest eigenvalue:", max(e))
        writeCalculations(self.textWidget, "Largest eigenvalue: " + str(max(e)), False)
        writeCalculations(self.textWidget, "Smallest eigenvalue: " + str(min(e)), False)
        # print("Smallest eigenvalue:", min(e))
        # nx.draw(networkx_g)
        plt.hist(e, bins=100)  # histogram with 100 bins
        plt.xlim(0, 2)  # eigenvalues between 0 and 2
        plt.show()

    def get_betweenness_centrality(self):
        networkx_g = nx.Graph()
        g = self.graph
        for EI in g.Edges():
            networkx_g.add_edge(EI.GetSrcNId(), EI.GetDstNId())
        b = nx.betweenness_centrality(networkx_g)
        print_seperator(self.textWidget, "Betweenness Centrality", False)
        for v in networkx_g.nodes():
            writeCalculations(self.textWidget, str(v) + ": " + str(b[v]), False)

    def get_degree_centrality(self):
        networkx_g = nx.Graph()
        g = self.graph
        for EI in g.Edges():
            networkx_g.add_edge(EI.GetSrcNId(), EI.GetDstNId())
        b = nx.degree_centrality(networkx_g)
        print_seperator(self.textWidget, "Degree Centrality", False)
        for v in networkx_g.nodes():
            writeCalculations(self.textWidget, str(v) + ": " + str(b[v]), False)

    def get_closeness_centrality(self):
        networkx_g = nx.Graph()
        g = self.graph
        for EI in g.Edges():
            networkx_g.add_edge(EI.GetSrcNId(), EI.GetDstNId())
        b = nx.closeness_centrality(networkx_g)
        print_seperator(self.textWidget, "Closeness Centrality", False)
        for v in networkx_g.nodes():
            writeCalculations(self.textWidget, str(v) + ": " + str(b[v]), False)

    def degree_histogram(self):
        networkx_g = nx.Graph()
        g = self.graph
        for EI in g.Edges():
            networkx_g.add_edge(EI.GetSrcNId(), EI.GetDstNId())
        degree_sequence = sorted(nx.degree(networkx_g).values(),reverse=True) # degree sequence
        #print "Degree sequence", degree_sequence
        dmax = max(degree_sequence)

        plt.loglog(degree_sequence,'b-',marker='o')
        plt.title("Degree rank plot")
        plt.ylabel("degree")
        plt.xlabel("rank")

        # draw graph in inset
        plt.axes([0.45,0.45,0.45,0.45])
        Gcc = sorted(nx.connected_component_subgraphs(networkx_g), key=len, reverse=True)[0]
        pos = nx.spring_layout(Gcc)
        plt.axis('off')
        nx.draw_networkx_nodes(Gcc,pos,node_size=20)
        nx.draw_networkx_edges(Gcc,pos,alpha=0.4)

        # plt.savefig("degree_histogram.png")
        plt.show()

    def block_model(self):
        networkx_g = nx.Graph()
        g = self.graph
        for EI in g.Edges():
            networkx_g.add_edge(EI.GetSrcNId(), EI.GetDstNId())
        edge = nx.write_edgelist(networkx_g,"test.edgelist")
        edlst = nx.read_edgelist("test.edgelist")
        # Extract largest connected component into graph H
        temp = nx.connected_component_subgraphs(edlst)
        next_component = next(temp, None)
        # H=next_temp[0]
        # Makes life easier to have consecutively labeled integer nodes
        H=nx.convert_node_labels_to_integers(next_component)
        # Create parititions with hierarchical clustering
        partitions = create_hc(H)
        # Build blockmodel graph
        BM = nx.blockmodel(H,partitions)


        # Draw original graph
        pos=nx.spring_layout(H,iterations=100)
        fig=plt.figure(1,figsize=(6,10))
        ax=fig.add_subplot(211)
        nx.draw(H,pos,with_labels=False,node_size=10)
        plt.xlim(0,1)
        plt.ylim(0,1)

        # Draw block model with weighted edges and nodes sized by number of internal nodes
        node_size=[BM.node[x]['nnodes']*10 for x in BM.nodes()]
        edge_width=[(2*d['weight']) for (u,v,d) in BM.edges(data=True)]
        # Set positions to mean of positions of internal nodes from original graph
        posBM={}
        for n in BM:
            xy=numpy.array([pos[u] for u in BM.node[n]['graph']])
            posBM[n]=xy.mean(axis=0)
        ax=fig.add_subplot(212)
        nx.draw(BM,posBM,node_size=node_size,width=edge_width,with_labels=False)
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.axis('off')
        # plt.savefig('hartford_drug_block_model.png')
        plt.show()

    def rcm_matrix(self):
        # Cuthill-McKee ordering of matrices
        # The reverse Cuthill-McKee algorithm gives a sparse matrix ordering that
        # reduces the matrix bandwidth.
        # Requires NumPy
        # Copyright (C) 2011 by
        # Aric Hagberg <aric.hagberg@gmail.com>
        # BSD License
        networkx_g = nx.Graph()
        g = self.graph
        for EI in g.Edges():
            networkx_g.add_edge(EI.GetSrcNId(), EI.GetDstNId())
        rcm = list(nx.utils.reverse_cuthill_mckee_ordering(networkx_g))
        print_seperator(self.textWidget, "Cuthill-McKee ordering of matrices", False)
        writeCalculations(self.textWidget, "ordering->" + str(rcm), False)

        writeCalculations(self.textWidget, "unordered Laplacian matrix", False)
        A = nx.laplacian_matrix(networkx_g)
        x,y = numpy.nonzero(A)
        #print("lower bandwidth:",(y-x).max())
        #print("upper bandwidth:",(x-y).max())
        bandwith1 = (y-x).max()+(x-y).max()+1
        writeCalculations(self.textWidget, "bandwidth: " + str(bandwith1), False)
        writeCalculations(self.textWidget, str(A), False)

        B = nx.laplacian_matrix(networkx_g,nodelist=rcm)
        writeCalculations(self.textWidget, "low-bandwidth Laplacian matrix", False)
        x,y = numpy.nonzero(B)
        #print("lower bandwidth:",(y-x).max())
        #print("upper bandwidth:",(x-y).max())
        bandwith2 = (y-x).max()+(x-y).max()+1
        writeCalculations(self.textWidget, "bandwidth: " + str(bandwith2), False)
        writeCalculations(self.textWidget, str(B), False)

    def get_graph_properties(self):
        networkx_g = nx.Graph()
        g = self.graph
        for EI in g.Edges():
            networkx_g.add_edge(EI.GetSrcNId(), EI.GetDstNId())

        pathlengths=[]

        print_seperator(self.textWidget, "Network Properties", False)
        writeCalculations(self.textWidget, "source vertex {target:length, }", False)
        for v in networkx_g.nodes():
            spl=nx.single_source_shortest_path_length(networkx_g,v)
            # writeCalculations(self.textWidget, str(v)+" "+str(spl), False)
            for p in spl.values():
                pathlengths.append(p)

        writeCalculations(self.textWidget, "", False)

        writeCalculations(self.textWidget, "average shortest path length " +
                          str(sum(pathlengths)/len(pathlengths)), False)
        # histogram of path lengths
        dist = {}
        for p in pathlengths:
            if p in dist:
                dist[p] += 1
            else:
                dist[p] = 1

        writeCalculations(self.textWidget, "", False)
        writeCalculations(self.textWidget, "length #paths", False)

        verts = dist.keys()
        for d in sorted(verts):
            writeCalculations(self.textWidget, str(d) + " " + str(dist[d]), False)

        try:
            writeCalculations(self.textWidget, "radius: " + str(nx.radius(networkx_g)), False)
            writeCalculations(self.textWidget, "diameter: " + str(nx.diameter(networkx_g)), False)
            writeCalculations(self.textWidget, "eccentricity: " + str(nx.eccentricity(networkx_g)), False)
            writeCalculations(self.textWidget, "center: " + str(nx.center(networkx_g)), False)
            writeCalculations(self.textWidget, "periphery: " + str(nx.periphery(networkx_g)), False)
            writeCalculations(self.textWidget, "density: " + str(nx.density(networkx_g)), False)
        except Exception as e:
            writeCalculations(self.textWidget, str(e) + "\nPlease use finite networks", TRUE)


class Graph:
    def __init__(self, g, Node_List, canvas, threshold):
        self.g = g
        self.total_community = []
        self.weights = []
        self.keys = []
        self.Node_List = Node_List
        self.canvas = canvas
        self.threshold = threshold

    def max_weight(self):
        self.weights = []
        self.keys = []
        for node, connections in self.g.items():
            self.keys.append(node)
            self.weights.append(len(connections))
        try:
            # print len(self.keys),self.keys
            maximum_weight = max(self.weights)
            temp = self.weights.index(maximum_weight)

            return self.keys[temp]
        except:
            return None

    def find_community(self):
        while len(self.g.keys()) <> 0:
            node = self.max_weight()
            self.delete(node)
            self.c = Community(node, self.threshold)
            while self.c.members(self.g):
                node = self.max_weight()
                self.c.add_descendants(node)
                self.delete(node)
            self.total_community.append(self.c)
            print self.c
        print "Completed finding community"

    def change_color(self):

        for node in self.Node_List:
            item = node.itemNo
            self.canvas.itemconfig(item, fill="white")

        for community in self.total_community:
            r = lambda: rd.randint(0, 255)
            a = ('#%02X%02X%02X' % (r(), r(), r()))

            if None not in community.descendants:
                parent = community.parent

                i = self.Node_List[parent].itemNo
                self.canvas.delete(i)
                self.Node_List[parent].drawRectangle(self.canvas)  # ,6*len(community.descendants)

                i = self.Node_List[parent].itemNo
                self.canvas.itemconfig(i, fill=a)

                for foll in community.descendants:
                    # print parent,foll
                    i = self.Node_List[foll].itemNo
                    self.canvas.itemconfig(i, fill=a)

            else:
                pass
                # print "sagar ", community

    def delete(self, node):
        try:
            self.g.pop(node)
        except:
            pass
            # print "Cannot delete None Key"

    def followers(self, node):
        return self.g[node]


class Community:
    def __init__(self, parent, threshold):
        self.parent = parent
        self.descendants = []
        self.threshold = threshold

    def __str__(self):
        string = "Parent: " + str(self.parent) + " Decendents: " + str(self.descendants)
        return string

    def add_descendants(self, node):
        self.descendants.append(node)

    def members(self, g):
        global threshold
        if g.keys() <> 0 and len(self.descendants) < self.threshold:
            return True
        return False


def maindraw():
    pass

# for block model
def create_hc(H):
    """Creates hierarchical cluster of graph G from distance matrix"""
    path_length=nx.all_pairs_shortest_path_length(H)
    distances=numpy.zeros((len(H),len(H)))
    for u,p in path_length.items():
        for v,d in p.items():
            distances[u][v]=d
    # Create hierarchical cluster
    Y=distance.squareform(distances)
    Z=hierarchy.complete(Y)  # Creates HC using farthest point linkage
    # This partition selection is arbitrary, for illustrive purposes
    membership=list(hierarchy.fcluster(Z,t=1.15))
    # Create collection of lists for blockmodel
    partition=defaultdict(list)
    for n,p in zip(list(range(len(H))),membership):
        partition[p].append(n)
    return list(partition.values())
