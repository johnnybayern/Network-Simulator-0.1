#-*-coding=utf-8-*-
from __future__ import division
__author__='Lu'
import csv
import snap
import random

def initalNet_from_csv(csvpath):
	edges=[]
	nodes=[]
	try:
		reader=csv.reader(open(csvpath,'rb'))
	except IOError:
		print 'no csv!'
	for line in reader:
		nodes+=[int(line[0]),int(line[1])]
		edges+=[(int(line[0]),int(line[1]))]
	g = snap.TUNGraph.New()
	for node in nodes:
		if not g.IsNode(node):
			g.AddNode(node)
		else:
			pass
	for edge in edges:
		g.AddEdge(edge[0],edge[1])
	return g
def getNodeInfoTable(g):
	nodeInfo=[]
	nodeInfo_dict={}
	for node in g.Nodes():
		neibs=[]
		for nth in range(node.GetDeg()):
			neib=node.GetNbrNId(nth)
			neibs+=[neib]
		nodeInfo+=[[node.GetId(),neibs,node.GetDeg(),0]]
		nodeInfo_dict[node.GetId()]=[neibs,node.GetDeg(),0]
	print  nodeInfo,nodeInfo_dict
	return nodeInfo,nodeInfo_dict
def getQofComm(comm,g):
	NIdV = snap.TIntV()
	for i in comm:
		NIdV.Add(i)
	SubG = snap.GetSubGraph(g, NIdV)
	l_in=SubG.GetEdges()
	l_out=0
	for n in comm:
		n_g=g.GetNI(n)
		for nth in range(0,n_g.GetDeg()):
			key=n_g.GetNbrNId(nth)
			# print str(key)+'key'
			if not key in comm:
				l_out+=1
			else:
				pass
	q=l_in/(l_in+l_out)
	return q
def getMaxQofN(nq_list):
	max_q=0
	n=0
	for v in nq_list:
		if v[1]>max_q:
			max_q=v[1]
			n=v[0]
	print n,max_q
	return n,max_q
def getNbrsofC(comm,nodeInfo_dict,q_dict):
	nbrs=[]
	print 'comm',comm
	for i in comm:
		print i,nodeInfo_dict[i]

		nbrs+=nodeInfo_dict[i][0]
		# print nbrs
	nbrs_set=set(nbrs)

	for i in q_dict.keys():
		if q_dict[i]>1 and i in nbrs_set:
			nbrs_set.remove(i)
	for i in comm:
		if i in nbrs_set:
			nbrs_set.remove(i)
	node_remove=[]
	for i,v in nodeInfo_dict.items():
		for j in nbrs_set:
			if i==j and v[2]>0:
				node_remove+=[j]
	for i in node_remove:
		nbrs_set.remove(i)
	return nbrs_set
def getCurrentandRestNet(n_list,g):
	# print n_list,'list'
	lis_g=[]
	NIdV1 = snap.TIntV()
	for i in g.Nodes():
		lis_g+=[i.GetId()]
	lis_rest=[i for i in lis_g if i not in n_list]
	# print lis_rest,'rest'
	for i in lis_rest:
		NIdV1.Add(i)
	SubG1 = snap.GetSubGraph(g, NIdV1)
	NIdV2 = snap.TIntV()
	for i in n_list:
		# print i
		NIdV2.Add(i)
	SubG2 = snap.GetSubGraph(g, NIdV2)
	return SubG1,SubG2
def fastDetCommunity(nodeInfo,nodeInfo_dict,net):
	g=net
	count = 1
	Communities={}
	while True:
		c=[]
		s=[]
		nodeInfo=[]
		q_dict={}#存储节点出现的次数
		for n,v in nodeInfo_dict.items():
			if v[2]==0:
				nodeInfo+=[n]
		if g.GetEdges()==0:#社区没有边退出循环
			break
		else:
			pass
		if len(nodeInfo)==0:#判断是否所有节点都被赋予社区
			break
		else:
			i_choose = random.choice(nodeInfo)
			c = [i_choose]
			nodeInfo_dict[i_choose][2]=count
			# print 'c',c
			s=getNbrsofC(c,nodeInfo_dict,q_dict)
			if len(s)==0:
				del nodeInfo_dict[i_choose]
				continue
			# s=nodeInfo_dict[i_choose][0]		
			# q_dict={}#存储节点出现的次数
		# print 'nodeinfo',nodeInfo_dict
		while True:	
			nq_list=[]#存储节点id和它对应的q
			q_c=getQofComm(c,g)
			for nbr in s:
				# nq_list=[]
				# q_dict={}
				c_t=c[:]
				c_t.append(nbr)
				# print 'ct',c_t
				q_i=getQofComm(c_t,g)
				nq_list+=[[nbr,q_i]]
				if q_i < q_c:
					if q_dict.has_key(nbr):
						q_dict[nbr]+=1
					else:
						q_dict[nbr]=1
					# if q_dict[nbr]<2:
					# 	continue
					# else:
					# 	s.remove(nbr)
				else:
					pass
			n_t,q_t=getMaxQofN(nq_list)
			# print 'qt',q_t,'nt',n_t
			if q_t>q_c:
				c+=[n_t]
				
				q_c=q_t
				nodeInfo_dict[n_t][2]=count
				s=getNbrsofC(c,nodeInfo_dict,q_dict)#更新s
			else:
				g,comm=getCurrentandRestNet(c,g)
				Communities[count]=comm
				count+=1
				break
	return Communities
	# for c_n,c in Communities.items():
	# 	snap.SaveEdgeList(c, 'fastdetComm//fastCommunity_'+str(c_n)+'.txt', 'Save as tab-separated list of edges')

if __name__ == '__main__':
	csvpath='net.csv'
	net=initalNet_from_csv(csvpath)
	# net=snap.GenRndGnm(snap.PUNGraph, 100, 80)
	nodeInfo,nodeInfo_dict=getNodeInfoTable(net)
	fastDetCommunity(nodeInfo,nodeInfo_dict,net)