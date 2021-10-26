import csv
from numpy import array,zeros,linspace
from nested import*
from igraph import Graph,plot
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np


def net_2_mat(net):
	p = sorted(list(set([i[0] for i in net])))
	a = sorted(list(set([i[1] for i in net])))
	mat = zeros([len(p),len(a)])
	for i in net:
		r = p.index(i[0])
		c = a.index(i[1])
		mat[r][c] = 1
	return mat



def nodf_an(mat):
	nodf = NODF(mat)
	null = array([NODF(EE(mat.copy())) for i in range(1000)])
	p_value = (null>=nodf).sum()/1000.0
	z = (nodf-mean(null))/std(null)
	return nodf,mean(null),p_value,z


from colormaps import plasma
def plot_net(w_mat,fname):
	net = []
	r,c = w_mat.shape
	for i in range(r):
		for j in range(c):
			if w_mat[i][j]>0:
				net.append(['p_'+str(i),'a_'+str(j),w_mat[i][j]])
	g = Graph.TupleList(net,directed=True,weights='w')
	g.vs['type'] = [1 if 'p' in i else 0 for i in g.vs['name']]
	deg = g.degree()
	deg0 = [deg[i] for i in range(len(deg)) if g.vs['type'][i]==0]
	deg1 = [deg[i] for i in range(len(deg)) if g.vs['type'][i]==1]
	deg0 = [i/float(max(deg0)) for i in deg0]
	deg1 = [i/float(max(deg1)) for i in deg1]
	g.vs['size'] = 12
	sc0,sc1 = 0, 0
	for i in range(len(g.vs)):
		if g.vs['type'][i] == 0:
			g.vs[i]['size']*=deg0[sc0]
			sc0+=1
		if g.vs['type'][i] == 1:
			g.vs[i]['size']*=deg1[sc1]
			sc1+=1
		g.vs[i]['size']+=5
	lay = g.layout_bipartite()
	n0 = len([i for i in g.vs['type'] if i==0])
	n1 = len([i for i in g.vs['type'] if i==1])
	x0 = linspace(0,max(n0,n1),n0)
	x1 = linspace(0,max(n0,n1),n1)
	sc0,sc1 = 0,0
	for i in range(len(lay)):
		if lay[i][1] == 1: ##note that the type 0 nodes have layout_y = 1
			lay[i][0] = x0[sc0]
			lay[i][1] = 0
			sc0+=1
		else:
			lay[i][0] = x1[sc1]
			lay[i][1] = 1
			sc1+=1
	g.vs['color'] = ['darkgreen' if 'p' in i else 'blue' for i in g.vs['name']]
	g.es['width'] = [1+(i*2) for i in g.es['w']]
	#std_w = [int(round(255*(i-min(g.es['w']))/float(max(g.es['w'])-min(g.es['w'])))) for i in g.es['w']]
	#g.es['color'] = [plasma(i) for i in std_w]
	g.es['arrow_size'] = 0.5
	visual_style = {}
	visual_style["layout"] = lay
	visual_style["bbox"] = (600, 300)
	visual_style["margin"] = 50
	visual_style["vertex_frame_color"] = '#FFFFFF'
	plot(g,fname,**visual_style)







####huang # weights by preference of interactions
h = [i for i in csv.reader(open('huang.csv','r'))]
h_mat = array([map(float,i[1:-1]) for i in h[1:]])
h_mat = h_mat.T #algae in row
h_mat = h_mat/h_mat.sum() #normalize by total interaction per animal

plot_net(h_mat,'huang1.pdf')


h_mat = array([map(float,i[1:-1]) for i in h[1:]])
h_mat = h_mat.T #algae in row
h_mat = h_mat/h_mat.sum(0) #normalize by total interaction per animal

plot_net(h_mat,'huang2.pdf')


h_mat = array([map(float,i[1:-1]) for i in h[1:]])
h_mat = h_mat.T #algae in row
h_mat = h_mat/h_mat.sum(0) #normalize by total interaction per animal
h_mat = 0.15*(h_mat>0.01)###used 0.5 for plotting reasons
plot_net(h_mat,'huang3.pdf')


h_mat = array([map(float,i[1:-1]) for i in h[1:]])
h_mat = h_mat.T #algae in row
h_mat = h_mat/h_mat.sum(0)
out = open('conn_vs_tre.csv','w')
out.write('tre,algae,amph,size,occ,conn\n')
for tre in arange(0,1.0,0.001):
	h_mat_ = 1*(h_mat>tre)
	#h_mat_ = delete(h_mat_,where(h_mat_.sum(0)==0),1)
	#h_mat_ = delete(h_mat_,where(h_mat_.sum(1)==0),0)
	occ = h_mat_.sum()
	size = float(h_mat_.size)
	algae,amph = sum(h_mat_.sum(1)>0),sum(h_mat_.sum(0)>0)
	conn = occ/size
	#nodf,null_nodf,p,z = nodf_an(h_mat_)
	out.write(','.join(map(str,[tre,algae,amph,size,occ,conn]))+'\n')
	print (tre)


out.close()




h_mat = array([map(float,i[1:-1]) for i in h[1:]])
h_mat = h_mat.T #algae in row
h_mat = h_mat/h_mat.sum()
out = open('conn_vs_tre_tot.csv','w')
out.write('tre,algae,amph,size,occ,conn\n')
for tre in arange(0,0.1,0.0001):
	h_mat_ = 1*(h_mat>tre)
	#h_mat_ = delete(h_mat_,where(h_mat_.sum(0)==0),1)
	#h_mat_ = delete(h_mat_,where(h_mat_.sum(1)==0),0)
	occ = h_mat_.sum()
	size = float(h_mat_.size)
	algae,amph = sum(h_mat_.sum(1)>0),sum(h_mat_.sum(0)>0)
	conn = occ/size
	#nodf,null_nodf,p,z = nodf_an(h_mat_)
	out.write(','.join(map(str,[tre,algae,amph,size,occ,conn]))+'\n')
	print (tre,conn)


out.close()


#####plot networks threshold
h_mat = array([map(float,i[1:-1]) for i in h[1:]])
h_mat = h_mat.T #algae in row
h_mat = h_mat/h_mat.sum(0) #normalize by total interaction per animal
h_mat = 0.5*(h_mat>0.295)###used 0.5 for plotting reasons
plot_net(h_mat,'huang_tre_amph.pdf')


h_mat = array([map(float,i[1:-1]) for i in h[1:]])
h_mat = h_mat.T #algae in row
h_mat = h_mat/h_mat.sum() #normalize by total interaction per animal
h_mat = 0.5*(h_mat>0.0014)###used 0.5 for plotting reasons
plot_net(h_mat,'huang_tre_algae.pdf')

#############nakamoto

p = [i for i in csv.reader(open('plants.csv','r'))]
a = [i for i in csv.reader(open('animals.csv','r'))]


yyy = set([i[0] for i in p[1:]])

nets = []
bm = []
for y in yyy:
	mat_p = array([map(float,i[2:]) for i in p[1:] if i[0]==y])
	p_names = p[0][2:]
	bm_dict = dict(zip(p_names,mat_p.sum(0)))
	bm.append(bm_dict)
	mat_p = 1*(mat_p>1) #treshold for plant biomass
	mat_p = mat_p.T #sites in column, species in rows
	mat_a = array([map(float,i[2:]) for i in a[1:] if i[0]==y])
	a_names = a[0][2:]
	mat_a = 1*(mat_a>4) #treshold for animal individuals
	mat_a = mat_a.T
	null_a = [EF(mat_a) for i in range(1000)]
	null_p = [EF(mat_p) for i in range(1000)]
	net = []
	for i in range(mat_a.shape[0]):
		for j in range(mat_p.shape[0]):
			obs = (mat_a[i]*mat_p[j]).sum()
			null = [(null_a[k][i]*null_p[k][j]).sum() for k in range(1000)]
			p_value = (null>=obs).sum()/1000.0
			net.append(['p_'+str(j),'a_'+str(i),p_value]) #plant,animal
	nets.append(net)


nets[0],nets[1] = [i for i in nets[0] if i[2]<0.01],[i for i in nets[1] if i[2]<0.01]
n_mat_1,n_mat_2 = net_2_mat(nets[0]),net_2_mat(nets[1])
out = open('nakamoto_mat_1.csv','w')
for i in n_mat_1:
	out.write(','.join(map(str,i))+'\n')

out.close()

out = open('nakamoto_mat_2.csv','w')
for i in n_mat_2:
	out.write(','.join(map(str,i))+'\n')

out.close()

g1 = Graph.TupleList(nets[0],directed=True)
g2 = Graph.TupleList(nets[1],directed=True)



g1.vs['type'] = [1 if 'p' in i else 0 for i in g1.vs['name']]
g2.vs['type'] = [1 if 'p' in i else 0 for i in g2.vs['name']]
spp_14 = set(g1.vs['name'])
spp_15 = set(g2.vs['name'])
both = spp_14&spp_15
net_2014 = set(['_'.join(map(str,i)) for i in nets[0] if (i[0] in both and i[1] in both)])
net_2015 = set(['_'.join(map(str,i)) for i in nets[1] if (i[0] in both and i[1] in both)])
len(net_2015&net_2014)/float(len(net_2014))

net_2014 = set(['_'.join(map(str,i)) for i in nets[0]])
net_2015 = set(['_'.join(map(str,i)) for i in nets[1]])
len(net_2015&net_2014)/float(len(net_2014))

#plot first individual and then both

glob_net = nets[0]+nets[1]
g = Graph.TupleList(glob_net,directed=True)
g.es['weight'] = 1
g.simplify(combine_edges="sum")

g.vs['type'] = [1 if 'p' in i else 0 for i in g.vs['name']]
g.vs['color'] = ['darkgreen' if 'p' in i else 'blue' for i in g.vs['name']]
g.es['color'] = ['lightgrey' if i==1 else 'black' for i in g.es['weight']]

g.es['width'] = 2
g.vs['size'] = 15#[20 if i in both else 10 for i in g.vs['name']]
col = []
for i in g.vs['name']:
	if i in both and 'p' in i:
		col.append('darkgreen')
	elif	i in both and 'a' in i:
		col.append('darkblue')
	elif i in spp_14 and 'p' in i:
		col.append('palegreen')
	elif	i in spp_14 and 'a' in i:
		col.append('lightblue')
	elif i in spp_15 and 'p' in i:
		col.append('lime')
	elif	i in spp_15 and 'a' in i:
		col.append('darkorchid')

g.vs['color'] = col
g.es['arrow_size'] = 0
lay = g.layout_bipartite()
for i in range(len(lay)):
	vname = g.vs[i]['name']
	if 'p' in vname:
		lay[i][1] = 1
	else:
		lay[i][1] = 0


visual_style = {}
visual_style["layout"] = lay
visual_style["bbox"] = (600, 300)
visual_style["margin"] = 10

plot(g,'nakamoto.pdf',**visual_style)


lay_dict = dict([[g.vs[i]['name'],lay[i]] for i in range(len(lay))])






year = ['2014','2015']
for net_n in range(2):
	g = Graph.TupleList(nets[net_n],directed=True)
	g.es['weight'] = 1
	g.simplify(combine_edges="sum")
	g.vs['type'] = [1 if 'p' in i else 0 for i in g.vs['name']]
	g.vs['color'] = ['darkgreen' if 'p' in i else 'darkblue' for i in g.vs['name']]
	g.es['color'] = ['lightgrey' if i==1 else 'black' for i in g.es['weight']]
	g.es['width'] = 2
	g.vs['size'] = 15#[20 if i in both else 10 for i in g.vs['name']]
	lay = g.layout_bipartite()
	for i in range(len(lay)):
		vname = g.vs[i]['name']
		lay[i] = lay_dict[vname]
	visual_style = {}
	visual_style["layout"] = lay
	visual_style["bbox"] = (600, 300)
	visual_style["margin"] = 10
	plot(g,'nakamoto_'+year[net_n]+'.pdf',**visual_style)


################chemello milazzo
h = [i for i in csv.reader(open('chemello_milazzo.csv','r'))]
h_mat = array([map(float,i) for i in h])
h_mat = h_mat.T #algae in row
h_mat*=0.1
plot_net(h_mat,'chemello_milazzo.pdf')


####Nestedness analysis



###load matrices/networks:
#huang
h = [i for i in csv.reader(open('huang.csv','r'))]
h_mat = array([map(float,i[1:-1]) for i in h[1:]])
h_mat = h_mat.T #algae in row
h_mat_1 = h_mat/h_mat.sum(0) #normalize by total interaction per animal
h_mat_1 = 1*(h_mat_1>0.295)###used 0.5 for plotting reasons
h_mat_2 = h_mat/h_mat.sum() #normalize by total interaction per animal
h_mat_2 = 1*(h_mat_2>0.0014)###used 0.5 for plotting reasons

#nakamoto
n_mat_1 = [i for i in csv.reader(open('nakamoto_mat_1.csv','r'))]
n_mat_1 = array([map(float,i) for i in n_mat_1])
n_mat_2 = [i for i in csv.reader(open('nakamoto_mat_2.csv','r'))]
n_mat_2 = array([map(float,i) for i in n_mat_2])


cm_file = [i for i in csv.reader(open('chemello_milazzo.csv','r'))]
cm_mat = array([map(float,i) for i in cm_file])
cm_mat = cm_mat.T #algae in row

all_mats = [h_mat_1,h_mat_2,cm_mat]

def nodf_an(mat):
	nodf = NODF(mat)
	null = array([NODF(FE(mat.copy())) for i in range(1000)])#keep number of invertebrate x macrophyte constant
	p_value = (null>=nodf).sum()/1000.0
	z = (nodf-mean(null))/std(null)
	return nodf,p_value,z


def plot_mat(m,file_name,size_fact=200):
	new_data = np.zeros(np.array(m.shape) * size_fact )
	for j in range(m.shape[0]):
		for k in range(m.shape[1]):
			new_data[j * size_fact: (j+1) * size_fact, k * size_fact: (k+1) * size_fact] = m[j, k]
	new_data = abs(1-new_data)
	plt.imsave(file_name, new_data, cmap=cm.gray)


for mat in all_mats:
	print(nodf_an(mat))


names = ['h_1','h_2','cm']
i = 0
for mat in all_mats:
	pm = pack(mat)
	plot_mat(pm,names[i]+'_packed.png')
	i+=1



#####make fully nested and random for visual comparison

def nest_mat(R,C,nes):
	mat = zeros([R,C])
	#mat.dtype = 'int16'
	for i in range(R):
		if C>i:
			for j in range(C-i):
				mat[i][j] = 1
		else:
			mat[i][0] = 1
	for i in range(int(round(R*C*(1-nes)))):
		rr,rc = randrange(R),randrange(C)
		mat[rr][rc] = abs(1-mat[rr][rc])
	return mat



def max_min_nest(mat_,reps=10000,mode='max'):
	mat = mat_.copy()
	R,C = mat.shape
	nodf0 = NODF(pack(mat))
	for i in range(reps):
		rr1,rc1 = randrange(R),randrange(C)
		rr2,rc2 = randrange(R),randrange(C)
		if mat[rr1][rc1]!=mat[rr2][rc2]:
			mat[rr1][rc1],mat[rr2][rc2] = mat[rr2][rc2],mat[rr1][rc1]
			new_nodf = NODF(pack(mat))
			if mode == 'max':
				if new_nodf>=nodf0 and 0 not in list(mat.sum(0))+list(mat.sum(1)):
					nodf0 = new_nodf
					print nodf0
				else:
					 mat[rr1][rc1],mat[rr2][rc2] = mat[rr2][rc2],mat[rr1][rc1]
			elif mode == 'min':
				if new_nodf<=nodf0 and 0 not in list(mat.sum(0))+list(mat.sum(1)):
					nodf0 = new_nodf
					print nodf0
				else:
					 mat[rr1][rc1],mat[rr2][rc2] = mat[rr2][rc2],mat[rr1][rc1]
			else:
				if 0 not in list(mat.sum(0))+list(mat.sum(1)):
					nodf0 = new_nodf
					print nodf0
				else:
					 mat[rr1][rc1],mat[rr2][rc2] = mat[rr2][rc2],mat[rr1][rc1]
	return mat


cm_nest = max_min_nest(cm_mat)
cm_min = max_min_nest(cm_mat,mode='min')
cm_rand = max_min_nest(cm_mat,mode='rand')


all_mats = [cm_mat,cm_rand,cm_nest,cm_min]
names = ['cm','cm_rand','cm_nest','cm_min']
i = 0
for mat in all_mats:
	pm = pack(mat)
	plot_mat(pm,names[i]+'.png')
	i+=1

i = 0
for mat in all_mats:
	pm = pack(mat)
	print (names[i],nodf_an(pm))
	i+=1


plot_mat(pack(cm_nest),'cm_nest.png')

#'cm', NODF = 61.2 (Z = 3.4; p<0.001)
#'cm_rand', NODF = 48.2 (Z = 0.64; p = 0.24)
#'cm_nest', NODF = 82.0 (Z = 5.0; p<0.001)
#'cm_min', NODF = 4.7 (Z = -17.6; p = 1.0)



#coextinction comparison
cm_file = [i for i in csv.reader(open('chemello_milazzo.csv','r'))]
cm_mat = array([map(float,i) for i in cm_file])
cm_mat = cm_mat.T #algae in row

rob_file = [i for i in csv.reader(open('robertson.csv','r'))]
rob_mat = array([map(float,i) for i in rob_file])
cm_mat = cm_mat.T #plants in row



def coe(mat_, mode = 'random'):
	r,c = mat_.shape
	div0 = float(sum(mat_.sum(0)>0))
	pl0 = float(sum(mat_.sum(1)>0))
	res = []
	for rep in range(100):
		res.append([1.0,1.0])
		mat = mat_.copy()
		while mat.sum()>0:
			if mode == 'random':
				nz_row = where(mat.sum(1)>0)[0]
			elif mode == 'best':
				row_sum = [i for i in mat.sum(1) if i>0]
				nz_row = where(mat.sum(1)==min(row_sum))[0]
			elif mode == 'worst':
				row_sum = [i for i in mat.sum(1) if i>0]
				nz_row = where(mat.sum(1)==max(row_sum))[0]
			rr = sample(nz_row,1)[0]
			mat[rr]*=0
			div1 = sum(mat.sum(0)>0)
			pl1 = sum(mat.sum(1)>0)
			res.append([pl1/pl0,div1/div0])
		print (rep)
	return res



res = []
for mode in ['random','best','worst']:
	res+=[[mode]+i for i in coe(cm_mat,mode = mode)]


out = open('coe_cm.csv','w')
for i  in res:
	out.write(','.join(map(str,i))+'\n')


out.close()



res = []
for mode in ['random','best','worst']:
	res+=[[mode]+i for i in coe(rob_mat,mode = mode)]


out = open('rob_cm.csv','w')
for i  in res:
	out.write(','.join(map(str,i))+'\n')


out.close()




