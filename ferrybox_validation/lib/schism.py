#!/usr/bin/env python
"""
SCHISM setup class
"""

__author__ = "Richard Hofmeister and Benjamin Jacob"
__copyright__   = "Copyright 2018 - 03\2021 Helmholtz-Zentrum Geesthacht"
__copyright__   = "Copyright 04\2021 - Helmholtz-Zentrum Hereon GmbH"

__license__ = "GNU GPL v2.0 (https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)"
__maintainer__ = "Benjamin Jacob"
__email__ = "benjamin.jacob@hzg.de"
__status__ = "Development"

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import PolyCollection
import datetime as dt
import xarray as xr
import glob
from mpl_toolkits.axes_grid1 import make_axes_locatable

class schism_setup(object):
  
  def __init__(self,hgrid_file='hgrid.gr3',ll_file='hgrid.ll',vgrid_file='vgrid.in'):
      self.hgrid_file=hgrid_file
      
      # parse hgrid file
      f = open(hgrid_file)
      self.description = f.readline().rstrip()
      dat = f.readline().split()
      self.nelements = int(dat[0])
      self.nnodes = int(dat[1])
      
      n=[]
      x=[]
      y=[]
      d=[]
      self.depthsdict={}
      self.xdict={}
      self.ydict={}
      for nn in range(self.nnodes):
        dat=f.readline().split()
        n.append(int(dat[0]))
        x.append(float(dat[1]))
        y.append(float(dat[2]))
        d.append(float(dat[3]))
        self.depthsdict[n[-1]] = d[-1]
        self.ydict[n[-1]] = y[-1]
        self.xdict[n[-1]] = x[-1]
      self.inodes = n
      self.x = x
      self.y = y
      self.depths = d
      
      n=[]
      nv = []
      nvdict = {}
      for nn in range(self.nelements):
        dat=f.readline().split()
        n.append(int(dat[0]))
        nvnum = int(dat[1])
        nv.append([ int(ii) for ii in dat[2:2+nvnum]])
        nvdict[n[-1]] = nv[-1]
      self.ielement = n
      self.nv = nv
      self.nvdict = nvdict

      # tri only grid 
      tris=[]
      quads=[]
      faces2=[]
      self.nvplt2nvp=[]	
      for nr,elem in enumerate(self.nv):
        if len(elem)==3:
          tris.append(elem)
          faces2.append(elem[:3])
          self.nvplt2nvp.append(nr)
        else: # split quad into tris
          quads.append(elem)
          faces2.append([ elem[0] ,elem[1],elem[2]])
          faces2.append([elem[0],elem[2],elem[3]])
          self.nvplt2nvp.append(nr)
          self.nvplt2nvp.append(nr)
      self.nvplt=np.asarray(faces2)-1 
      self.nvplt2nvp=np.asarray(self.nvplt2nvp)
      self.tris=np.asarray(tris)-1
      self.quads=np.asarray(quads)-1	  
	  
      # get resolution of elements
      res = {}
      dmin = {}	
      #import numpy as np
      for el in self.nvdict:
        inodes = self.nvdict[el]
        x = [self.xdict[ii] for ii in inodes]
        y = [self.ydict[ii] for ii in inodes]
        res[el] = (np.sqrt((x[1]-x[0])**2 + (y[1]-y[0])**2) + np.sqrt((x[2]-x[1])**2 + (y[2]-y[1])**2) + np.sqrt((x[0]-x[2])**2 + (y[0]-y[2])**2))/3.0
        dmin[el] = min(np.array( [np.sqrt((x[1]-x[0])**2 + (y[1]-y[0])**2), np.sqrt((x[2]-x[1])**2 + (y[2]-y[1])**2) , np.sqrt((x[0]-x[2])**2 + (y[0]-y[2])**2)]))
      self.resolution_by_element = dict(res)
      self.min_sidelength_by_element = dict(res)

      # compute sides
      # first get sequence array
      self.nx = {}
      self.nx[3] = [[1,2],[2,0],[0,1]]
      self.nx[4] = [[1,2],[2,3],[3,0],[0,1]]

      # get inverse nv - neighbouring elements
      self.node_neighbour_elements = { i:[] for i in self.inodes }
      self.n_node_neighbours = { i:0 for i in self.inodes }
      for i,nv in zip(self.ielement,self.nv):
        for nnode in nv:
          self.n_node_neighbours[nnode] += 1
          self.node_neighbour_elements[nnode].append(i)
      # find neighbouring elements around element (ic3)
      self.element_sides={}
      self.side_nodes={}
      for i,nv in zip(self.ielement,self.nv):
        isides = []
        # loop around element for existing sides
        for iloop,(ind1,ind2) in enumerate(self.nx[len(nv)]):
          iside = 0
          nd1,nd2 = nv[ind1],nv[ind2]
          for checkelement in self.node_neighbour_elements[nd1]:
            if (checkelement != i) and (nd2 in self.nvdict[checkelement]):
              iside = checkelement
              break
          isides.append(iside)
        self.element_sides[i] = isides
      # count sides
      self.nsides = 0
      element_ids = sorted(self.element_sides.keys())
      #element_ids.sort() python2 syntax
      for i in element_ids:
        for ii,iside in enumerate(self.element_sides[i]):
          if iside==0 or i<iside:
            self.nsides += 1
            iinds = self.nx[len(self.element_sides[i])][ii]
            self.side_nodes[self.nsides] = [self.nvdict[i][iinds[0]],self.nvdict[i][iinds[1]]]

      # average resolution_by_element on nodes
      self.resolution_by_nodes = {}
      for inode in self.node_neighbour_elements:
        elids = self.node_neighbour_elements[inode]
        self.resolution_by_nodes[inode] = np.asarray([self.resolution_by_element[ii] for ii in self.node_neighbour_elements[inode]]).mean()

      self.num_bdy_segments = int(f.readline().split()[0])
      self.num_bdy_nodes = int(f.readline().split()[0])
      self.bdy_segments=[]
      self.bdy_nodes=[]
      for iseg in range(self.num_bdy_segments):
        nlines = int(f.readline().split()[0])
        self.bdy_segments.append( [int(f.readline()) for nn in range(nlines)] )
        self.bdy_nodes.extend(self.bdy_segments[-1])

      self.num_land_segments = int(f.readline().split()[0])
      self.num_land_nodes = int(f.readline().split()[0])
      self.land_nodes=[]
      self.island_segments=[]
      self.land_segments=[]
      for iseg in range(self.num_land_segments):
        dat = f.readline().split()
        nlines = int(dat[0])
        if int(dat[1])==0:
          self.land_segments.append( [int(f.readline()) for nn in range(nlines)] )
          self.land_nodes.extend(self.land_segments[-1])
        elif int(dat[1])==1:
          self.island_segments.append( [int(f.readline()) for nn in range(nlines)] )
          self.land_nodes.extend(self.island_segments[-1])
     
      f.close()

      # parse hgrid.ll file
      try:
        f = open(ll_file)
        line = f.readline().rstrip()
        dat = f.readline().split()
        ll_nelements = int(dat[0])
        ll_nnodes = int(dat[1])
      
        nll = []
        lon=[]
        lat=[]
        self.londict={}
        self.latdict={}
        for nn in range(self.nnodes):
          dat=f.readline().split()
          nll.append(int(dat[0]))
          lon.append(float(dat[1]))
          lat.append(float(dat[2]))
          self.londict[nll[-1]]=lon[-1]
          self.latdict[nll[-1]]=lat[-1]
        self.ill = nll
        self.lon = lon
        self.lat = lat
        f.close()

        # compute minimum side length in 
        # element wise local projection (for cfl computation)
        res = {}
        dmin = {}
        #import numpy as np

        for el in self.ielement:
            inodes = self.nv[el-1]
            x,y=self.cpp_proj_elemNodes(inodes)
            # convert lon lat to local projection for element
            res[el] = (np.sqrt((x[1]-x[0])**2 + (y[1]-y[0])**2) + np.sqrt((x[2]-x[1])**2 + (y[2]-y[1])**2) + np.sqrt((x[0]-x[2])**2 + (y[0]-y[2])**2))/3.0
            dmin[el] = min(np.array( [np.sqrt((x[1]-x[0])**2 + (y[1]-y[0])**2), np.sqrt((x[2]-x[1])**2 + (y[2]-y[1])**2) , np.sqrt((x[0]-x[2])**2 + (y[0]-y[2])**2)]))
        self.cpp_resolution_by_element = res
        self.min_cpp_sidelength_by_element = dmin
        
      except:
        print('  no hgrid.ll available')

      try:
        self.parse_vgrid(vgrid_file=vgrid_file)
      except:
        print('  no vgrid.in available')
        #self.parse_vgrid(vgrid_file=vgrid_file)
        pass
      self.node_tree_xy = None
      self.node_tree_latlon = None
      self.element_tree_xy = None
      self.element_tree_latlon = None
      self.element_depth = None

  def parse_vgrid(self,vgrid_file='vgrid.in'):
      #import numpy as np
      f = open(vgrid_file)
      line=f.readline()
      if '!' in line:
        line=line.split('!')[0]	
      first = int(line)
      if first == 1: # unsructured vertical
        znum = int(f.readline())
        self.znum = znum
        a = {}
        self.bidx = {}
        for line in f.readlines():
          sigma1d = -9999.*np.ones((znum,))
          data = line.split()
          self.bidx[int(data[0])] = int(data[1])
          sigma1d[int(data[1])-1:] = np.asarray([float(ii) for ii in data[2:]])
          a[int(data[0])] = np.ma.masked_equal(sigma1d,-9999.)
        f.close()
        self.vgrid = a
      else: # sigma - z assume z
        self.znum,self.n_zlevels,self.h_s=f.readline().split(' ')[:3]
        self.znum,self.n_zlevels,self.h_s=int(self.znum),int(self.n_zlevels),np.float(self.h_s)
        f.readline()
        f.close()
        M=np.loadtxt(vgrid_file,skiprows=3,comments='!',max_rows=self.n_zlevels)
        if len(M.shape) > 1:
          self.zlevels=M[:,1]
        #else:		  
        #  self.zlevels=M
        self.slevels=np.loadtxt(vgrid_file,skiprows=5+self.n_zlevels,comments='!')[:,1]
        self.vgrid={}
        if self.n_zlevels > 1: #s_z
          self.vgrid_notes='vgrid, vertical  SZ levels are represented here as sigma levels with  resptec to total depths.'		
          comb=np.hstack((self.zlevels[:-1],self.h_s*self.slevels)) #z
          mask_sonly=np.hstack((np.ones(self.n_zlevels-1,bool),np.zeros(len(self.slevels),bool))) #z
          for i,d in enumerate(self.depths): #pure slevel when shallower?
            if d < self.h_s: # pure sigma ?		   
              #combi=np.hstack((self.zlevels[:-1],s.depths[i]*self.slevels))
			  #combi=np.hstack((self.zlevels[:-1],s.depths[i]*self.slevels)) #express as total depths part
              #combi=np.hstack((self.zlevels[:-1],self.h_s*self.slevels)) #z
              #self.vgrid[i+1]=np.ma.masked_array(combi,mask=comb<-d) # stores z
              #self.vgrid[i+1]=np.ma.masked_array(combi/d,mask=comb<-d) # expressed as sigma with respect to total depth#
              self.vgrid[i+1]=np.ma.masked_array(np.hstack((-np.ones(self.n_zlevels-1),self.slevels)),mask=mask_sonly)
            else:
              #self.vgrid[i+1]=np.ma.masked_array(comb,mask=comb<-d)
			  #self.vgrid[i+1]=np.ma.masked_array(comb/d,mask=comb<-d)
             # self.vgrid[i+1]=np.ma.masked_array(np.hstack((-np.ones(self.n_zlevels),self.slevels)),mask=mask_sonly)
			  
              combi=np.hstack((self.zlevels[:-1],self.h_s*self.slevels)) #z
              #self.vgrid[i+1]=np.ma.masked_array(combi,mask=comb<-d) # stores z
              self.vgrid[i+1]=np.ma.masked_array(combi/d,mask=comb<-d) # expressed as sigma with 
			  
        else: # pure sigma
          for i,d in enumerate(self.depths):
            #self.vgrid[i+1]=np.ma.masked_array(self.depths[i]*self.slevels,mask=False) # not times depths
            self.vgrid[i+1]=np.ma.masked_array(self.slevels,mask=False) # not times depths

			
  def dump_hgridll(self,filename='hgrid_new.ll'):
    f = open(filename,'w')
    f.write('%s\n'%filename)
    f.write('%d %d\n'%(self.nelements,self.nnodes))
    # write nodes
    for n,x,y,d in zip(self.inodes,self.lon,self.lat,self.depths):
      f.write('%d %0.12f %0.12f %0.12f\n'%(n,x,y,d))

    # write elements
    for n,nv in zip(self.ielement,self.nv):
      f.write('%d %d '%(n,len(nv)))
      for nvi in nv:
        f.write('%d '%nvi)
      f.write('\n')
    # write boundaries
    f.close()


  def dump_hgridgr3(self,filename='hgrid_new.gr3',comment='grid_transfer'):
      f = open(filename,'w')
      f.write('%s - %s\n'%(filename,comment))
      f.write('%d %d\n'%(self.nelements,self.nnodes))
      # write nodes
      for n,x,y,d in zip(self.inodes,self.x,self.y,self.depths):
          f.write('%d %0.12f %0.12f %0.12f\n'%(n,x,y,d))

      # write elements
      for n,nv in zip(self.ielement,self.nv):
          f.write('%d %d '%(n,len(nv)))
          for nvi in nv:
              f.write('%d '%nvi)
          f.write('\n')

      # write boundaries
      f.write('%i = Number of open boundaries\n'%self.num_bdy_segments)
      f.write('%i = Total number of open boundary nodes\n'%self.num_bdy_nodes)

      for iseg,seg in enumerate(self.bdy_segments):
          f.write('%i 0 = Number of nodes for open boundary %i\n'%(len(seg),iseg+1))
          for inode in range(len(seg)):			
              f.write('%i\n'%(seg[inode]))

      f.write('%i = number of land boundaries\n'%self.num_land_segments)
      f.write('%i = Total number of land boundary nodes\n'%self.num_land_nodes)

      for iseg,seg in enumerate(self.land_segments):
          f.write('%i 0 = Number of nodes for land boundary %i\n'%(len(seg),iseg+1))
          for inode in range(len(seg)):			
              f.write('%i\n'%(seg[inode]))

      for iseg,seg in enumerate(self.island_segments):
          f.write('%i 1 = Number of nodes for island boundary %i\n'%(len(seg),iseg+1))
          for inode in range(len(seg)):			
              f.write('%i\n'%(seg[inode]))
      f.close()

  def dump_vgrid(self,filename='vgrid_new.in'):
      with open(filename,'w') as f:
          f.write('{:8d}\n{:8d}'.format(1,len(self.vgrid[1])))
          nlevels=len(self.vgrid[1])
          for node in self.vgrid.keys():
              n=self.vgrid[node].count()
              level_start=nlevels-n+1
              outarray=self.vgrid[node][-n:]
              f.write('\n')
              f.write('{:8d} {:8d}'.format(node,level_start))
              for val in outarray:
                  f.write(' {:12f}'.format(val))


#  def get_bdy_latlon(self):
#      bdylon = [ self.lon[ii-1] for ii in self.bdy_nodes ]
#      bdylat = [ self.lat[ii-1] for ii in self.bdy_nodes ]
#      return (bdylon,bdylat)
  def get_bdy_latlon(self):
      bdnodes=[]
      for land,ocean in list(zip(self.land_segments,self.bdy_segments)):		
        bdnode.append(ocean)
        bdnode.append(land[1:])
      bdnodes=np.hstack(bdnodes)		
      bdylon=np.asarray(self.lon)[bdnodes-1]
      bdylat=np.asarray(self.lat)[bdnodes-1]
      return (bdylon,bdylat)

	  
  def plot_domain_boundaries(self,append=0,nr=0,latlon=True,plot_legend=True,obd_clr='r',lcol = (0.6,0.6,0.6),islcol = (0.3,0.3,0.3)):
      from pylab import figure,plot,show,legend,xlabel,ylabel
      if append==0:	
	      f = figure()

      if latlon==True:	
	      x,y=self.lon,self.lat
	      labelx='longitude [degE]'		  
	      labely='latitude [degN]'
      else:	
	      x,y=self.x,self.y	 
	      labelx='x [m]'		  		  
	      labely='y [m]'		  		  		  
		  
      label='open boundary'
      for seg in self.bdy_segments:
        lon = [ x[ii-1] for ii in seg ]
        lat = [ y[ii-1] for ii in seg ]
        plot( lon,lat, obd_clr+'o-', ms=0.5, label=label, markeredgecolor='r')
        label=''

      if nr:		
        for i,seg in enumerate(self.bdy_segments):
          lon = np.asarray([ x[ii-1] for ii in seg ]).mean()
          lat = np.asarray([ y[ii-1] for ii in seg ]).mean()
          plt.text(lon,lat,'bd '+str(i))
		
		
      label='land boundary'
      for seg in self.land_segments:
        lon = [ x[ii-1] for ii in seg ]
        lat = [ y[ii-1] for ii in seg ]
        plot( lon,lat, 'o-',color=lcol, ms=0.1, label=label, markeredgecolor=lcol)
        label=''
 
      label='island boundary'
      for seg in self.island_segments:
        lon = [ x[ii-1] for ii in seg ]
        lat = [ y[ii-1] for ii in seg ]
        plot( lon,lat, 'o-',color=lcol, ms=0.1, label=label, markeredgecolor=islcol)
        label=''

      if plot_legend==True:		
        legend(loc='lower right',frameon=False)
      xlabel(labelx)
      ylabel(labely)
      show()

  def plot_mesh(self,tri_color='k',quad_color='m',linewidth=0.2,latlon=True):	  
      if latlon==True:  
          xy=np.c_[self.lon,self.lat]
      else:
          xy=np.c_[self.x,self.y]
      tripc = PolyCollection(xy[self.tris,:3],facecolors='none',edgecolors=tri_color,linewidth=linewidth) #, **kwargs)
      phtri=plt.gca().add_collection(tripc)
      if(len(self.quads)):
          quadpc = PolyCollection(xy[self.quads,:4],facecolors='none',edgecolors=quad_color,linewidth=linewidth) #, **kwargs)
          phquad=plt.gca().add_collection(quadpc)	  
		  

  def signed_area(self,nodelist):
      x1,y1 = (self.xdict[nodelist[0]],self.ydict[nodelist[0]])
      x2,y2 = (self.xdict[nodelist[1]],self.ydict[nodelist[1]])
      x3,y3 = (self.xdict[nodelist[2]],self.ydict[nodelist[2]])
      return ((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2.0
  
    
  def proj_area(self,nodelist):
      #import numpy as np
      x1,y1 = (self.londict[nodelist[0]],self.latdict[nodelist[0]])
      x2,y2 = (self.londict[nodelist[1]],self.latdict[nodelist[1]])
      x3,y3 = (self.londict[nodelist[2]],self.latdict[nodelist[2]])
      #
      rad=np.pi/180
      rlambda=(x1+x2+x3)/3.0*rad # x,y correct order
      phi=(y1+y2+y3)/3.0*rad
      #
      r=6378206.4
      x1=r*(x1*rad-rlambda)*np.cos(phi)
      x2=r*(x2*rad-rlambda)*np.cos(phi)
      x3=r*(x3*rad-rlambda)*np.cos(phi)
      #
      y1=r*y1*rad
      y2=r*y2*rad
      y3=r*y3*rad
      return ((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2.0  


  def dump_gr3(self,filename,const=0.0,comment='gr3 by create_ic.py'):
      f = open(filename,'w')
      f.write('%s\n'%comment)
      f.write('%d %d\n'%(self.nelements,self.nnodes))
      for i,x,y in zip(self.inodes,self.x,self.y):
        f.write('%d %0.12f %0.12f %0.12f\n'%(i,x,y,const))
      f.close()

  def dump_gr3_spat_var(self,filename,nodevalues,comment='gr3 by create_ic.py'):
      f = open(filename,'w')
      f.write('%s\n'%comment)
      f.write('%d %d\n'%(self.nelements,self.nnodes))
      for i,x,y,value in zip(self.inodes,self.x,self.y,nodevalues):
        f.write('%d %0.12f %0.12f %0.12f\n'%(i,x,y,value))
      f.close()

  def read_gr3(self,filename):
      if not hasattr(self,'gr3'):
        self.gr3={}
      self.gr3[filename.split('.')[0]] = np.loadtxt(filename,skiprows=2,max_rows=self.nnodes)[:,-1]

  def read_reg(self,filename):
      if not hasattr(self,'reg'):
        self.reg={}
      self.reg[filename.split('.')[0]] = np.loadtxt(filename,skiprows=3)
	  

  def dump_tvd_prop(self):
      f = open('tvd.prop','w')
      for i in self.ielement:
        f.write('%d 1\n'%i)
      f.close()


  def init_node_tree(self,latlon=True):
    #print('  build node tree')
    from scipy.spatial import cKDTree
    if latlon:
      self.node_tree_latlon = cKDTree(list(zip(self.lon,self.lat)))
    else:
      self.node_tree_xy = cKDTree(list(zip(self.x,self.y)))

  def init_element_tree(self,latlon=True):
    """
    build element tree for xy or latlon coordinates
    (default: latlon=True)
    """
    from scipy.spatial import cKDTree
    
    #print('  schism.py: build element tree')
    if self.element_depth == None:
      self.element_depth={}
      calc_depths = True
    else:
      calc_depths = False

    if latlon:
      self.element_lon={}
      self.element_lat={}
      for el in self.nvdict:
        self.element_lon[el] = sum([self.londict[idx] for idx in self.nvdict[el]])/len(self.nvdict[el])
        self.element_lat[el] = sum([self.latdict[idx] for idx in self.nvdict[el]])/len(self.nvdict[el])
        if calc_depths:
          self.element_depth[el] = sum([self.depthsdict[idx] for idx in self.nvdict[el]])/len(self.nvdict[el])
      self.element_tree_latlon = cKDTree(list(zip(self.element_lon.values(),self.element_lat.values())))
      self.element_tree_ids = list(self.element_lon.keys())
    else:
      self.element_x={}
      self.element_y={}
      for el in self.nvdict:
        self.element_x[el] = sum([self.xdict[idx] for idx in self.nvdict[el]])/len(self.nvdict[el])
        self.element_y[el] = sum([self.ydict[idx] for idx in self.nvdict[el]])/len(self.nvdict[el])
        if calc_depths:
          self.element_depth[el] = sum([self.depthsdict[idx] for idx in self.nvdict[el]])/len(self.nvdict[el])
      self.element_tree_xy = cKDTree(list(zip(self.element_x.values(),self.element_y.values())))
      self.element_tree_ids = list(self.element_x.keys())


  def find_nearest_node(self,x,y,latlon=True):
    """
    find nearest node for given coordinate,
    returns the node id
    """
    ridx=-1
    if latlon:
      if self.node_tree_latlon==None:
        self.init_node_tree(latlon=True)
      d,idx = self.node_tree_latlon.query((x,y),k=1)
      ridx = self.ill[idx]
    else:
      if self.node_tree_latlon==None:
        self.init_node_tree(latlon=False)
      d,idx = self.node_tree_latlon.query((x,y),k=1)
      ridx = self.inodes[idx]
    return ridx

  def find_nearest_element(self,x,y,latlon=True,mindepth=3.0,numcells=5):
    """
    give coordinates and find nearest element,
    returns the element id
    """
    ridx=-1
    numcells=numcells
    while ridx<0:
      if latlon:
        if self.element_tree_latlon==None:
          self.init_element_tree(latlon=True)
        d,idx = self.element_tree_latlon.query((x,y),k=numcells)
      else:
        if self.element_tree_xy==None:
          self.init_element_tree(latlon=False)
        d,idx = self.element_tree_xy.query((x,y),k=numcells)
      depths = np.array([self.element_depth[self.element_tree_ids[searchid]] for searchid in idx])
      maxidx = np.argmax(depths)
      if depths[maxidx]>mindepth:
        ridx = self.element_tree_ids[idx[maxidx]]
        #print('  info: found cell with depth>%0.1fm in %d nearest cells'%(mindepth,numcells))
      else:
        # continue iterating:
        numcells = numcells*2
        if numcells > 10000:
          print('  warning: cannot find cell with depth>%0.1fm'%(mindepth))
          break
    return ridx



  def cpp_proj_elemNodes(self,nodelist):
      # create cntral point porjection for element of lon lat        
      #import numpy as np
      x1,y1 = (self.londict[nodelist[0]],self.latdict[nodelist[0]])
      x2,y2 = (self.londict[nodelist[1]],self.latdict[nodelist[1]])
      x3,y3 = (self.londict[nodelist[2]],self.latdict[nodelist[2]])
      #
      rad=np.pi/180
      rlambda=(x1+x2+x3)/3.0*rad # x,y correct order
      phi=(y1+y2+y3)/3.0*rad
      #
      r=6378206.4
      x1=r*(x1*rad-rlambda)*np.cos(phi)
      x2=r*(x2*rad-rlambda)*np.cos(phi)
      x3=r*(x3*rad-rlambda)*np.cos(phi)
      #
      y1=r*y1*rad
      y2=r*y2*rad
      y3=r*y3*rad
      
      return [x1,x2,x3],[y1,y2,y3]    



  def compute_cfl(self,dt,grid='gr3'):

      
    #import numpy as np
    
    if grid=='cpp':
        dxs=np.asarray(list(self.min_cpp_sidelength_by_element.values()))
    elif grid=='gr3':
        dxs=np.asarray(list(self.min_sidelength_by_element.values()))
    else:
        print('specify grid as either gr3 or cpp')
        return

    g=9.81
    cfl_by_elements=np.zeros(self.nelements)    
    
    # calc depth at nodes
    maxDepth=np.zeros(self.nelements)
    
    for ielement in self.ielement:
        maxDepth[ielement-1]=max([self.depths[idx-1] for idx in self.nv[ielement-1]])
        
    
    
    #for ielem in self.ielement:
    cfl_by_elements=np.sqrt(g*maxDepth)*dt/dxs
    
    return cfl_by_elements

  def write_hotstart(self,tr_nd,filename='hotstart.nc',time=0.0,iths=0,ifile=0,elev=0.0,uvel=0.0,vvel=0.0):
    """
    write hotstart.nc from given tracer concentrations on nodes
    tracer concentrations on sides and elements will be interpolated
    tr_nd.shape has to be (nelements,nvrt,ntracers)
    """
    inum,znum,ntracers = tr_nd.shape
    if ((inum,znum) != (self.nnodes,self.znum)):
      print('  shape(tr_nd) = (%d,%d) while setup requires (%d,%d)'%(inum,znum,self.nnodes,self.znum))

    import netCDF4
    nc = netCDF4.Dataset(filename,'w',format='NETCDF4_CLASSIC')
    nc.createDimension('node',self.nnodes)
    nc.createDimension('elem',self.nelements)
    nc.createDimension('side',self.nsides)
    nc.createDimension('nVert',self.znum)
    nc.createDimension('ntracers',ntracers)
    nc.createDimension('one',1)
    nc.createDimension('three',3)

    v = nc.createVariable('time','f8',('one',))
    v[:] = time
    v = nc.createVariable('iths','i',('one',))
    v[:] = iths 
    v = nc.createVariable('ifile','i',('one',))
    v[:] = ifile
    v = nc.createVariable('idry_e','i',('elem',))
    v[:] = 0
    v = nc.createVariable('idry_s','i',('side',))
    v[:] = 0
    v = nc.createVariable('idry','i',('node',))
    v[:] = 0
    v = nc.createVariable('eta2','f8',('node',))
    v[:] = elev
    v = nc.createVariable('we','f8',('elem','nVert'))
    v[:] = 0.0
    v = nc.createVariable('su2','f8',('side','nVert'))
    v[:] = uvel
    v = nc.createVariable('sv2','f8',('side','nVert'))
    v[:] = vvel
    v = nc.createVariable('q2','f8',('node','nVert'))
    v[:] = 0.0
    v = nc.createVariable('xl','f8',('node','nVert'))
    v[:] = 0.0
    v = nc.createVariable('dfv','f8',('node','nVert'))
    v[:] = 0.0
    v = nc.createVariable('dfh','f8',('node','nVert'))
    v[:] = 0.0
    v = nc.createVariable('dfq1','f8',('node','nVert'))
    v[:] = 0.0
    v = nc.createVariable('dfq2','f8',('node','nVert'))
    v[:] = 0.0
    nc.sync()

	# add for new schism
    v = nc.createVariable('nsteps_from_cold','i',('one',))
    v[:] = 0.0#
    v = nc.createVariable('cumsum_eta','i',('one',))
    v[:] = 0.0
    nc.sync()

    # write tracer concentrations on nodes
    v = nc.createVariable('tr_nd','f8',('node','nVert','ntracers'))
    v[:] = tr_nd
    nc.sync()
    v = nc.createVariable('tr_nd0','f8',('node','nVert','ntracers'))
    v[:] = tr_nd
    nc.sync()

    # write tracer concentrations on elements
    v = nc.createVariable('tr_el','f8',('elem','nVert','ntracers'))
    for ie in self.nvdict:
      inds = self.nvdict[ie]
      tr_coll = [tr_nd[ind-1] for ind in inds]
      tr_el = np.mean(tr_coll,axis=0)
      v[ie-1] = tr_el
    nc.sync()
    nc.close()

	
  def compute_element_areas():
    """ compute  element areas (quads are splite before) """  
    # Spatial Average
    A=[]
    for i in range(self.nvplt.shape[0]):
      nodes=s.nvplt[i,:]+1
      A.append(s.proj_area(nodes))	
    return A

	  
  def write_bdy_netcdf(self,filename,time,data,frcbdnodes=[]):
      """
      write boundary data for schism setup
      """
      import netCDF4

      if self.num_bdy_nodes==0:
        print('  setup has no open boundaries')
        return

      datadims = len(data.shape)
      if datadims==2:
        tnum,nbdy = data.shape
        znum = 1
        ncom = 1
      elif datadims==3:
        tnum,nbdy,ncom = data.shape
        znum = 1
      elif datadims==4:
        tnum,nbdy,znum,ncom = data.shape

      if len(frcbdnodes)>0:
        nbdndes=len(frcbdnodes)
      else:
        nbdndes=self.num_bdy_nodes		

      nc = netCDF4.Dataset(filename,'w',format='NETCDF3_CLASSIC')
      nc.createDimension('time',None)
      nc.createDimension('nOpenBndNodes',nbdndes)
      nc.createDimension('nLevels',znum)
      nc.createDimension('nComponents',ncom)
      nc.createDimension('one',1)
      
      v = nc.createVariable('time_step','f4',('one',))
      v.long_name = 'time step in seconds'
      v[:] = time[1]-time[0]

      v = nc.createVariable('time','f8',('time',))
      v.long_name = 'simulation time in seconds'
      v[0:tnum] = time

      v = nc.createVariable('time_series','f4',('time','nOpenBndNodes','nLevels','nComponents'))
      if datadims==2:
        v[0:tnum,0:nbdy,0,0] = data
      elif datadims==3:
        v[0:tnum,0:nbdy,0,0:ncom] = data
      elif datadims==4:
        #v[0:tnum] = data
        v[0:tnum,0:nbdy,0:znum,0:ncom] = data
	
      nc.close()
      return


  # plot functions 
  def plotAtnodes(self,nodevalues,cmap=plt.cm.viridis,mask=None,latlon=True,extend='neither',ax=None):
      """
      visualisation routine plotting triangles at nodes (quads are splitted)
      """
	        # cartopy projection and features
      if ax == None: # else ax is given for eg sub axe üöots
        ax = plt.axes()
	  
      if latlon==True:	        
        ph=ax.tripcolor(self.lon,self.lat,self.nvplt,nodevalues,shading='flat',mask=mask,cmap=cmap)
      else:	              
        ph=ax.tripcolor(self.x,self.y,self.nvplt,nodevalues,shading='flat',mask=mask,cmap=cmap)	  
      #ch=plt.colorbar(ph,extend=extend)	
      # adjust height of colorbar to fit plot axes
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05, axes_class=plt.Axes)
      ch=plt.colorbar(ph,cax=cax,extend=extend)

      return ph,ch,ax

  def plotAtelems(self,values,cmap=plt.cm.viridis,mask=None,latlon=True,extend='neither'):
      """
      visualisation routine plotting at triangles (quads are splitted)
      """
      if latlon==True:	        	  
        x,y=self.lon,self.lat
      else:	              
        x,y=self.x,self.y
		
      if len(values)==self.nnodes:  	
          ph=plt.tripcolor(x,y,self.nvplt,facecolors=values[self.nvplt[:,:3]].mean(axis=1),shading='flat',cmap=cmap)# shading needs gouraud to allow correct update
      elif len(values)==len(self.nvplt):
          ph=plt.tripcolor(x,y,self.nvplt,facecolors=values,shading='flat',mask=mask,cmap=cmap)# shading needs
      ch=plt.colorbar(extend=extend)
      plt.tight_layout()
      return ph,ch

	  
  # plot functions - using basemap
  def plotAtnodesGeo(self,values,cmap=plt.cm.jet,mask=None,proj='merc',offset=0.1,stock_image=False,extend='both',region_limit=None,drycolor='grey',ax=None):
      """	
      visualisation routine plotting data at nodes (quads are splitted) and use cartopy map to draw a map
      valid projections are merc:=mercator and stere:=stereographic      plotAtnodesGeo(self,values,cmap=plt.cm.jet,mask=None,proj='merc',offset=0.1,stock_image=False,extend='both',region_limit=(lonmin,lonmax,latmin,latmax) or None,drycolor='grey',ax= geoaxis handle for subbplots with cartopy):
	  """

      #cmap=plt.cm.jet # use colormap
      cmap.set_bad(drycolor,1) # set dry area value
	  
      lon,lat = np.asarray(self.lon),np.asarray(self.lat)
      from mpl_toolkits.axes_grid1 import make_axes_locatable # allow axis adjustment of colorbar	  
      # cartopy
      import cartopy.crs as ccrs
      import cartopy.feature as cfeature
      from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter


      ### cartopy ########
      if proj=='merc':
          proj=ccrs.Mercator()  # define Prijection
      ## load higher resolutions coastline assets
      else:
          proj=ccrs.PlateCarree()  # define Prijection
		  
      land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',edgecolor='face',                                        facecolor=cfeature.COLORS['land'])
      ocean_10m = cfeature.NaturalEarthFeature('physical', 'ocean', '10m', edgecolor='face',facecolor=[0,0,1])
      #zoom_extend=(np.min(s.lon)-0.1, np.max(s.lon)+0.02, np.min(s.lat)-0.04, np.max(s.lat)+0.1)
	  
      if region_limit==None:
        zoom_extend=(np.min(self.lon)-offset,np.max(self.lon)+offset, np.min(self.lat)-offset, np.max(self.lat)+offset)
      else:      
        zoom_extend=region_limit
      ### project coordinates of schism class used in internal plot functions to map projection
      print(zoom_extend)
      outproj=proj.transform_points(ccrs.Geodetic(),lon,lat)
      self.projx,self.projy=outproj[:,0],outproj[:,1]
      ####################

      # cartopy projection and features
      if ax == None: # else ax is given for eg sub axe üöots
          ax = plt.axes(projection=proj)
		
		
      ax.set_extent(zoom_extend)
      ax.add_feature(land_10m,zorder=-2)

      if (proj!='merc') & stock_image:
          ax.stock_img()
      #ph=plt.tripcolor(self.projx,self.projy,self.nvplt,nodevalues,shading='flat',mask=mask,cmap=cmap)

      if len(values)==self.nnodes: #test ax instead plt here 	
          ph=ax.tripcolor(self.projx,self.projy,self.nvplt,facecolors=np.ma.masked_array(values[self.nvplt[:,:3]].mean(axis=1),mask=mask),shading='flat',cmap=cmap)# shading needs gouraud to allow correct update
      elif len(values)==len(self.nvplt):
          ph=ax.tripcolor(self.projx,self.projy,self.nvplt,facecolors=np.ma.masked_array(values,mask=mask),shading='flat',cmap=cmap)# shading needs
      #ch=plt.colorbar(extend=extend)
      plt.tight_layout()

	  
      # scale colorbar to 1 and 99% quantile
      vmin=np.round(np.quantile(values,0.01),1)
      vmax=np.round(np.quantile(values,0.99),1)
      ph.set_clim((vmin,vmax))
      
      # Set tick labels
      #xticks=np.unique(np.round(np.linspace((lon.min()),(lon.max()),8),1))
      #yticks=np.unique(np.round(np.linspace((lat.min()),(lat.max()),8),1))
      ax.set_extent(zoom_extend)
      longrange=np.max(self.lon)-np.min(self.lon)
      latrange=np.max(self.lat)-np.min(self.lat)
      ratio=latrange/longrange
      if ratio > 1:
        nxticks=int(np.floor(8/ratio))
        nyticks=8
      else:
        nxticks=8
        nyticks=int(np.floor(8*ratio))
      
      xticks=np.unique(np.round(np.linspace((lon.min()),(lon.max()),nxticks),1))
      yticks=np.unique(np.round(np.linspace((lat.min()),(lat.max()),nyticks),1))


      ax.set_xticks(xticks, crs=ccrs.PlateCarree())
      ax.set_yticks(yticks, crs=ccrs.PlateCarree())
      lon_formatter = LongitudeFormatter(number_format='.1f',degree_symbol='',dateline_direction_label=True)
      lat_formatter = LatitudeFormatter(number_format='.1f',degree_symbol='')
      ax.xaxis.set_major_formatter(lon_formatter)
      ax.yaxis.set_major_formatter(lat_formatter)
      
      # adjust height of colorbar to fit plot axes
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05, axes_class=plt.Axes)
      ch=plt.colorbar(ph,cax=cax,extend=extend)
      
      return ph,ch,ax
	  
  def find_parent_tri(self,xq,yq,dThresh=1000,latlon=False):
    """ parents,ndeweights=find_parent_tri(xq,yq,dThresh=1000)
    	find parent for coordinates xq,yq within triangulation tris,xun,yun.
    	return: parent triangle ids (quads are splitted into triangles in a previous step) and barycentric weights of triangle coordinates.
    """    
    if latlon:
      xun=np.asarray(self.lon)
      yun=np.asarray(self.lat)
      dThresh=np.min((dThresh,1))
    else:
      xun=np.asarray(self.x)
      yun=np.asarray(self.y)
	
    #% Distance threshold for Point distance
    dThresh=dThresh**2
    tris=np.asarray(self.nvplt)
    trinr=np.arange(tris.shape[0])
    
    trisX,trisY=xun[tris],yun[tris]
    #% orthogonal of side vecotrs
    SideX=np.diff(trisY[:,[0, 1, 2, 0]],axis=1)
    SideY=-np.diff(trisX[:,[0, 1, 2, 0]],axis=1)
    
    p=np.stack((xq,yq),axis=1)
    parent=-1*np.ones(len(p),int)
    for ip in range(len(p)):
      
      dx1=(p[ip,0]-trisX[:,0])
      dy1=(p[ip,1]-trisY[:,0])
      subind=(dx1*dx1+dy1*dy1) < dThresh # preselection
      subtris=trinr[subind]
      
      #% dot products
      parenti=(subtris[ (dx1[subind]*SideX[subind,0] + dy1[subind]*SideY[subind,0] <= 0) \
      			   & ((p[ip,0]-trisX[subind,1])*SideX[subind,1] + (p[ip,1]-trisY[subind,1])*SideY[subind,1] <= 0) \
      				 & ( (p[ip,0]-trisX[subind,2])*SideX[subind,2] + (p[ip,1]-trisY[subind,2])*SideY[subind,2] <= 0) ][:])
      if len(parenti):
        parent[ip]=parenti
      
    # tri nodes
    xabc=xun[tris[parent]]
    yabc=yun[tris[parent]]
    
    # barycentric weights
    divisor=(yabc[:,1]-yabc[:,2])*(xabc[:,0]-xabc[:,2])+(xabc[:,2]-xabc[:,1])*(yabc[:,0]-yabc[:,2])
    w1=((yabc[:,1]-yabc[:,2])*(xq-xabc[:,2])+(xabc[:,2]-xabc[:,1])*(yq-yabc[:,2]))/divisor
    w2=((yabc[:,2]-yabc[:,0])*(xq-xabc[:,2])+(xabc[:,0]-xabc[:,2])*(yq-yabc[:,2]))/divisor
    w3=1-w1-w2
    
    return parent,np.stack((w1,w2,w3)).transpose() 	  
	  

class schism_output():
    import netCDF4
    nc = None

    def __init__(self,filename):
      """
      read output filename and initialize grid
      """
      import netCDF4
      from netcdftime import utime
      self.nc = netCDF4.Dataset(filename)
      self.ncv = self.nc.variables
      self.lon = self.ncv['SCHISM_hgrid_node_x'][:]
      self.lat = self.ncv['SCHISM_hgrid_node_y'][:]
      self.nodeids = np.arange(len(self.lon))
      self.nv = self.ncv['SCHISM_hgrid_face_nodes'][:,:3]-1
      self.time = self.ncv['time'][:] # s
      self.ut = utime(self.ncv['time'].units)
      self.dates = self.ut.num2date(self.time)
      self.node_tree_latlon = None

    def init_node_tree(self,latlon=True):
      """
      build a node tree using cKDTree
      for a quick search for node coordinates
      """
      from scipy.spatial import cKDTree
      if latlon:
        self.node_tree_latlon = cKDTree(list(zip(self.lon,self.lat)))
      else:
        self.node_tree_xy = cKDTree(list(zip(self.x,self.y)))

    def find_nearest_node(self,x,y,latlon=True):
      """
      find nearest node for given coordinate,
      returns the node id
      """
      ridx=-1
      if latlon:
        if self.node_tree_latlon==None:
          self.init_node_tree(latlon=True)
        d,idx = self.node_tree_latlon.query((x,y),k=1)
        ridx = self.nodeids[idx]
      else:
        if self.node_tree_latlon==None:
           self.init_node_tree(latlon=False)
        d,idx = self.node_tree_latlon.query((x,y),k=1)
        ridx = self.inodes[idx]
      return ridx

	  
	  
class schism_output2():
    nc = None

    def __init__(self,ncdir):
      """
      read output files in folder and initialize grid
      """
      import xarray as xr      
      import glob
      schismfiles=[] 
      for iorder in range(8): # check for schout_nc files until 99999
        schismfiles+=glob.glob(ncdir+'schout_'+'?'*iorder+'.nc')
        nrs=[int(file[file.rfind('_')+1:file.index('.nc')]) for file in schismfiles]
        schismfiles=list(np.asarray(schismfiles)[np.argsort(nrs)])
        nrs=list(np.asarray(nrs)[np.argsort(nrs)])


      try:
        #self.nc = xr.open_mfdataset(schismfiles)
        self.nc =xr.concat( [xr.open_dataset(fname).chunk() for fname in schismfiles], dim='time')
		
      except:
        print('error loading, trying now leaving out last stack '+ schismfiles[-1])	  
        #self.nc = xr.open_mfdataset(schismfiles[:-1])
        self.nc =xr.concat( [xr.open_dataset(fname).chunk() for fname in schismfiles[:-1]], dim='time')
      self.ncv = self.nc.variables
      self.lon = self.ncv['SCHISM_hgrid_node_x'][:].values
      self.lat = self.ncv['SCHISM_hgrid_node_y'][:].values
      self.nodeids = np.arange(len(self.lon))
      self.nv = self.ncv['SCHISM_hgrid_face_nodes'][:,:3]-1
      self.time = self.ncv['time'][:].values # s
      #self.ut = utime(self.ncv['time'].units)
      self.dates = self.time #self.ut.num2date(self.time)
      self.node_tree_latlon = None

    def init_node_tree(self,latlon=True):
      """
      build a node tree using cKDTree
      for a quick search for node coordinates
      """
      from scipy.spatial import cKDTree
      if latlon:
        self.node_tree_latlon = cKDTree(list(zip(self.lon,self.lat)))
      else:
        self.node_tree_xy = cKDTree(list(zip(self.x,self.y)))

    def find_nearest_node(self,x,y,latlon=True):
      """
      find nearest node for given coordinate,
      returns the node id
      """
      ridx=-1
      if latlon:
        if self.node_tree_latlon==None:
          self.init_node_tree(latlon=True)
        d,idx = self.node_tree_latlon.query((x,y),k=1)
        ridx = self.nodeids[idx]
      else:
        if self.node_tree_latlon==None:
           self.init_node_tree(latlon=False)
        d,idx = self.node_tree_latlon.query((x,y),k=1)
        ridx = self.inodes[idx]
      return ridx	  

	  
class bp_file():
  def __init__(self,bpfile,lonlat=False):
      """ read bp file """  
      trans=np.loadtxt(bpfile,skiprows=2)[:,1:3]
      self.coords=list(zip(trans[:,0],trans[:,1]))
      if lonlat==False:
        self.x,self.y=trans[:,0],trans[:,1]
        coords=list(zip(self.x,self.y))
        p=np.asarray(coords)
        dl=np.sqrt((np.diff(p,axis=0)**2).sum(axis=1))
        self.l=np.hstack((0,np.cumsum(dl)))
      else:        		        
        self.lon,self.lat=trans[:,0],trans[:,1]
	  
	  
class bp_transect_nn():
 
    def __init__(self,s,ds,x,y,latlon=False,nn=False):
      """ bp_transect(s,ds,x,y,lonlat=False,nn=False) s:=schism_setup(), ds:= xarray handle ,x,y coordinates
	  Constraut transect object with track tangent normal and tangents, 
	  and along across projection function. If nn = True limit input coordinates
	  to those with different nearest neighbouts, only, ie to 
	  get xarray output access as nearest neighbour for transect coordinates. """
      self.x,self.y=x,y	  
      s.init_node_tree(latlon=latlon)
      coords=list(zip(x,y))
      p=np.asarray(coords)
      dd,self.nn=s.node_tree_xy.query(coords)
      self.nn,iunq=np.unique(self.nn,return_index=True)
      self.x,self.y=self.x[iunq],self.y[iunq]
      self.ds=ds.sel(nSCHISM_hgrid_node=self.nn)
      if nn:
          p=p[iunq,:]
      # projection for velocties
      dl=np.sqrt((np.diff(p,axis=0)**2).sum(axis=1))
      #dl=np.sqrt(np.diff(x)**2+np.diff(y)**2)
      self.l=np.hstack((0,np.cumsum(dl)))
      z=ds['zcor']
      self.l=np.tile(self.l,[z.shape[-1],1]).T
      tangent=np.diff(p,axis=0)/np.tile(dl,(2,1)).T
      nor=np.vstack((tangent[:,1],-tangent[:,0] )).T
      tangent=np.vstack((tangent[0,:],tangent))
      nor=np.vstack((nor[0,:],nor))
      self.tangent=np.tile(tangent,(z.shape[-1],1,1)).swapaxes(0,1)
      self.nor=np.tile(nor,(z.shape[-1],1,1)).swapaxes(0,1)	  
    def proj_hvel_along_accros(self,transect,hvel):
      """ Projet hvel of transect to along and across direction """	
      ualong=(hvel*transect.tangent[:,:,:]).sum(axis=2) 	
      uacross=(hvel*transect.nor[:,:,:]).sum(axis=2) 	
      return ualong,uacross	        


	  
class bp_transect():
 
    def __init__(self,s,ds,x,y,latlon=False,nn=False):
      """ bp_transect(s,ds,x,y,lonlat=False,nn=False) s:=schism_setup(), ds:= xarray handle ,x,y coordinates
	  Constraut transect object with track tangent normal and tangents, 
	  and along across projection function. If nn = True limit input coordinates
	  to those with different nearest neighbouts, only, ie to 
	  get xarray output access as nearest neighbour for transect coordinates. """
      self.x,self.y=x,y	  
      s.init_node_tree(latlon=latlon)
      coords=list(zip(x,y))
      p=np.asarray(coords)
      if latlon==False:	  
       dd,self.nn=s.node_tree_xy.query(coords)
      else:	  
       dd,self.nn=s.node_tree_latlon.query(coords)
      
      if nn:
       self.nn,iunq=np.unique(self.nn,return_index=True)
       self.x,self.y=self.x[iunq],self.y[iunq]
       p=p[iunq,:]
       self.ds=ds.sel(nSCHISM_hgrid_node=self.nn) #nearest neighbour

      else:		
       parents,weights=s.find_parent_tri(x,y,dThresh=1000,latlon=latlon)	# get weights  
       punq_nodes=np.unique(s.nvplt[parents]) # subslect parent nodes
       self.ds=ds.sel(nSCHISM_hgrid_node=punq_nodes)
       self.parents=s.nvplt[parents]
       #aunq,ainb,bina=intersect1d(punq_nodes, gr.faces[parents], assume_unique=False, return_indices=True)
       for nr,node in enumerate(punq_nodes):
        self.parents[self.parents==node]=nr
        #self.weights=[np.tile(np.tile(weights[:,i],(self.ds.chunks['time'][0],1)),(self.ds['zcor'].shape[2],1,1)).swapaxes(0,1).swapaxes(1,2) for i in range(3)]       
       self.weights=[np.tile(np.tile(weights[:,i],(self.ds['zcor'].shape[0],1)),(self.ds['zcor'].shape[2],1,1)).swapaxes(0,1).swapaxes(1,2) for i in range(3)]
       
    def interp_to_trans(self,varname):
      """ perform barycentric interpolation along sigmalyers for variable """
      ncomponents=1+(len(self.ds[varname].shape)==4)
      if ncomponents==2:
        out=(np.zeros(( (self.ds['zcor'].shape[0],) + self.x.shape + (self.ds['zcor'].shape[-1],ncomponents) )))
        for icomp in range(ncomponents):
            out[:,:,:,icomp]=self.ds[varname][:,self.parents[:,0],:,icomp]*self.weights[0]+self.ds[varname][:,self.parents[:,1],:,icomp]*self.weights[1]+self.ds[varname][:,self.parents[:,2],:,icomp]*self.weights[2]
        return out
      else:
       return self.ds[varname][:,self.parents[:,0],:]*self.weights[0]+self.ds[varname][:,self.parents[:,1],:]*self.weights[1]+self.ds[varname][:,self.parents[:,2],:]*self.weights[2]
       
      # projection for velocties
      dl=np.sqrt((np.diff(p,axis=0)**2).sum(axis=1))
      #dl=np.sqrt(np.diff(x)**2+np.diff(y)**2)
      self.l=np.hstack((0,np.cumsum(dl)))
      z=ds['zcor']
      self.l=np.tile(self.l,[z.shape[-1],1]).T
      tangent=np.diff(p,axis=0)/np.tile(dl,(2,1)).T
      nor=np.vstack((tangent[:,1],-tangent[:,0] )).T
      tangent=np.vstack((tangent[0,:],tangent))
      nor=np.vstack((nor[0,:],nor))
      self.tangent=np.tile(tangent,(z.shape[-1],1,1)).swapaxes(0,1)
      self.nor=np.tile(nor,(z.shape[-1],1,1)).swapaxes(0,1)	  
	  
	  
    def proj_hvel_along_accros(self,transect,hvel):
      """ Projet hvel of transect to along and across direction """	
      ualong=(hvel*transect.tangent[:,:,:]).sum(axis=2) 	
      uacross=(hvel*transect.nor[:,:,:]).sum(axis=2) 	
      return ualong,uacross	        

class bp_transect_newio():
 
    def __init__(self,s,acc,x,y,latlon=False,nn=False):
      """ bp_transect(s,ds,x,y,lonlat=False,nn=False) s:=schism_setup(), ds:= xarray handle ,x,y coordinates
	  Constraut transect object with track tangent normal and tangents, 
	  and along across projection function. If nn = True limit input coordinates
	  to those with different nearest neighbouts, only, ie to 
	  get xarray output access as nearest neighbour for transect coordinates. """
      #from IPython import embed; embed()
      self.x,self.y=x,y	  
      s.init_node_tree(latlon=latlon)
      coords=list(zip(x,y))
      p=np.asarray(coords)
      if latlon==False:	  
       dd,self.nn=s.node_tree_xy.query(coords)
      else:	  
       dd,self.nn=s.node_tree_latlon.query(coords)

      # projection for velocties
      dl=np.sqrt((np.diff(p,axis=0)**2).sum(axis=1))
      #dl=np.sqrt(np.diff(x)**2+np.diff(y)**2)
      self.l=np.hstack((0,np.cumsum(dl)))
      #z=ds['zcor']
	   
	   
	   
       # restrict to on transect data
      class cacc():
       def __init__(self,acc):
        self.vardict=acc.vardict
        self.get=acc.get

      self.acc=cacc(acc) 
      self.acc.ds={}
      if nn:
       self.nn,iunq=np.unique(self.nn,return_index=True)
       self.x,self.y=self.x[iunq],self.y[iunq]
       p=p[iunq,:]
	   
       self.acc.ds={}
       for key in np.unique(list(acc.vardict.values())):
        if type(acc.ds[key])==dict:
         self.acc.ds[key]=acc.ds[key][key].sel(nSCHISM_hgrid_node=self.nn) #nearest neighbour			
       else:
         self.acc.ds[key]=acc.ds[key].sel(nSCHISM_hgrid_node=self.nn) #nearest neighbour

      else:		
       parents,weights=s.find_parent_tri(x,y,dThresh=1000,latlon=latlon)	# get weights  
       punq_nodes=np.unique(s.nvplt[parents]) # subslect parent nodes
       self.parents=s.nvplt[parents]
       for key in np.unique(list(acc.vardict.values())):
         if type(acc.ds[key])==dict:
           self.acc.ds[key]=acc.ds[key][key].sel(nSCHISM_hgrid_node=punq_nodes) #nearest neighbour			
         else:
           self.acc.ds[key]=acc.ds[key].sel(nSCHISM_hgrid_node=punq_nodes) #nearest neighbour	   
	   
       for nr,node in enumerate(punq_nodes):
         self.parents[self.parents==node]=nr
        #self.weights=[np.tile(np.tile(weights[:,i],(self.ds.chunks['time'][0],1)),(self.ds['zcor'].shape[2],1,1)).swapaxes(0,1).swapaxes(1,2) for i in range(3)]  
       dsz=self.acc.ds['zCoordinates']['zCoordinates']		
       self.weights=[np.tile(np.tile(weights[:,i],(dsz.shape[0],1)),(dsz.shape[2],1,1)).swapaxes(0,1).swapaxes(1,2) for i in range(3)]
       
	   
      z=self.acc.ds['zCoordinates']['zCoordinates']		
      self.l=np.tile(self.l,[z.shape[-1],1]).T
      tangent=np.diff(p,axis=0)/np.tile(dl,(2,1)).T
      nor=np.vstack((tangent[:,1],-tangent[:,0] )).T
      tangent=np.vstack((tangent[0,:],tangent))
      nor=np.vstack((nor[0,:],nor))
      self.tangent=np.tile(tangent,(z.shape[-1],1,1)).swapaxes(0,1)
      self.nor=np.tile(nor,(z.shape[-1],1,1)).swapaxes(0,1)	  	   
	   
	   
    def interp_to_trans(self,varname):
      """ perform barycentric interpolation along sigmalyers for variable within xarray """
      ncomponents=1+(len(self.acc.ds[varname][varname].shape)==4)
      dsz=self.acc.ds['zCoordinates']['zCoordinates']
      dsv=self.acc.ds[varname][varname]
      if ncomponents==2:
        out=(np.zeros(( (dsz.shape[0],) + self.x.shape + (dsz.shape[-1],ncomponents) )))
        for icomp in range(ncomponents):
            out[:,:,:,icomp]=dsv[:,self.parents[:,0],:,icomp]*self.weights[0]+dsv[:,self.parents[:,1],:,icomp]*self.weights[1]+dsv[:,self.parents[:,2],:,icomp]*self.weights[2]
        return out
      else:
       return  dsv[:,self.parents[:,0],:]*self.weights[0]+ dsv[:,self.parents[:,1],:]*self.weights[1]+ dsv[:,self.parents[:,2],:]*self.weights[2]

    def interp_to_trans2(self,varname,ti=0):
      """ perform barycentric interpolation along sigmalyers for variable for one timestep """
      #from IPython import embed; embed()
      ncomponents=1+(len(self.acc.ds[varname][varname].shape)==4)
      dsz=self.acc.ds['zCoordinates']['zCoordinates']
      dsv=self.acc.ds[varname][varname]
      v=dsv[ti,:].values
      if ncomponents==2:
        out=(np.zeros(( self.x.shape + (dsz.shape[-1],ncomponents) )))
        for icomp in range(ncomponents):
            vi=dsv[self.parents[:,0],icomp].values
            out[:,:,icomp]=v[self.parents[:,0],:,icomp]*self.weights[0][0,:]+v[self.parents[:,1],:,icomp]*self.weights[1][0,:]+v[self.parents[:,2],:,icomp]*self.weights[2][0,:]
        return out
      else:
       return  v[self.parents[:,0],:]*self.weights[0][0,:]+ v[self.parents[:,1],:]*self.weights[1][0,:]+ v[self.parents[:,2],:]*self.weights[2][0,:]	   
      
    def proj_hvel_along_accros(self,transect,hvel):
      """ Projet hvel of transect to along and across direction """	
      ualong=(hvel*transect.tangent[:,:,:]).sum(axis=2) 	
      uacross=(hvel*transect.nor[:,:,:]).sum(axis=2) 	
      return ualong,uacross	        

	  
	  
class sources:		
	"""	load source fluxes for schism setup """
	
	def __init__(self,s,fname='source_sink.in',comments='!'):
		self.i_src=np.asarray(np.loadtxt(fname,comments=comments)[1:-1],int)-1
		#self.elem=np.asarray(s.nv)[self.i_src,:]-1	
		#self.x=np.asarray(s.x)[self.elem].mean(axis=1)
		#self.y=np.asarray(s.y)[self.elem].mean(axis=1)
		# also for quads:
		self.x=[]
		self.y=[]
		x=np.asarray(s.x)
		y=np.asarray(s.y)
		for ielem in self.i_src:
			elem=np.asarray(s.nv[ielem])
			self.x.append(x[elem-1].mean())
			self.y.append(y[elem-1].mean())
		self.x=np.asarray(self.x)	
		self.y=np.asarray(self.y)	
		
		
		with open(fname) as f:
			n=int(f.readline().split(comments)[0])
			#self.names=[f.readline().split(comments)[1].split(',')[0].split('river')[1] for i in range(n)]
			self.names=[f.readline().split(comments)[1] for i in range(n)]
			
		# volume
		try:
			self.vsource=np.asarray(np.loadtxt('vsource.th',comments=comments),int)
			self.vtime=self.vsource[:,0]
			self.vsource=self.vsource[:,1:]
			
			# concentration i.e temperature
			self.msource=np.asarray(np.loadtxt('msource.th',comments=comments),int)
			self.mtime=self.msource[:,0]
			self.msource=self.msource[:,1:]	
		except:
			pass
	def plot_sources(self):
		names=[name.split(',')[0] for name in self.names]
		nsrc=int(self.vsource.shape[-1])
		plt.subplot(2,2,1)
		plt.plot(self.vtime/86400,self.vsource)
		plt.ylabel('Q [m^3/s]')
		plt.subplot(2,2,2)
		plt.plot(self.vtime/86400,self.msource[:,:nsrc])
		plt.ylabel('T [deg C]')			
		plt.subplot(2,2,3)
		plt.plot(self.vtime/86400,self.msource[:,nsrc:])
		plt.ylabel('S [psu]')
		plt.subplot(2,2,4)
		plt.plot(self.vtime/86400,self.vsource*np.nan)
		plt.legend(names,ncol=2)
		plt.tight_layout()
		plt.show()


class schism_station_output:
	def __init__(self,fname='station.in',output='outputs',comments='!'):
		with open('station.in') as f:
			#line=f.readline().replace('\n','')
			lines=f.readlines()
			active_flag,var_names=lines[0].replace('\n','').split('for')
			active_flag=np.asarray([int(val) for val in active_flag.split('!')[0].split(' ')[:9]])
			var_names=[name.replace(' ','') for name in var_names.split(',')]
			nstations=int(lines[1].split('!')[0])
			
			line=lines[2]
			num,txt=line.split('!')[:2]
			self.coords=np.asarray([np.float(i) for i in num.split(' ')[1:4]])
			self.stations=[txt.split(' ')[1]]
			for line in lines[3:]:
				num,txt=line.replace('\n','').split('!')
				self.coords=np.vstack((self.coords,np.asarray([np.float(i) for i in num.split(' ')[1:4]])))
				self.stations.append(txt[1:])
			
			# coords to lon lat based on nearset neighbour
			if (self.coords[:,1]>90).sum() > 0:
				s=schism_setup()
				s.init_node_tree(latlon=False)
				nn=s.node_tree_xy.query(list(zip(self.coords[:,0],self.coords[:,1])))[1]
				lon,lat=np.asarray(s.lon),np.asarray(s.lat)	
				self.coords[:,0]=lon[nn]
				self.coords[:,1]=lat[nn]
				
			self.station_out={name:0 for name in var_names}
			for nr in np.where(active_flag)[0]:
				values=np.loadtxt(output+'/staout_{:d}'.format(nr+1))
				self.station_out[var_names[nr]]=values[:,1:]
			self.station_out['time']=values[:,0]
			
			with open('param.nml') as f:
				for line in f.readlines():
					if 'start_year' in line:
						year=int(line.split('=')[1].split('!')[0])
					elif 'start_month' in line:
						month=int(line.split('=')[1].split('!')[0])
					elif 'start_day' in line:
						day=int(line.split('=')[1].split('!')[0])
					elif 'start_hour' in line:
						hour=np.float(line.split('=')[1].split('!')[0])
						minute=int((hour-np.floor(hour))*60)
						hour=int(np.floor(hour))
						break
			refdate=dt.datetime(year,month,day,hour,minute)
			self.time=refdate+dt.timedelta(seconds=self.station_out['time'][0])*np.arange(1,len(self.station_out['time'])+1)
			
	def plot_station(self,istat=0):		
		if type(istat)==int:
			i=0				
		else: # get from name
			for i,name in enumerate(stations):
				if istat.lower() in name.lower():
					break
			
		plt.plot(self.time,
		self.station_out['elev'][:,i])
		plt.title(self.stations[i])
		plt.grid()
		plt.show()
		plt.gcf().autofmt_xdate()
		plt.tight_layout()

class schism_outputs_by_variable():
   def __init__(self,ncdir='./outputs/',max_stack=-1,varlist=None):
      """ output class for xarray access to schism  by variable output
			schism_outputs_by_variable(ncdir='./outputs/',max_stack=-1). ncdir is netcdf output directory max_stack is the highest number of stacks(this is identical with the highest stack_number only if all stacks from _1 in the folder). If varfiles==None, If varfiles=[out2d] only those files with a matching pattern will be taken"""
      self.ncdir=ncdir

      #from IPython import embed; embed()
	  
	  
      nr_nc0=np.sort(glob.glob('{:s}/out2d*.nc'.format(ncdir)))[0].split('/')[-1].split('_')[-1]
      
      varfiles=glob.glob('{:s}*_{:s}'.format(ncdir,nr_nc0))
      if varlist!=None and type(varlist)==list:   
         varfiles=[var[var.rindex('/')+1:var.rindex('_')] for var in varfiles if var[var.rindex('/')+1:var.rindex('_')] in varlist]
      else:
         varfiles=[var[var.rindex('/')+1:var.rindex('_')] for var in varfiles]
      
      files=dict.fromkeys(varfiles)
      self.ds=dict.fromkeys(varfiles)
      for key in files.keys():
         files[key]=np.hstack([np.sort(glob.glob('{:s}{:s}_{:s}.nc'.format(ncdir,key,'?'*iorder))) for iorder in range(1,6)])
         #self.ds[key]=xr.open_mfdataset(files[key][:max_stack]) # this crashes
         #self.ds[key]= [xr.open_dataset(file).chunk() for file in files[key][:max_stack]]
		 self.ds[key]=xr.concat([xr.open_dataset(file).chunk() for file in files[key][:max_stack]],dim='time')
		#s.nc=xr.concat(all_dsets, dim='time')

      exclude=[]#exclude=['time','SCHISM_hgrid', 'SCHISM_hgrid_face_nodes', 'SCHISM_hgrid_edge_nodes', 'SCHISM_hgrid_node_x',     'SCHISM_hgrid_node_y', 'bottom_index_node', 'SCHISM_hgrid_face_x', 'SCHISM_hgrid_face_y',     'ele_bottom_index', 'SCHISM_hgrid_edge_x', 'SCHISM_hgrid_edge_y', 'edge_bottom_index',     'sigma', 'dry_value_flag', 'coordinate_system_flag', 'minimum_depth', 'sigma_h_c', 'sigma_theta_b',      'sigma_theta_f', 'sigma_maxdepth', 'Cs'] # exclude for plot selection
      self.vardict={} # variable to nc dict relations
      self.varlist=list(self.vardict.keys())	  
      vector_vars=[] # stack components for convenience
      for nci_key in self.ds.keys():
          for vari in self.ds[nci_key].keys():
              if vari not in exclude:
                  self.vardict[vari]=nci_key	
                  self.varlist+=[nci_key]
	  
      # merge vector components
              if vari[-1] =='Y': 
                  vector_vars.append(vari[:-1])
                  varY=vari
                  vari_vec=vari[:-1]
                  varX=vari_vec+'X'
                  self.varlist+=[vari_vec]
                  self.vardict[vari_vec]=vari_vec
                  #self.ds[vari_vec] ={vari_vec: xr.concat([self.ds[self.vardict[varX]][varX], self.ds[self.vardict[varY]][varY]], dim='ivs')}

      for vari_vec in vector_vars:			
          varX=vari_vec+'X'	  
          varY=vari_vec+'Y'	  
          self.ds[vari_vec] ={vari_vec: xr.concat([self.ds[self.vardict[varX]][varX], self.ds[self.vardict[varY]][varY]], dim='ivs')}
				  
      #self.varlist=list(self.vardict.keys())	  
      # add velocity vecotr # moved to above part
      #if	'horizontalVelX' in self.varlist:
      #    
      #    self.varlist+=['hvel']
      #    self.vardict['hvel']='hvel'
      #    self.ds['hvel'] ={'hvel': xr.concat([self.ds[self.vardict['horizontalVelX']]['horizontalVelX'], self.ds[self.vardict['horizontalVelY']]['horizontalVelY']], dim='ivs')}
      #	
      #if	'windSpeedX' in self.varlist:			
      #    self.varlist+=['wind']
      #    self.vardict['wind']='wind'
      #    self.ds['wind'] ={'wind': xr.concat([self.ds[self.vardict['windSpeedX']]['windSpeedX'], self.ds[self.vardict['windSpeedY']]['windSpeedY']], dim='ivs')}
      self.varlist=list(np.unique(np.sort(self.varlist)))			

   def get(self,varname):
      """ Return handle to xarray dataset"""	
      return self.ds[self.vardict[varname]][varname]	  
		
class param:		
	"""	functions for param.in for reading and editing. Operates in local directory """
	import os
	def __init__(self,fname='param.nml',comments='!'):
		#self.param_in=np.asarray(np.loadtxt(fname,comments=comments)[1:-1],int)-1		
		
		if '/' in fname:
			islash=fname.rindex('/')
			self.dir=fname[:islash+1]		
		else:
			self.dir='./'
			
		f=open(self.dir+'param.nml')	
		self.lines=f.readlines()
		f.close()
		
	def get_parameter(self,param='dt'):
		""" read parameter from param.in"""
		
		for line in self.lines:
			if param+' =' in line:
				param= line.split('=')[1].split('!')[0]
				try:
					param=float(param)
				except:
					param=str(param)
				break
		return param

	def set_parameter(self,params,values,outname='param.nml',outdir='./'):
		"""set_parameter(self,params,values,outname='param.nml',outdir='./') change parameters in param.in """
		if outname=='param.nml':
			try:
				os.rename('param.nml','param.nml.bkp')
			except:
				pass
		
		if type(params) == str:
			params=[params,]
			values=[values,]
		fout=open(outdir+outname,'w') 
		for line in self.lines:
			for param,value in zip(params,values):
				if param+' =' in line:
					line=' {:s} = {:.0f} !'.format(param,value)+line.split('!')[1]+'\n'
					values.remove(value)
					params.remove(param)
			fout.write(line)		

		fout.close()		
		# load updated param.nml
		print('updated param.nml has been loaded and will be accessed by get_parameters')	
		f=open(outdir+outname)	
		self.lines=f.readlines()
		f.close()
			
	def set_time_step(self,dt,outname='param.nml',outdir='./',rnday=None):
		""" set_time_step(self,dt,outname='param.nml',outdir='./',rnday=None) updates time step (dt) and all time related parameters (nspool,ihfskip,nhot_write,nspool_sta) to maintain output timings. rnday != None changes simulate length to specified Nr of days. dt new time step """
		params=['dt','nspool','ihfskip','nhot_write','nspool_sta','wtiminc']
		values=[self.get_parameter(param=param) for param in params]
		values=[dt]+list(np.asarray(values[1:])*values[0]/dt)
		if rnday != None:
				params.append('rnday')
				values.append(rnday)
		self.set_parameter(params=params,values=values,outname=outname,outdir='./')		
		
		
		
if __name__ == '__main__':

    from pylab import *
    setup = schism_setup('hgrid.gr3',ll_file='hgrid.gr3')
    # plot domain
    setup.plot_domain_boundaries()

    #triplot(setup.x,setup.y,asarray(setup.nv)-1)

    show()

    if False:
      # read elevation boundaries
      t,e = setup.bdy_array('elev2D.th')
      figure()
      plot(t[:],e[:,50])
      show()
