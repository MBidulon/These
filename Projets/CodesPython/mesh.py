import inspect;
import numpy as np;

class UnionFind:
    """A UnionFind structure for computing connected components."""
    def __init__(self,label):
        """Initialize a node with a given label."""
        self.parent=self;
        self.label=label;
        self.rank=0;

    def Find(self):
        """Return the representant ancestor."""
        if self.parent == self:
            return self;
        else:
            self.parent=self.parent.Find();
            return self.parent;

    def Union(self,y):
        """Merge oneself with another element."""
        selfRoot=self.Find();
        yRoot=y.Find();
        if selfRoot.rank > yRoot.rank:
            yRoot.parent=selfRoot;
        elif selfRoot.rank< yRoot.rank:
            selfRoot.parent=yRoot;
        elif selfRoot != yRoot:
            yRoot.parent=selfRoot;
            selfRoot.rank=selfRoot.rank+1;


class Mesh:
    """A mesh structure matching the medit format."""
    def __init__(self,meshFile):
        """Initialize a mesh object with a mesh file in the medit format."""
        self.MeshVersionFormatted=None;
        self.Dimension=None;
        self.Connectivity=dict();
        """Connectivity structure, essentially all the fields of the input .mesh."""
        self.Boundaries=dict();
        self.debug=0;
        self.meshFile=meshFile;
        self.graphVertices=None;
        """A variable containing a graph of all the vertex neighbors. Useful when computing connected components."""
        if self.debug>0:
            print("Loading "+meshFile);
        if meshFile=="":
            return;
        self.Boundaries=dict();
        """A dictionary of all the boundaries of the mesh."""
        self.Connectivity=dict();
        f=open(meshFile,"r");
        lines=[line for line in f.readlines() if line.rstrip()];
        f.close();
        self.MeshVersionFormatted=int(lines[0].split()[1]);
        try:
            self.Dimension=int(lines[1].split()[1]);
            position=2;
        except:
            self.Dimension=int(lines[2].split()[0]);
            position=3;
        startMesh=False;
        while startMesh==False:
            key=lines[position].strip();
            if key=="Vertices" or key=="END" or key=='End"':
                startMesh=True;
                position=position-1
            position=position+1
        while position<len(lines):
            key=lines[position].strip();
            if key=="END" or key=='End':
                break;
            nentries=int(lines[position+1]);
            if self.debug>0:
                print(f"Found {nentries} entries in section {key}");
            position=position+2;
            self.Connectivity[key]=[];
            if key=="Vertices":
                for i in range(nentries):
                    self.Connectivity[key].append([ float(x) for x in lines[position+i].split()[0:self.Dimension]]);
                    self.Connectivity[key][i].append(int(lines[position+i].split()[self.Dimension]));
            else:
                self.Connectivity[key]=[ [int(x) for x in lines[position+i].split()] for i in range(nentries) ];
            position=position+nentries;
        self.getBoundaries();
    
    def getEdgeLength(self,edgeNumber):
        """Return the length of the edge edgeNumber."""
        edge=self.Connectivity['Edges'][edgeNumber];
        pt1=self.Connectivity['Vertices'][edge[0]-1];
        pt2=self.Connectivity['Vertices'][edge[1]-1];
        length=np.sqrt(sum([(x-y)**2 for (x,y) in zip(pt1[0:-1],pt2[0:-1])]));
        if self.debug>10:
            print("Length of edge "+str(edgeNumber)+": \n"\
                  +"pt1="+pt1.__str__()+"; pt2="+pt2.__str__()+"\n length="+str(length)+"\n");
        return length;

    def getDistance(self,pt1,pt2):
        """Compute distance between two vertices pt1 and pt2 given by index number.
        Warning: the numbering is the one of the .mesh file (starting from 1, not the 
        one of python).
        """
        point1=self.Connectivity['Vertices'][pt1-1];
        point2=self.Connectivity['Vertices'][pt2-1];
        length=np.sqrt(sum([(x-y)**2 for (x,y) in zip(point1[0:-1],point2[0:-1])]));
        if self.debug>10:
            print("Length of edge "+str(edgeNumber)+": \n"\
                  +"pt1="+pt1.__str__()+"; pt2="+pt2.__str__()+"\n length="+str(length)+"\n");
        return length;

    def getBoundaries(self): 
        """Compute arrays containing various boundaries of the mesh."""
        if 'Edges' in self.Connectivity:
            key='Edges';
        else:
            key='Tetrahedra';
        for (k,e) in enumerate(self.Connectivity[key]):
            if not  e[-1] in self.Boundaries:
                self.Boundaries[e[-1]]=[];
                if self.debug>0:
                    print("Found boundary index "+str(e[-1]));
            self.Boundaries[e[-1]].append(k);

    def connectedComponent(self):
        """Compute the connected component of the triangles of the mesh.
        Returns a UnionFind structure whose groups are connected components of triangles sharing
        the same label."""
        nt=len(self.Connectivity['Triangles']);
        if self.graphVertices is None:
            self.graphVertices=[set() for x in range(self.nv())];
            for (k,tri) in enumerate(self.Connectivity['Triangles']):
                for pt in tri[:-1]:
                    self.graphVertices[pt-1].add(k);
        def neighbors(i):
                neigh={x for j in self.Connectivity['Triangles'][i][:-1] for x in self.graphVertices[j-1] \
                       if self.Connectivity['Triangles'][i][3]==self.Connectivity['Triangles'][x][3]};
                return neigh;
        graph=[UnionFind(i) for i in range(nt)];
        for tri in range(nt):
           try:
               neighboors=set(neighbors(tri));
           except:
               import pdb;
               pdb.set_trace();
           for n in neighboors:
               graph[tri].Union(graph[n]);
        return graph;


    def nt(self):
        """Number of triangles"""
        return len(self.Connectivity['Triangles']);

    def nv(self):
        """Number of vertices"""
        return len(self.Connectivity['Vertices']);

    def trunc(self,region):
        M=Mesh("");
        M.Dimension=self.Dimension;
        M.MeshVersionFormatted=self.MeshVersionFormatted;
        correspondences=dict();
        for key in self.Connectivity:
            correspondences[key]=dict();
        if self.Dimension==2:
            countVertices=0;
            countTriangles=0;
            fields={'Triangles','Vertices','Edges'};
            for field in fields:
                M.Connectivity[field]=[];
            for triNum,tri in enumerate(self.Connectivity['Triangles']):
                if tri[-1]==region:
                    if self.debug>3:
                        print(f"Adding triangle {tri}");
                    for pt in tri[0:3]:
                        if not pt in correspondences['Vertices']:
                            correspondences['Vertices'][pt]=countVertices+1;
                            M.Connectivity['Vertices'].append(self.Connectivity['Vertices'][pt-1]);
                            countVertices+=1;
                            if self.debug>3:
                                print(f"Added vertex {pt} ==> {countVertices}.");
                    M.Connectivity['Triangles'].append([correspondences['Vertices'][tri[0]],correspondences['Vertices'][tri[1]],\
                        correspondences['Vertices'][tri[2]],region]);
                    if self.debug>3:
                        print(f"Added triangle {triNum+1} ==> {len(M.Connectivity['Triangles'])}");
            if self.debug>2:
                print(f"Added {countTriangles} triangles.");
                print(f"Added {countVertices} vertices.");
            counterEdges=0;
            for edgeNum,edge in enumerate(self.Connectivity['Edges']):
                if edge[0] in correspondences['Vertices'] and edge[1] in correspondences['Vertices']:
                    counterEdges+=1;
                    M.Connectivity['Edges'].append([correspondences['Vertices'][edge[0]],correspondences['Vertices'][edge[1]],\
                            edge[2]]);
                    if self.debug>3:
                        print(f"Adding boundary edge {edgeNum+1}==> {counterEdges}");
            if self.debug>2:
                print(f"Added {counterEdges} edges.");
        return M;

    def save(self,meshFile):
        """Save the mesh in the medit file format."""
        f=open(meshFile,"w");
        f.write("MeshVersionFormatted "+str(self.MeshVersionFormatted)+"\n");
        f.write("\n\n");
        f.write("Dimension "+str(self.Dimension)+"\n\n\n");
        f.write("Vertices\n");
        f.write(str(len(self.Connectivity['Vertices']))+"\n");
        for elem in self.Connectivity["Vertices"]:
            f.write(" ".join(map(str,elem))+"\n");
        for (key,values) in self.Connectivity.items():
            if key!="Vertices":
                f.write("\n\n");
                f.write(key+"\n");
                f.write(str(len(values))+"\n");
                for elem in values:
                    f.write(" ".join(map(str,elem))+"\n");
        f.write("\n\n");
        f.write("END\n");
        f.close(); 
        self.meshFile=meshFile;
    
    def saveFormat2(self,meshFile):
        """Save the mesh in the medit file format."""
        f=open(meshFile,"w");
        f.write("MeshVersionFormatted "+str(self.MeshVersionFormatted)+"\n");
        f.write("\n\n");
        f.write("Dimension "+str(self.Dimension)+"\n\n\n");
        f.write("Vertices\n");
        f.write(str(len(self.Connectivity['Vertices']))+"\n");
        for elem in self.Connectivity["Vertices"]:
            f.write(" ".join(map(str,elem))+"\n");
        f.write("\n\n");
        f.write("Triangles\n");
        f.write(str(len(self.Connectivity['Triangles']))+"\n");
        for elem in self.Connectivity["Triangles"]:
            f.write(" ".join(map(str,elem))+"\n");
        f.write("\n\n");
        f.write("END\n");
        f.close();
        self.meshFile=meshFile;
                    
    def plot(self,colormap=None,XLIM=None,YLIM=None,figsizeX=None,figsizeY=None,resolution=None,doNotPlot=False,edgeColor='gray',axis='off',\
            boundary=None,linewidth=0.3):
        """Plot the mesh with matplotlib library."""
        import matplotlib as mp;
        import matplotlib.pyplot as plt
        coords=list(map(list,zip(*self.Connectivity['Vertices'])));
        x=np.asarray(coords[0]);
        y=np.asarray(coords[1]);
        tri=np.asarray([[x-1 for x in t[:-1]] for t in self.Connectivity['Triangles']]);
        import random;
        colorIndex=lambda x: random.randint(0,10);
        colors=np.asarray([colorIndex(t[-1]) for t in self.Connectivity['Triangles']]);
        if colormap=='bw':
            colormap=[[1,1,1],[0,0,0]];
            edgeColor='none';
        if colormap=='dim':
            colormap=[[0.9,0.9,1],[1,1,0.9]];
            edgeColor='gray';
        if colormap is None:
            mapp=np.array([[0.5,1,0.5],[1,1,0.5]]);
            cmap=mp.colors.ListedColormap(mapp);
        try:
            cmap=mp.colors.ListedColormap(colormap);
        except:
            pass;
        if figsizeX==None:
            if resolution==None:
                fig, ax = plt.subplots()
            else:
                fig, ax = plt.subplots(dpi=resolution)
        else:
            if resolution==None:
                fig, ax = plt.subplots(figsize=(figsizeX,figsizeY))
            else:
                fig, ax = plt.subplots(figsize=(figsizeX,figsizeY),dpi=resolution)
        ax.tripcolor(x,y,tri,facecolors=colors,edgecolors=edgeColor,cmap=cmap,linewidth=linewidth);
        ax.margins(0);
        #fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
        ax.set_aspect('equal');
        if axis=='off':
            ax.tick_params(axis='both',which='both',length=0);
            plt.setp(ax.get_xticklabels(), visible=False)
            plt.setp(ax.get_yticklabels(), visible=False)
        if not XLIM is None:
            ax.set_xlim(XLIM);
        if not YLIM is None:
            ax.set_ylim(YLIM);
        lines=[];
        if boundary:
            for e in self.Boundaries[10]:
                pts=self.Connectivity['Edges'][e][:-1];
                pt1=self.Connectivity['Vertices'][pts[0]-1];
                pt2=self.Connectivity['Vertices'][pts[1]-1];
                lines.append([(pt1[0],pt1[1]),(pt2[0],pt2[1])]);
            lc = mp.collections.LineCollection(lines,linewidths=0.5,colors=boundary)
            ax.add_collection(lc);
        fig.tight_layout();
        ax.autoscale_view();
        if not doNotPlot:
            plt.show();
        return fig,ax;
#        return (fig,ax);

class P1Function:
    """A structure for P1 function defined on a medit mesh, based on the medit .sol format."""

    def __init__(self,m,phi=None):
        """Initialize the P1 function. Arguments are:
        m         :  input mesh
        phi       :  Either: 
                        - a .sol file of a scalar function
                        - a list or a numpy.ndarray of function values for each of the mesh
                          vertices with the order given by m.Connectivity['Vertices']
                        - a lambda function x : f(x). The values at each node of coordinates
                            x[0],x[1] will be computed accordingly.
        """
        
        self.mesh=m;
        self.MeshVersionFormatted=self.mesh.MeshVersionFormatted;
        self.Dimension=self.mesh.Dimension;
        self.values=[];
        self.isSolid=[0 for i in range(self.mesh.nt())];
        if isinstance(phi,np.ndarray):
            phi=list(phi);
        if isinstance(phi,list):
            if len(phi)!=self.mesh.nv():
                raise Exception("Error : the number of function values \
                                "+str(len(phi))+"!="+str(self.mesh.nv())+"number of \
                                mesh nodes");
            self.values=phi;
        if phi is None:
            self.values=[0 for i in range(self.mesh.nv())];
        if isinstance(phi,str):
            f=open(phi,"r");
            lines=[line for line in f.readlines() if line.rstrip()];
            f.close();
            if phi.endswith('.sol'):
                self.MeshVersionFormatted=int(lines[0].split()[1]);
                self.Dimension=int(lines[1].split()[1]);
                nv=int(lines[3]);
                self.values=[float(x) for x in lines[5:5+nv]];
            elif phi.endswith('.gp'):
                self.MeshVersionFormatted=1;
                nv=int(lines[0]);
                self.values=[float(x) for line in lines[1:] for x in line.split()];
            if len(self.values)!=self.mesh.nv() or nv!=self.mesh.nv():
                raise Exception("Error : the number of function values "+\
                        str(nv)+"!="+str(self.mesh.nv())+" number of mesh nodes");
        if inspect.isfunction(phi):
            for pt in self.mesh.Connectivity["Vertices"]:
                self.values.append(phi(pt));
        self.values=np.asarray(self.values);
        self.getIsSolid();


    def getIsSolid (self):
        print("in the function")
        nt=len(self.mesh.Connectivity['Triangles']);
        for tri in range(nt):
            pt1=self.mesh.Connectivity['Triangles'][tri][0];
            pt2=self.mesh.Connectivity['Triangles'][tri][1];
            pt3=self.mesh.Connectivity['Triangles'][tri][2];
            if ((self.values[pt1-1]<0) | (self.values[pt2-1]<0) | (self.values[pt3-1]<0)):
                self.isSolid[tri]=1;
    
    def testIsSolid(self,cmap='jet',figsizeX=None, figsizeY=None, resolution=None, doNotPlot=False,tickFormat='2.1f',bcColor='k',bcLineWidth=0.5,niso=49,XLIM=None,YLIM=None):
        """Plot a P1 function with matplotlib."""
        import matplotlib as mp;
        import numpy as np
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        coords=list(map(list,zip(*self.mesh.Connectivity['Vertices'])));
        x=np.asarray(coords[0]);
        y=np.asarray(coords[1]);
        triang=np.asarray([[x-1 for x in t[:-1]] for t in self.mesh.Connectivity['Triangles']]);

        z=[0 for i in range(self.mesh.nv())];
        nt=len(self.mesh.Connectivity['Triangles']);
        for tri in range(nt):
            for pt in range(3):
                pt2=self.mesh.Connectivity['Triangles'][tri][pt];
                if (z[pt2-1]==0):
                    z[pt2-1]=self.isSolid[tri];
        fig, ax = plt.subplots();
        plot=plt.tricontourf(x,y,triang,z,niso,cmap='jet');
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        ax.margins(0);
        cbar=fig.colorbar(plot,cax=cax);
        cbar.set_ticks(np.linspace(min(z),max(z),5));
        cbar.set_ticklabels([format(x,tickFormat) for x in np.linspace(min(z),max(z),5)]);
        ax.set_aspect('equal');
        ax.tick_params(axis='both',which='both',length=0);
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        if not XLIM is None:
            ax.set_xlim(XLIM);
        if not YLIM is None:
            ax.set_ylim(YLIM);
        if not doNotPlot:
            plt.show();


    def connectedComponent(self):
        """Compute the connected component of the triangles of the mesh.
        Returns a UnionFind structure whose groups are connected components of triangles sharing
        the same label."""
        nt=len(self.mesh.Connectivity['Triangles']);
        if self.mesh.graphVertices is None:
            self.mesh.graphVertices=[set() for x in range(self.mesh.nv())];
            for (k,tri) in enumerate(self.mesh.Connectivity['Triangles']):
                for pt in tri[:-1]:
                    self.mesh.graphVertices[pt-1].add(k);
        def neighbors(i):
                neigh=set();
                for j in self.mesh.Connectivity['Triangles'][i][:-1]: #j correspond au numero de point freefem
                    for x in self.mesh.graphVertices[j-1]: #x correspond au numero de triangle python
                        if self.isSolid[i]==self.isSolid[x] :
                            neigh.add(x);
                return neigh
        graph=[UnionFind(i) for i in range(nt)];
        for tri in range(nt):
            graph[tri].label=-1;
        compConInt=0;
        compConExt=0;
        for tri in range(nt):
            try:
                neighboors=set(neighbors(tri));
            except:
                import pdb;
                pdb.set_trace();
            for n in neighboors:
                graph[tri].Union(graph[n]);
        for tri in range(nt):
            parent=graph[tri].Find();
            if(parent.label==-1):
                parent.label=self.isSolid[tri]*(2*compConInt)+(1-self.isSolid[tri])*(2*compConExt+1);
                compConInt+=self.isSolid[tri];
                compConExt+=(1-self.isSolid[tri]);
            graph[tri].label=parent.label;
        for tri in range(nt):
            self.mesh.Connectivity['Triangles'][tri][-1]=graph[tri].label;
        print("nb comp intérieures = ",compConInt);
        print("nb comp extérieures = ",compConExt);
        a=[0,0];
        a[0]=compConInt;
        a[1]=compConExt;
        return(a);


    def plot(self,cmap='jet',figsizeX=None, figsizeY=None, resolution=None, doNotPlot=False,tickFormat='2.1f',bcColor='k',bcLineWidth=0.5,niso=49,XLIM=None,YLIM=None):
        """Plot a P1 function with matplotlib."""
        import matplotlib as mp;
        import numpy as np
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        coords=list(map(list,zip(*self.mesh.Connectivity['Vertices'])));
        x=np.asarray(coords[0]);
        y=np.asarray(coords[1]);
        triang=np.asarray([[x-1 for x in t[:-1]] for t in self.mesh.Connectivity['Triangles']]);
        z=self.values
        
        vMax=z.max(axis=0);
        vMin=z.min(axis=0);
        vMax=(vMax-vMin)/18.;
        visoV2=np.array([-1,0,vMin+vMax,vMin+2*vMax,vMin+3*vMax,vMin+4*vMax,vMin+5*vMax,vMin+6*vMax,vMin+7*vMax,vMin+8*vMax,vMin+9*vMax,vMin+10*vMax,vMin+11*vMax,vMin+12*vMax,vMin+13*vMax,vMin+14*vMax,vMin+15*vMax,vMin+16*vMax,vMin+17*vMax,vMin+18*vMax,vMin+19*vMax]);
        
        if figsizeX==None:
            if resolution==None:
                fig, ax = plt.subplots()
            else:
                fig, ax = plt.subplots(dpi=resolution)
        else:
            if resolution==None:
                fig, ax = plt.subplots(figsize=(figsizeX,figsizeY))
            else:
                fig, ax = plt.subplots(figsize=(figsizeX,figsizeY),dpi=resolution)
        
#        plot=plt.tricontour(x,y,triang,z,visoV2,niso,cmap='jet');
        plot=plt.tricontour(x,y,triang,z,niso,cmap='jet');
        #plot=ax.tripcolor(x,y,z,cmap=cmap,shading='gouraud');
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        ax.margins(0);
        #cbar=fig.colorbar(plot,cax=cax);
        #cbar.set_ticks(np.linspace(min(z),max(z),5));
        #cbar.set_ticklabels([format(x,tickFormat) for x in np.linspace(min(z),max(z),5)]);
        ax.set_aspect('equal');
        ax.tick_params(axis='both',which='both',length=0);
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        if not XLIM is None:
            ax.set_xlim(XLIM);
        if not YLIM is None:
            ax.set_ylim(YLIM);
        if not doNotPlot:
            plt.show();

    def save(self,fileName):
       """Save in the medit file format .sol."""
       f=open(fileName,"w");
       f.write("MeshVersionFormatted "+str(self.MeshVersionFormatted)+"\n\n");
       f.write("Dimension "+str(self.Dimension)+"\n\n");
       f.write("SolAtVertices\n");
       f.write(str(self.mesh.nv())+"\n");
       f.write("1 1\n");
       for val in self.values:
           f.write(str(val)+"\n");
       f.write("\n\nEnd"); 
       

        

class P1Vector:
    """A structure for P1 vector functions on medit mesh, based on the medit .sol format."""
    mesh=None;
    values=[];
    MeshVersionFormatted=None;
    Dimension=None;

    def __init__(self,m,phi):
        """Create a P1Vector function from an array."""
        self.mesh=m;
        self.Dimension=m.Dimension;
        self.MeshVersionFormatted=m.MeshVersionFormatted;
        try:
            if isinstance(phi[0],str):
                self.values=[];
                for i in range(len(phi)):
                    self.values.append(P1Function(m,phi[i]).values);
                self.values=np.asarray(self.values).transpose();
            else:
                self.values=np.array(phi).transpose();
        except:
            raise Exception(f"Error: unknown vector entry {phi}");

    def save(self,fileName):
       """Save in the medit .sol file format."""
       f=open(fileName,"w");
       f.write("MeshVersionFormatted "+str(self.MeshVersionFormatted)+"\n\n");
       f.write("Dimension "+str(self.Dimension)+"\n\n");
       f.write("SolAtVertices\n");
       f.write(str(self.mesh.nv())+"\n");
       f.write("1 2\n");
       for val in self.values:
           f.write(" ".join([str(x) for x in val])+"\n");
       f.write("\n\nEnd"); 


