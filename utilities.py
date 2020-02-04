#----------------------------------------------------------------------
# Copyright (C) 2017, All rights reserved
#
# JENSEN HUGHES
#
# 3610 Commerce Blvd, Suite 817
#
# Baltimore, MD 21227
#
# http://www.jensenhughes.com
#
# JENSEN HUGHES. Copyright Information
#
#----------------------------------------------------------------------
#======================================================================
# DESCRIPTION:
# This script has utilities used to visualize polygons from a series of points.
#
# ATTRIBUTION:
# The following subroutines were based on work posted to stackoverflow:
#   pnt_in_cvex_hull
#   simplify
#   simplifyFaces
#   is_adjacent
#   is_coplanar
#   reorder
#   get_distance
# The original question was posted here:
#   https://stackoverflow.com/questions/49098466/plot-3d-convex-closed-regions-in-matplot-lib
# The original author was Nikferrari:
#   https://stackoverflow.com/users/7122932/nikferrari
# The answerer was Paul Brodersen:
#   https://stackoverflow.com/users/2912349/paul-brodersen

import networkx as nx
import scipy.spatial as scsp
import numpy as np
import mpl_toolkits.mplot3d as a3
import matplotlib.pyplot as plt
import matplotlib.colors as pltc

def kalmanFilter(z, Q=1e-5):

    # intial parameters
    #n_iter = 50
    sz = z.shape[0] # size of array
    #x = -0.37727 # truth value (typo in example at top of p. 13 calls this z)
    #z = np.random.normal(x,0.1,size=sz) # observations (normal about x, sigma=0.1)

    #Q = 1e-5 # process variance

    # allocate space for arrays
    xhat=np.zeros(sz)      # a posteri estimate of x
    P=np.zeros(sz)         # a posteri error estimate
    xhatminus=np.zeros(sz) # a priori estimate of x
    Pminus=np.zeros(sz)    # a priori error estimate
    K=np.zeros(sz)         # gain or blending factor
    
    R = 0.5**2 # estimate of measurement variance, change to see effect
    
    # intial guesses
    xhat[0] = z[0]
    P[0] = 1.0
    
    for k in range(1,sz):
        # time update
        xhatminus[k] = xhat[k-1]
        Pminus[k] = P[k-1]+Q
    
        # measurement update
        K[k] = Pminus[k]/( Pminus[k]+R )
        xhat[k] = xhatminus[k]+K[k]*(z[k]-xhatminus[k])
        P[k] = (1-K[k])*Pminus[k]
    
    return xhat

def pts2polygons(groups):
    '''
    Build polygons from series of points.
    '''
    polygons = []
    
    for group in groups:
        linkedPolygons = []
        for pts in group:
            try:
                linkedPolygons.append(scsp.ConvexHull(pts))
            except:
                print("Failed points:")
                print(pts)
        polygons.append(linkedPolygons)
        
    return polygons, len(polygons)

def ConvexHull(pts):
    p = scsp.ConvexHull(pts)
    return p

def pnt_in_cvex_hull(hull, pnt):
    '''
    This subroutine checks to see if a point is contained within a polygon.
    '''
    new_hull = scsp.ConvexHull(np.concatenate((hull.points, [pnt])))
    if np.array_equal(new_hull.vertices, hull.vertices): 
        return True
    return False

def simplify(triangles):
    """
    Simplify an iterable of triangles such that adjacent and coplanar triangles form a single face.
    Each triangle is a set of 3 points in 3D space.
    """
    flatten = lambda l: [item for sublist in l for item in sublist]
    
    # create a graph in which nodes represent triangles;
    # nodes are connected if the corresponding triangles are adjacent and coplanar
    G = nx.Graph()
    G.add_nodes_from(range(len(triangles)))
    for ii, a in enumerate(triangles):
        for jj, b in enumerate(triangles):
            if (ii < jj): # test relationships only in one way as adjacency and co-planarity are bijective
                if is_adjacent(a, b):
                    if is_coplanar(a, b, np.pi / 180.):
                        G.add_edge(ii,jj)

    # triangles that belong to a connected component can be combined
    components = list(nx.connected_components(G))
    simplified = [set(flatten(triangles[index] for index in component)) for component in components]

    # need to reorder nodes so that patches are plotted correctly
    reordered = [reorder(face) for face in simplified]

    return reordered

def simplifyFaces(faces,verts):
    triangles = []
    for s in faces:
        sq = [
            (verts[s[0], 0], verts[s[0], 1], verts[s[0], 2]),
            (verts[s[1], 0], verts[s[1], 1], verts[s[1], 2]),
            (verts[s[2], 0], verts[s[2], 1], verts[s[2], 2])
        ]
        triangles.append(sq)
    new_faces = simplify(triangles)
    return new_faces

def is_adjacent(a, b):
    return len(set(a) & set(b)) == 2 # i.e. triangles share 2 points and hence a side

def is_coplanar(a, b, tolerance_in_radians=0):
    import sympy as syp
    a1, a2, a3 = a
    b1, b2, b3 = b
    plane_a = syp.Plane(syp.Point3D(a1), syp.Point3D(a2), syp.Point3D(a3))
    plane_b = syp.Plane(syp.Point3D(b1), syp.Point3D(b2), syp.Point3D(b3))
    if not tolerance_in_radians: # only accept exact results
        return plane_a.is_coplanar(plane_b)
    else:
        angle = plane_a.angle_between(plane_b).evalf()
        angle %= np.pi # make sure that angle is between 0 and np.pi
        return (angle - tolerance_in_radians <= 0.) or \
            ((np.pi - angle) - tolerance_in_radians <= 0.)

def reorder(vertices):
    """
    Reorder nodes such that the resulting path corresponds to the "hull" of the set of points.

    Note:
    -----
    Not tested on edge cases, and likely to break.
    Probably only works for convex shapes.

    """
    if len(vertices) <= 3: # just a triangle
        return vertices
    else:
        # take random vertex (here simply the first)
        reordered = [vertices.pop()]
        # get next closest vertex that is not yet reordered
        # repeat until only one vertex remains in original list
        vertices = list(vertices)
        while len(vertices) > 1:
            idx = np.argmin(get_distance(reordered[-1], vertices))
            v = vertices.pop(idx)
            reordered.append(v)
        # add remaining vertex to output
        reordered += vertices
        return reordered


def get_distance(v1, v2):
    v2 = np.array(list(v2))
    difference = v2 - v1
    ssd = np.sum(difference**2, axis=1)
    return np.sqrt(ssd)

def getPlotColors(numberOfGroups):
    pcs = []
    for i in range(0,numberOfGroups):
        pcs.append(pltc.rgb2hex(np.random.rand(3)))
    return pcs

def maxValuePlot(times,mPts,names,namespace,fs=16,lw=3,pcs=None,vName=''):
    '''  mPts rows correlated to times, columns correlated to different groups. '''
    numberOfGroups = mPts.shape[1]
    if pcs is None: pcs = getPlotColors(numberOfGroups)
    fig = plt.figure(figsize=(12,12))
    for i in range(0,numberOfGroups):
        plt.plot(times,mPts[:,i],color=pcs[i],label=names[i],linewidth=lw)
    plt.legend(fontsize=fs)
    plt.xlabel('time (s)',fontsize=fs)
    plt.ylabel('%s'%(vName),fontsize=fs)
    plt.tick_params(labelsize=fs)
    plt.tight_layout()
    plt.savefig('%s_maxTPlot.png'%(namespace),dpi=300)
    plt.show()
    return fig

def maxValueCSV(times, mPts, names, namespace):
    '''  mPts rows correlated to times, columns correlated to different groups. '''
    numberOfGroups = mPts.shape[1]
    header = 'Time,'
    for i in range(0,numberOfGroups):
        name = names[i].replace(',','_')
        header = header+name+','
    header = header[:-1]+'\n'
    data = np.append(np.reshape(times,(times.shape[0],1)),mPts,axis=1)
    np.savetxt('%s.csv'%(namespace),data,delimiter=',',header=header)
    return '%s.csv'%(namespace)    

def failureTimesCSV(times,mPts,names,thresholds,namespace):
    ''' mPts rows correlated to times, columns correlated to different groups. '''
    numberOfGroups = mPts.shape[1]
    ftimes = []
    header = ''
    for i in range(0,numberOfGroups):
        (name,threshold) = (names[i],thresholds[i])
        inds = np.where(mPts[:,i] >= threshold)[0]
        if inds.shape[0] > 0:
            tind = np.min(inds)
            ftime = times[tind]
        else:
            ftime = -1
        name = name.replace(',','_')
        ftimes.append(ftime)
        header = header+name+','
    ftimes = np.array(ftimes)
    ftimes = np.reshape(ftimes,(1,ftimes.shape[0]))
    header = header[:-1]+'\n'
    np.savetxt('%s_failuretimes.csv'%(namespace),ftimes,delimiter=',',header=header)
    return '%s_failuretimes.csv'%(namespace)

def polygonVisual(polygons,namespace,pcs=None,fs=16,fig=None,ax=None,
                  limits=[0,15,0,8,0,5]):
    '''
    polygons structure is the following:
        list of groups
          list of linked polygons
            scsp.ConvexHull(points)
    '''
    numberOfGroups = len(polygons)
    if pcs is None: pcs = getPlotColors(numberOfGroups)
    if fig is None: fig = plt.figure(figsize=(12,12))
    if ax is None: ax = a3.Axes3D(fig)
    for i in range(0,numberOfGroups):
        group = polygons[i]
        for p in group:
            pc = pcs[i]
            faces = simplifyFaces(p.simplices, p.points)
            for sq in faces:
                f = a3.art3d.Poly3DCollection([sq])
                f.set_color(pc)
                f.set_edgecolor('k')
                f.set_alpha(0.1)
                ax.add_collection3d(f)
        plt.xlim(limits[0],limits[1])
        plt.ylim(limits[2],limits[3])
        ax.set_zlim(limits[4],limits[5])
    plt.xlabel('x (m)',fontsize=fs)
    plt.ylabel('y (m)',fontsize=fs)
    ax.set_zlabel('z (m)',fontsize=fs)
    plt.tick_params(labelsize=fs)
    #ax.autoscale_view()
    #plt.tight_layout()
    plt.savefig('%s_polyvisual.png'%(namespace),dpi=300)
    
    return fig, ax

def smvVisual(obstructions,surfaces,namespace,fs=16,fig=None,ax=None,
              limits=[0,15,0,8,0,5]):
    if fig is None: fig = plt.figure(figsize=(12,12))
    if ax is None: ax = a3.Axes3D(fig)
    
    for obst in obstructions:
        pts, colors = getPtsFromObst(obst,surfaces)
        print(pts)
        print(colors)
        for pt, color in zip(pts,colors):
            f = a3.art3d.Poly3DCollection(pt)
            f.set_color(color)
            f.set_edgecolor('k')
            #f.set_alpha(1)
            ax.add_collection3d(f)
    plt.xlim(limits[0],limits[1])
    plt.ylim(limits[2],limits[3])
    ax.set_zlim(limits[4],limits[5])
    plt.xlabel('x (m)',fontsize=fs)
    plt.ylabel('y (m)',fontsize=fs)
    ax.set_zlabel('z (m)',fontsize=fs)
    plt.tick_params(labelsize=fs)
    plt.savefig('%s_smvvisual.png'%(namespace),dpi=300)
    
    return fig, ax

def buildSMVgeometry(file):
    with open(file,'r') as f:
        lines = f.readlines()
    inds = []
    for i in range(0,len(lines)):
        if lines[i][0] != ' ' and lines[i][0] != '-':
            inds.append(i)
    surfaces = []
    obstructions = []
    for ind in inds:
        if 'SURFACE' in lines[ind]:
            sname = ' '.join(lines[ind+1].split())
            (Tign,eps) = (lines[ind+2].split()[0],lines[ind+2].split()[1])
            (stype,t_width,t_height) = (lines[ind+3].split()[0],lines[ind+3].split()[1],lines[ind+3].split()[2])
            (c1,c2,c3,c4) = (lines[ind+3].split()[3],lines[ind+3].split()[4],lines[ind+3].split()[5],lines[ind+3].split()[6])
            surfaces.append([sname,Tign,eps,stype,t_width,t_height,c1,c2,c3,c4])
        if 'OBST' in lines[ind] and '_OBST' not in lines[ind]:
            nObst = int(lines[ind+1].split()[0])
            for i in range(0,nObst):
                obst = [float(x) for x in lines[ind+i+2].split()]
                obstructions.append(obst)
    return surfaces, obstructions

def getPtsFromObst(obst,surfaces):
    pts = []
    colors = []
    pts = np.array([[obst[0],obst[2],obst[4]],
                   [obst[0],obst[2],obst[5]],
                   [obst[1],obst[2],obst[5]],
                   [obst[1],obst[2],obst[4]],
                   [obst[0],obst[3],obst[5]],
                   [obst[1],obst[3],obst[5]],
                   [obst[1],obst[3],obst[4]],
                   [obst[0],obst[3],obst[4]]])
    
    # y-negative surface
    #pts.append([(obst[0],obst[2],obst[4]),(obst[1],obst[2],obst[4]),
    #            (obst[1],obst[2],obst[5]),(obst[0],obst[2],obst[5])])
    surf = surfaces[int(obst[7])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    # y-positive surface
    #pts.append([(obst[0],obst[3],obst[4]),(obst[1],obst[3],obst[4]),
    #            (obst[1],obst[3],obst[5]),(obst[0],obst[3],obst[5])])
    surf = surfaces[int(obst[8])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    # x-negative surface
    #pts.append([(obst[0],obst[2],obst[4]),(obst[0],obst[2],obst[5]),
    #            (obst[0],obst[3],obst[5]),(obst[0],obst[3],obst[4])])
    surf = surfaces[int(obst[9])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    # x-positive surface
    #pts.append([(obst[1],obst[2],obst[4]),(obst[1],obst[2],obst[5]),
    #            (obst[1],obst[3],obst[5]),(obst[1],obst[3],obst[4])])
    surf = surfaces[int(obst[10])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    # z-negative surface
    #pts.append([(obst[0],obst[2],obst[4]),(obst[1],obst[2],obst[4]),
    #            (obst[1],obst[3],obst[4]),(obst[0],obst[3],obst[4])])
    surf = surfaces[int(obst[11])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    # z-positive surface
    #pts.append([(obst[0],obst[2],obst[5]),(obst[1],obst[2],obst[5]),
    #            (obst[1],obst[3],obst[5]),(obst[0],obst[3],obst[5])])
    surf = surfaces[int(obst[12])]
    colors.append((float(surf[6]),float(surf[7]),float(surf[8]),float(surf[9])))
    
    return pts, colors
    
