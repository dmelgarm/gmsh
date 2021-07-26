'''
Some tools to deal with preparing and extracting info from gmsh inpout and output files
DMM 01/2015
'''




def xyz2gmsh(fout,lon,lat,depth,coord_type='UTM',projection_zone=None,point_offset=0):
    '''
    This is for taking Slab 1.0 grids and converting them to something gmsh can 
    deal with. It takesLat,lon,depth numpy arrays converts them to UTM and then 
    make a gmsh input file
    '''
    f=open(fout,'w')
    if coord_type=='ECEF':
        x,y,z=llz2xyz(lon,lat,depth)
    if coord_type=='UTM':
        x,y=llz2utm(lon,lat,projection_zone)
        z=depth.copy()
        x=x/1000
        y=y/1000
    for k in range(len(x)):
        #Point(1000) = { 0.99975300, 0.00129500, 0.00000000, 0.01};
        line='Point('+str(k+1+point_offset)+') = {%.6f, %.6f, %.6f, 0.1};\n' %(x[k],y[k],z[k])
        f.write(line)
    f.close()


def gmsh2ascii(fout,msh_in,utm_zone='10T',flip_lon=False,flip_strike=False):
    '''
    This prepares the output of gmsh for something usable by mudpy.
    Read triangle information from gmsh msh file and extract node coordinates, 
    coordinates of the centroid, area, strike and dip and convert back to lat/lon
    '''
    from numpy import c_,savetxt
    nodes,elements=read_msh(msh_in)
    #Convert elements to node-coordinate format
    ncoords=nodes2coords(nodes,elements)
    #Get centroids
    centroids=get_centroids(ncoords)
    #Get area
    areas=get_areas(ncoords)
    #get vertex lengths
    L=get_vertex_lengths(ncoords)
    #Get strike and dip
    strike,dip=get_geometry(ncoords)
    if flip_strike:
        strike += 180
    #Convert from UTM to lat,lon
    ncoords,centroids=coords2latlon(ncoords,centroids,zone=utm_zone)
    if flip_lon==True:
        centroids[:,1]=centroids[:,1]-360
        ncoords[:,1]-=360
        ncoords[:,4]-=360
        ncoords[:,7]-=360
    #And done, save
    out=c_[centroids,ncoords[:,1:],L.mean(1),areas,strike,dip]
    savetxt(fout,out,fmt='%i\t%12.6f\t%12.6f\t%12.8f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%12.6f\t%9.2f\t%9.2f\t%9.2f\t%9.2f',header='  fault No. , centroid(lon,lat,z[km]) , node2(lon,lat,z[km]) , node3(lon,lat,z[km]) , mean vertex length(km) , area(km^2) , strike(deg) , dip(deg)')
    
def make_mudpy_fault(faultin,faultout,vr=2.6,rise_time_all=5.0,flip_strike=False):
    '''
    Convert to mudpy format fault
    '''
    from numpy import genfromtxt,ones,c_,savetxt
    f=genfromtxt(faultin)
    #Fault numbers
    n=f[:,0]
    #Coordinates, use the centroid
    lon=f[:,1]
    lat=f[:,2]
    depth=-f[:,3] #In positive km
    if flip_strike==True:
        strike=f[:,15]-180
    else:
        strike=f[:,15]
    dip=f[:,16]
    triangle=0.5*ones(len(f))
    if rise_time_all==None:
        rise_time=f[:,13]/vr  #Mean vertex length divided by rupture speed
    else:
        rise_time=rise_time_all*ones(len(lon))
    length=(f[:,14]**0.5)*1000  #Square root of the area (in m)
    #Put together and write file
    out=c_[n,lon,lat,depth,strike,dip,triangle,rise_time,length,length]
    h='No,lon,lat,depth(km),strike,dip,type,rise time(s),length(km),width(km)'
    savetxt(faultout,out,fmt='%i\t%10.6f\t%10.6f\t%12.8f\t%.2f\t%.2f\t%6.3f\t%6.3f\t%10.2f\t%10.2f',header=h)
    
def get_perimeter(gmsh):
    '''
    Get the eprimeter of the fault model
    '''
    from numpy import array,r_,float64,zeros,c_,savetxt
    from pyproj import Proj
    
    #Read all points and their coordinates
    f=open(gmsh,'r')
    x=array([])
    y=array([])
    z=array([])
    while True:
        l=f.readline()
        if 'Point' in l: #Line is a point
            x=r_[x,float64(l.split('=')[-1].split(',')[0].replace('{',''))]
            y=r_[y,float64(l.split('=')[-1].split(',')[1])]
            z=r_[z,float64(l.split('=')[-1].split(',')[2].replace('}',''))]
        else:
            break
    f.close()
    #Iquique
    #i1=array([1,180,338,499,661,825,992,1163,1336,1512,1691,1872,2055,2242])-1
    #i2=array([2242,2237,2233,2225,2216,2207,2197,2191,2179,2170,2159,2144,2131,2115,2107,2096,2085,2077,2069,2059])-1
    #i3=array([24,40,190,346,504,665,829,1166,1516,1875,2059])-1
    #i4=array([24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1])-1
    #i=r_[i3,i2[::-1],i1[::-1],i4[::-1]] #This si the perimeter
    #Maule
    #i1=array([ 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 19, 21, 24, 26, 29, 30, 33, 34, 35, 38, 41, 44, 46, 47, 49, 51, 53, 54, 55])-1
    #i2=array([55, 56, 385, 714, 1042, 1365, 1686, 2006, 2325, 2642, 2961, 3276, 3594, 3912])-1
    #i3=array([3912, 3918, 3931, 3943, 3956, 3976, 3992, 4004, 4023, 4044, 4062, 4079, 4092, 4107, 4127, 4142, 4155, 4169, 4182, 4194, 4205, 4217, 4229])-1
    #i4=array([4229, 3911, 3593, 3275, 2960, 2641, 2324, 2005, 1685, 1364, 1041, 713, 384, 1])-1
    #Coquimbo
    #Maule
    #i1=array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29])-1
    #i2=array([1, 179, 331, 481, 633, 788, 942, 1098, 1258, 1418, 1580, 1743, 1906, 2072])-1
    #i3=array([2072, 2061, 2049, 2037, 2024, 2010, 1991, 1972, 1958, 1947, 1936, 1924, 1915, 1908])-1
    #i4=array([29, 30, 181, 332, 482, 634, 789, 943, 1100, 1259, 1419, 1581, 1744, 1908])-1
    #i=r_[i1,i4,i3[::-1],i2[::-1]] #This si the perimeter
    # Cascadia
    i1=array([16039, 16702, 17371, 18042, 688, 1375, 2064, 3462, 4164, 4875, 5587, 6304, 7033, 8514, 10013, 10769, 11531, 12295, 13062, 13834, 14607, 15383])
    i2=array([16039, 16035, 16032, 16028, 16025, 16021, 16019, 16015, 16012, 16008, 16004, 16002, 16001, 15999, 15994, 15991, 15989, 15986, 15984, 15981, 15978, 15974, 15970, 15966, 15963, 15961, 15957, 15951, 15947, 15942, 15938, 15933, 15928, 15923, 15921, 15919, 15916, 15913, 15909, 15906, 15902, 15899, 15894, 15886, 15882, 15876, 15869, 15863, 15856, 15849, 15844, 15839, 15834, 15831, 15828, 15824, 15819, 15814, 15808, 15792, 15777, 15772, 15764, 15757, 15751, 15743, 15732, 15727, 15714, 15700, 15689, 15676, 15662, 15650, 15644, 15637, 15629, 15625, 15616, 15609, 15598, 15589, 15579, 15572, 15565, 15557, 15549, 15536, 15526, 15518, 15507, 15492, 15475, 15462, 15442, 15424, 15407, 15387])
    i3=array([15387, 16703, 17372, 4, 689, 1376, 2065, 2761, 3463, 4166, 4876, 5588, 6305, 7034, 7773, 8516, 9263, 10014, 10770, 11532, 12296, 13063, 13835, 14609])
    i4=array([14609, 14617, 14626, 14639, 14650, 14661, 14673, 14684, 14702, 14714, 14727, 14739, 14751, 14761, 14768, 14777, 14787, 14800, 14809, 14819, 14829, 14835, 14844, 14853, 14860, 14867, 14875, 14880, 14883, 14886, 14889, 14894, 14897, 14906, 14915, 14920, 14927, 14934, 14941, 14948, 14952, 14956, 14962, 14970, 14984, 14998, 15002, 15007, 15027, 15044, 15051, 15067, 15076, 15084, 15092, 15098, 15111, 15120, 15133, 15146, 15160, 15178, 15189, 15192, 15197, 15202, 15209, 15217, 15226, 15238, 15252, 15262, 15272, 15280, 15289, 15298, 15306, 15313, 15318, 15323, 15331, 15340, 15350, 15360, 15365, 15374, 15383])
    i1=i1-1
    i2=i2-1
    i3=i3-1
    i4=i4-1
    i=r_[i2,i3,i4,i1[::-1]]

    
    x=x[i]
    y=y[i]
    #Convert back to lat,lon
    zone='10T'
    p = Proj(proj='utm',zone=zone,ellps='WGS84')
    LL=p(x*1000,y*1000,inverse=True)
    lon=LL[0]
    lat=LL[1]
    #Write
    out=c_[lon,lat]
    savetxt('/Users/dmelgar/Cascadia/mesh/perimeter.xy',out,fmt='%10.6f\t%10.6f')

        
        
        
        
        

def llz2xyz(lon,lat,depth):
    '''
    Convert lat,lon to Earth centered Earth fixed coords
    '''
    from numpy import deg2rad,pi,sin,cos
    R=6371 #Radius of the Earth
    r=R+depth #Convert to spherical radius
    phi=(pi/2)-deg2rad(lat)
    theta=deg2rad(lon)
    #Apply conversions
    x=r*sin(theta)*cos(phi)
    y=r*sin(theta)*sin(phi)
    z=r*cos(theta)
    return x,y,z
    

def llz2utm(lon,lat,projection_zone='None'):
    '''
    Convert lat,lon to UTM
    '''
    from numpy import zeros,where,chararray
    import utm
    from pyproj import Proj
    from scipy.stats import mode
    
    x=zeros(lon.shape)
    y=zeros(lon.shape)
    zone=zeros(lon.shape)
    b=chararray(lon.shape)
    if projection_zone==None:
        #Determine most suitable UTM zone
        for k in range(len(lon)):
            #x,y,zone[k],b[k]=utm.from_latlon(lat[k],lon[k]-360)
            x,y,zone[k],b[k]=utm.from_latlon(lat[k],lon[k])
        zone_mode=mode(zone)
        i=where(zone==zone_mode)[0]
        letter=b[i[0]]
        z=str(int(zone[0]))+letter
    else:
        z=projection_zone
    p = Proj(proj='utm',zone=z,ellps='WGS84')
    x,y=p(lon,lat)
    return x,y
    
def read_msh(msh_in):
    '''
    Parse msh format and extract nodes and elements and their coordinates
    '''
    from numpy import zeros,float64

    #Parse mesh file
    f=open(msh_in)
    while True:
        line=f.readline()
        print(line)
        if '$Nodes' in line:
            print(line)
            num_nodes=int(f.readline())
            nodes=zeros((num_nodes,4))
            for k in range(num_nodes):
                L=f.readline()
                nodes[k,:]=float64(L.split())
        if '$Elements' in line:
            num_elements=int(f.readline())
            elements=zeros((num_elements,8))
            for k in range(num_elements):
                L=f.readline()
                elements[k,:]=float64(L.split())
            break
    return nodes,elements

    
def nodes2coords(nodes,elements):
    '''
    Associate node coordinates to each element
    '''
    from numpy import zeros,where

    ncoords=zeros((len(elements),10))
    for k in range(len(elements)):
        node1=int(elements[k,5])
        node2=int(elements[k,6])
        node3=int(elements[k,7])
        ncoords[k,0]=k+1
        #Get node 1 coordinates
        i=where(nodes[:,0]==node1)[0]
        ncoords[k,1:4]=nodes[i,1:4]
        #Get node 2 coordinates
        i=where(nodes[:,0]==node2)[0]
        ncoords[k,4:7]=nodes[i,1:4]
        #Get node 1 coordinates
        i=where(nodes[:,0]==node3)[0]
        ncoords[k,7:10]=nodes[i,1:4]
    return ncoords
    
    
def get_centroids(ncoords):
    '''
    From node coordinates get element centroid
    '''
    from numpy import arange,zeros
    
    centroids=zeros((len(ncoords),4))
    centroids[:,0]=arange(1,len(ncoords)+1)
    centroids[:,1]=(1./3)*(ncoords[:,1]+ncoords[:,4]+ncoords[:,7]) #x centroids
    centroids[:,2]=(1./3)*(ncoords[:,2]+ncoords[:,5]+ncoords[:,8]) #y centroids
    centroids[:,3]=(1./3)*(ncoords[:,3]+ncoords[:,6]+ncoords[:,9]) #z centroids
    return centroids


def get_areas(ncoords):
    '''
    Compute area of triangles area=||ABxAC||/2
    '''
    from numpy import zeros,cross,array
    from numpy.linalg import norm
    
    areas=zeros(len(ncoords))
    for k in range(len(ncoords)):
        #Form vertex vectors
        ABx=ncoords[k,4]-ncoords[k,1]
        ABy=ncoords[k,5]-ncoords[k,2]
        ABz=ncoords[k,6]-ncoords[k,3]
        AB=array([ABx,ABy,ABz])
        ACx=ncoords[k,7]-ncoords[k,1]
        ACy=ncoords[k,8]-ncoords[k,2]
        ACz=ncoords[k,9]-ncoords[k,3]
        AC=array([ACx,ACy,ACz])
        areas[k]=norm(cross(AB,AC))/2
    return areas
    
def get_vertex_lengths(ncoords):
    '''
    For every element compute the legnth of its 3 vertices
    '''
    from numpy import zeros,sqrt
    L=zeros((len(ncoords),3))
    L[:,0]=((ncoords[:,1]-ncoords[:,4])**2+(ncoords[:,2]-ncoords[:,5])**2+(ncoords[:,3]-ncoords[:,6])**2)**0.5
    L[:,1]=((ncoords[:,1]-ncoords[:,7])**2+(ncoords[:,2]-ncoords[:,8])**2+(ncoords[:,3]-ncoords[:,9])**2)**0.5
    L[:,2]=((ncoords[:,4]-ncoords[:,7])**2+(ncoords[:,5]-ncoords[:,8])**2+(ncoords[:,6]-ncoords[:,9])**2)**0.5
    return L
        
    
def get_geometry(ncoords):
    '''
    Get strike and dip of each triangle, this assumes coordinates are UTM or similar
    cartesian system, i.e. y is due north. This won't work for an Earth centered
    reference frame.
    '''
    
    from numpy import zeros,cross,array,arctan,arcsin,rad2deg,where,deg2rad,sin,cos
    from numpy.linalg import norm
    
    strike=zeros(len(ncoords))
    dip=zeros(len(ncoords))
    for k in range(len(ncoords)):
        #Form vertex vectors
        ABx=ncoords[k,4]-ncoords[k,1]
        ABy=ncoords[k,5]-ncoords[k,2]
        ABz=ncoords[k,6]-ncoords[k,3]
        AB=array([ABx,ABy,ABz])
        ACx=ncoords[k,7]-ncoords[k,1]
        ACy=ncoords[k,8]-ncoords[k,2]
        ACz=ncoords[k,9]-ncoords[k,3]
        AC=array([ACx,ACy,ACz])
        #Compute normal vector
        n=cross(AB,AC)
        #Extract equation of the plane parameters
        a=n[0]
        b=n[1]
        c=n[2]
        d=a*ncoords[k,1]+b*ncoords[k,2]+c*ncoords[k,3]
        #Compute strike
        strike[k]=rad2deg(arctan(-b/a))
        #Compute dip
        beta=deg2rad(strike[k]+90)
        m=array([sin(beta),cos(beta),0]) #Points in dip direction
        n=array([a,b,c]) #Normal to the plane
        dip[k]=abs(rad2deg(arcsin(m.dot(n)/(norm(m)*norm(n)))))
    #Arange strike from 0 to 360
    i=where(strike<0)[0]
    strike[i]=360+strike[i]
    return strike,dip
    
def coords2latlon(ncoords,centroids,zone):
    '''
    Convert from utm to lat.lon
    '''
    from pyproj import Proj
    

    p = Proj(proj='utm',zone=zone,ellps='WGS84')
    

    LL=p(ncoords[:,1]*1000,ncoords[:,2]*1000,inverse=True)
    ncoords[:,2]=LL[1] ; ncoords[:,1]=LL[0]
    LL=p(ncoords[:,4]*1000,ncoords[:,5]*1000,inverse=True)
    ncoords[:,5]=LL[1] ; ncoords[:,4]=LL[0]
    LL=p(ncoords[:,7]*1000,ncoords[:,8]*1000,inverse=True)
    ncoords[:,8]=LL[1] ; ncoords[:,7]=LL[0]
    ncoords[:,1]=360+ncoords[:,1]
    ncoords[:,4]=360+ncoords[:,4]
    ncoords[:,7]=360+ncoords[:,7]
    LL=p(centroids[:,1]*1000,centroids[:,2]*1000,inverse=True)
    centroids[:,2]=LL[1] ; centroids[:,1]=LL[0]
    centroids[:,1]=360+centroids[:,1]
    return ncoords,centroids
    
def get_utm_zone(lon,lat):
    import utm
    from numpy import zeros,mode,where
    #Determine most suitable UTM zone
    zone=zeros(len(lon))
    for k in range(len(lon)):
        x,y,zone[k],b[k]=utm.from_latlon(lat[k],lon[k]-360)
    zone_mode=mode(zone)
    print(zone)
    i=where(zone==zone_mode)[0]
    letter=b[i[0]]
    z=str(int(zone[0]))+letter
    return z
    
def plot_strike(centroids,strike):
    '''
    Quick scatter plot ot visualize strike values
    '''
    import matplotlib.pyplot as plt
    plt.figure()
    plt.scatter(centroids[:,1],centroids[:,2],marker='x')
    for k in range(len(strike)):
        plt.annotate(str(int(strike[k])), xy = (centroids[k,1],centroids[k,2]))
    plt.axis('equal')
    plt.show()
    
    
def plot_dip(centroids,dip):
    '''
    Quick scatter plot ot visualize dip values
    '''
    import matplotlib.pyplot as plt
    plt.figure()
    plt.scatter(centroids[:,1],centroids[:,2],marker='x')
    for k in range(len(dip)):
        plt.annotate(str(int(dip[k])), xy = (centroids[k,1],centroids[k,2]))
    plt.axis('equal')
    plt.show()
    

        
        
        
    
    
