import pyvista as pv
import numpy as np

plotter = pv.Plotter()
#File Names must be input here
Filenameroot=''
Filenamexyz=f'{Filenameroot}_static_o_DS2_DEN.xyz'
Filenamecoeff=f'{Filenameroot}_static-coeff'
Filenamecell=f'{Filenameroot}_static-cell'
Filenamegeo=f'{Filenameroot}_static-geo'

#=======================================================================
#=================Preferred Methods Of Scaling==========================
#=======================================================================
#Both of ScaleLobe and Scaleatom are controlled by sliders in the rendered chemical pressure
ScaleLobe=100 #Scales the size of the chemical pressure lobes. Default is 100
Scaleatom=1  #Can adjust for the size of atoms. Default is 1, 

unitcelllinewidth=5 #Default is 5

#Notes
# This code works for 4 element compounds
# This code does not work for compounds with more than 4 elements 
# For help with this code contact the Fredrickson Group at UW Madison

#Reading xyz from the xyz file.
def ReadTxt(Filenamexyz):
    with open(f"{Filenamexyz}", 'r') as f:
        lines = f.readlines()
        nameofatoms = []
        for i, line in enumerate(lines):
            line = line.strip().split()
            line= [var.split("'") for var in line]
            line=[element for var in line for element in var if element != '']
            if i == 0:
                num_atoms =int(line[0])
                xyz= np.zeros((num_atoms,3))
            elif i ==1:
                continue
            else:
                coords = line[1:]
                xyz[i-2, :] = np.array([float(coord) for coord in coords])
            nameofatoms.append(line[0])
        del(nameofatoms[0])
        return xyz,nameofatoms
xyz,nameofatoms=ReadTxt(Filenamexyz)
atom1=nameofatoms[0] 
atom2=0

#This finds and labels each atom in the system 
atoms={}
count=0
for thing in nameofatoms:
    if thing not in atoms.keys():
        atoms[thing] = count
        count+=1

#This converts the strings in the nameofatoms list to numbers in order to be read below.
#Converts atom1 to 0, atom2 to 1, atom3 to 2, and atom4 to 3
numberofatoms=[]
numberofatoms = [atoms[s] for s in nameofatoms]

#Getting coeff out of the coeff file
def ReadTxt2(Filenamecoeff):
    with open(f"{Filenamecoeff}", 'r') as f:
        lines = f.readlines()
        coeffvalues = []
        for i, line in enumerate(lines):
            line = line.strip().split()
            line= [var.split("'") for var in line]
            line=[element for var in line for element in var if element != '']
            del(line[0])
            line = [float(s) for s in line]
            coeffvalues.append(line[0])
        return coeffvalues
coeffvalues=ReadTxt2(Filenamecoeff)

#Changing Coeff to be a list of lists for the next step
CPcoeff_temporary = []
for i in range(len(coeffvalues) // 49):
    CPcoeff_temporary.append([])
    for j in range(49):
        CPcoeff_temporary[i].append(coeffvalues[49 * i + j])

coeffvalues = CPcoeff_temporary

#Getting cell out of the cell file
def ReadTxt3(Filenamecell):
    cellvalues = []
    with open(Filenamecell, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                raise ValueError("Each line must contain exactly 3 numbers")
            cellvalues.append([float(x) for x in parts])
    return np.array(cellvalues)
cellvalues= ReadTxt3(Filenamecell)

#getting the geo out of the file
def ReadTxt4(filename):
    geoelements = []
    geovalues = []
    with open(filename, 'r') as f:
        for i, line in enumerate(f, start=1):
            parts = line.strip().split()

            if len(parts) != 4:
                raise ValueError(f"Line {i} does not contain exactly 4 values")
            
            numbers = [float(parts[1]), float(parts[2]), float(parts[3])]

            geoelements.append(parts[0])
            geovalues.append(numbers)

        return geoelements, np.array(geovalues, dtype=float)
geoelements, geovalues=ReadTxt4(Filenamegeo)

#This generates the chemical pressure values
def Chemicalpressure(coeffvalues):
    Geo_Cps =[]
    for val in coeffvalues:
        scale=1    #Use the scales at the top of the code this scales unpredictably

        positive_surf = pv.ParametricEllipsoid(xradius = 1, yradius = 1, zradius = 1,)
        negative_surf = pv.ParametricEllipsoid(xradius = 1, yradius = 1, zradius = 1,)
        for i in range(len(positive_surf.points)):
            r=0
            x = positive_surf.points[i][0]
            y = positive_surf.points[i][1]
            z = positive_surf.points[i][2]
            theta = np.arccos(z)

            phi=0
            if((theta > 0) and (theta < 3.14159)): 
                if(x/np.sin(theta) > 1): 
                    phi = 0
                elif(x/np.sin(theta) < -1):
                    phi = np.pi
                else:
                    phi = np.arccos(x/np.sin(theta))

            if y < 0:
                phi = -phi
            #This is a Legendre Polynomial expansion. It goes to L=4 this can be expanded
            #for additional quality. 

            #L = 0
            r += scale * 0.5 * (1 / np.pi)**0.5 * val[0]

            #L = 1
            r += scale * 0.5 * (3 / np.pi)**0.5 * np.cos(theta) * val[1]
            r += -scale * (3 / (8 * np.pi))**0.5 * np.sin(theta) * np.cos(phi) * np.sqrt(2) * val[2]
            r += -scale * (3 / (8 * np.pi))**0.5 * np.sin(theta) * np.sin(phi) * np.sqrt(2) * val[3]

            #L = 2
            r += scale * 0.25 * (5 / np.pi)**0.5 * (3*np.cos(theta)**2 - 1) * val[4]
            r += -scale * 0.5 * (15 / np.pi)**0.5 * np.sin(theta) * np.cos(phi) * np.cos(theta) * val[5]
            r += -scale * 0.5 * (15 / np.pi)**0.5 * np.sin(theta) * np.sin(phi) * np.cos(theta) * val[6]
            r += scale * 0.25 * (15 / np.pi)**0.5 * np.sin(theta)**2 * np.cos(2*phi) * val[7]
            r += scale * 0.25 * (15 / np.pi)**0.5 * np.sin(theta)**2 * np.sin(2*phi) * val[8]

            #L = 3
            r += scale * 0.25 * (7 / np.pi)**0.5 * (5*np.cos(theta)**3 - 3*np.cos(theta)) * val[9]
            r += -scale * 0.125 * (21 / np.pi)**0.5 * np.sin(theta) * (5*np.cos(theta)**2 - 1) * np.cos(phi) * np.sqrt(2) * val[10]
            r += -scale * 0.125 * (21 / np.pi)**0.5 * np.sin(theta) * (5*np.cos(theta)**2 - 1) * np.sin(phi) * np.sqrt(2) * val[11]
            r += scale * 0.25 * (105 / (2*np.pi))**0.5 * np.sin(theta)**2 * np.cos(theta) * np.cos(2*phi) * np.sqrt(2) * val[12]
            r += scale * 0.25 * (105 / (2*np.pi))**0.5 * np.sin(theta)**2 * np.cos(theta) * np.sin(2*phi) * np.sqrt(2) * val[13]
            r += -scale * 0.125 * (35 / np.pi)**0.5 * np.sin(theta)**3 * np.cos(3*phi) * np.sqrt(2) * val[14]
            r += -scale * 0.125 * (35 / np.pi)**0.5 * np.sin(theta)**3 * np.sin(3*phi) * np.sqrt(2) * val[15]

            #L = 4
            r += scale * (3/16) * (1 / np.pi)**0.5 * (35*np.cos(theta)**4 - 30*np.cos(theta)**2 + 3) * val[16]
            r += -scale * (3/8) * (5 / np.pi)**0.5 * np.sin(theta) * (7*np.cos(theta)**3 - 3*np.cos(theta)) * np.cos(phi) * np.sqrt(2) * val[17]
            r += -scale * (3/8) * (5 / np.pi)**0.5 * np.sin(theta) * (7*np.cos(theta)**3 - 3*np.cos(theta)) * np.sin(phi) * np.sqrt(2) * val[18]
            r += scale * (3/8) * (5 / (2*np.pi))**0.5 * np.sin(theta)**2 * (7*np.cos(theta)**2 - 1) * np.cos(2*phi) * np.sqrt(2) * val[19]
            r += scale * (3/8) * (5 / (2*np.pi))**0.5 * np.sin(theta)**2 * (7*np.cos(theta)**2 - 1) * np.sin(2*phi) * np.sqrt(2) * val[20]
            r += -scale * (3/8) * (35 / np.pi)**0.5 * np.sin(theta)**3 * np.cos(theta) * np.cos(3*phi) * np.sqrt(2) * val[21]
            r += -scale * (3/8) * (35 / np.pi)**0.5 * np.sin(theta)**3 * np.cos(theta) * np.sin(3*phi) * np.sqrt(2) * val[22]
            r += scale * (3/16) * (35 / (2*np.pi))**0.5 * np.sin(theta)**4 * np.cos(4*phi) * np.sqrt(2) * val[23]
            r += scale * (3/16) * (35 / (2*np.pi))**0.5 * np.sin(theta)**4 * np.sin(4*phi) * np.sqrt(2) * val[24]

            r *= 4 * np.pi

            if r < 0:
                positive_surf.points[i] = np.zeros(3).copy()
                negative_surf.points[i]= np.array([ r * np.sin(theta) * np.cos(phi), 
                                                    r * np.sin(theta) * np.sin(phi),
                                                    r * np.cos(theta)]).copy()
            else:
                negative_surf.points[i] = np.zeros(3).copy()
                positive_surf.points[i]= np.array([ r * np.sin(theta) * np.cos(phi), 
                                                    r * np.sin(theta) * np.sin(phi),
                                                    r * np.cos(theta)]).copy()
        Geo_Cps.append([positive_surf,negative_surf])

    return Geo_Cps

Geo_Cps = Chemicalpressure(coeffvalues)
atom_data = []
#Graphing the spheres just the atom location
for i, coords in enumerate(xyz):
    Sphere = pv.Sphere()
    orig_pts = Sphere.points.copy()
    Sphere.points *=Scaleatom
    Sphere=Sphere.translate(coords, inplace=True)
    if numberofatoms[i] == 0:
        plotter.add_mesh(Sphere, color='blue')
    elif numberofatoms[i] == 1:
        plotter.add_mesh(Sphere, color='red')
    elif numberofatoms[i] == 2:
        plotter.add_mesh(Sphere, color='orange')
    else:  #If more than 4 atoms are needed increase the size of this section to the number of atoms needed
        plotter.add_mesh(Sphere, color='green')
    atom_data.append((Sphere, orig_pts, np.array(coords)))
    
#Needed for scaling
lobe_data = []

#Appling the chemical presure lobe
def apply_CP(geovalues, cellvalues, Geo_Cps, xyz):
    natoms_geo = len(geovalues)
    natoms_coords = len(xyz)
    for j1 in range(-4, 5):
        for j2 in range(-4, 5):
            for j3 in range(-4, 5):
                for k1 in range(natoms_geo):
                    for k2 in range(natoms_coords):
                        x1 = (geovalues[k1][0] + cellvalues[0][0]*j1 +cellvalues[1][0]*j2 +cellvalues[2][0]*j3)
                        y1 = (geovalues[k1][1] + cellvalues[0][1]*j1 + cellvalues[1][1]*j2 + cellvalues[2][1]*j3)
                        z1 = (geovalues[k1][2] + cellvalues[0][2]*j1 + cellvalues[1][2]*j2 + cellvalues[2][2]*j3)
                        x2 = xyz[k2][0]
                        y2 = xyz[k2][1]
                        z2 = xyz[k2][2]
                        dist = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5 #This filters out all objects that are on top of other objects should filter to the number of atoms past this point
                        if(dist < 0.1): 
                            for i in range(0,1):
                                ScaleXYZ=1  #This can be ignored do not change from 1, 
                                #this changes how much the chemical pressures scale if changed from 1 
                                #only use for bug fixes
                                xyztemp=[x2/ScaleXYZ,y2/ScaleXYZ,z2/ScaleXYZ]
                                ellipsoid1 = Geo_Cps[k1][1].copy() # negative lobe
                                ellipsoid2 = Geo_Cps[k1][0].copy() # positive lobe

                                orig_pts1 = ellipsoid1.points.copy()
                                orig_pts2 = ellipsoid2.points.copy()

                                ellipsoid1.points *= ScaleLobe  
                                ellipsoid2.points *= ScaleLobe
                            
                                ellipsoid1.translate(xyztemp, inplace=True)
                                ellipsoid2.translate(xyztemp, inplace=True)

                                lobe_data.append((ellipsoid1, orig_pts1,xyztemp))
                                lobe_data.append((ellipsoid2, orig_pts2,xyztemp))
                            
                                plotter.add_mesh(ellipsoid1, style='surface', color='black')  #Black in website code,  Negative
                                plotter.add_mesh(ellipsoid2, style='surface', color='white')  #White in website code,  Positive

CreateCPs = apply_CP(geovalues, cellvalues, Geo_Cps, xyz)
#Sliders that scaleatoms and scale the chemical pressure lobes
def Scale_lobes(value):
    for mesh, orig_pts, center in lobe_data:
        mesh.points = orig_pts *value +center
    plotter.render()

plotter.add_slider_widget(
    callback=Scale_lobes,
    rng = [1,300],
    value = ScaleLobe,
    title= 'Scaling Lobes',
    pointa=(0,0.1),
    pointb=(0.4,0.1),
)

def scale_atoms(value):
    for mesh, orig_pts, center in atom_data:
        mesh.points = orig_pts *value +center
    plotter.render()

plotter.add_slider_widget(
    callback=scale_atoms,
    rng = [0.01,3],
    value = Scaleatom,
    title= 'Scaling Atoms',
    pointa=(0,0.25),
    pointb=(0.4,0.25),
)

#Creates points for the unit cell
cellvalue0=cellvalues[0]
cellvalue1=cellvalues[1]
cellvalue2=cellvalues[2]

Zero= np.array([0,0,0])
cellvalue01 = cellvalue0 + cellvalue1
cellvalue02 = cellvalue0 + cellvalue2
cellvalue12 = cellvalue1 + cellvalue2
cellvalue012 = cellvalue0 + cellvalue1 + cellvalue2

unitcellpoints = np.array([Zero,cellvalue0,Zero,cellvalue1,cellvalue0,cellvalue01,cellvalue1,cellvalue01, #Bottom Face
                           Zero,cellvalue2,cellvalue0,cellvalue02,cellvalue1,cellvalue12,cellvalue01,cellvalue012, #Vertical Edges
                           cellvalue2,cellvalue02,cellvalue2,cellvalue12,cellvalue02,cellvalue012,cellvalue12,cellvalue012], #Top Face
                           dtype=float)

#Create lines for unit cell
plotter.add_lines(unitcellpoints, color = 'green', width = unitcelllinewidth,)   #The value for the width can be changed at the top of the code
plotter.show()
