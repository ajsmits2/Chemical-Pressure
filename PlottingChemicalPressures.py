import pyvista as pv
import numpy as np


plotter = pv.Plotter()
Filenamexyz='Ti_TiSe2_static_o_DS2_DEN.xyz'
Filenamecoeff='Ti_TiSe2_static-coeff'
filenamecell='Ti_TiSe2_static-cell'
filenamegeo='Ti_TiSe2_static-geo'

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

#Only works for 2atom coumpounds, need to fix This finds the second atom
for i, thing in enumerate(nameofatoms):
    if thing in nameofatoms != atom1:
        atom2=thing

#This converts the strings in the nameofatoms list to numbers in order to be read below.
Conversion = {str(atom1):0, str(atom2):1}
numberofatoms=[]
numberofatoms = [Conversion[s] for s in nameofatoms]

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
print(type(coeffvalues[0]))

#Getting cell out of the cell file
def ReadTxt3(filenamecell):
    cellvalues = []
    with open(filenamecell, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                raise ValueError("Each line must contain exactly 3 numbers")
            cellvalues.append([float(x) for x in parts])
    return np.array(cellvalues)
cellvalues= ReadTxt3(filenamecell)

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
geoelements, geovalues=ReadTxt4(filenamegeo)

#Variables: coeffvalues,cellvalues,geovalues,geoelements numberofatoms,xyz,coords

#Notes/questions of the java code to create the chem pressures

#---phi = Math.acos(x/np.sin(theta))
#---theta=math.acos(z)
#X Y and Z are directly pulled from the data


#Need to split my varible xyz which contains all of x and all of z into distict lists
#for this part
#Ask patrick if there is a easy way to split np arrays by collumn. 



#Scale is in the java script code to change the size keeping 1 now


import numpy as np

def Chemicalpressure(coeffvalues):
    Geo_Cps =[]
    for val in coeffvalues:
        scale=1

        positive_surf  = pv.ParametricEllipsoid(
        xradius = 1,
        yradius = 1,
        zradius = 1,
        )
        negative_surf = pv.ParametricEllipsoid(
        xradius = 1,
        yradius = 1,
        zradius = 1,
        )
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
                positive_surf.points[i] = np.zeros(3)
                negative_surf.points[i]= np.array([ r * np.sin(theta) * np.cos(phi), 
                                                    r * np.sin(theta) * np.sin(phi),
                                                    r * np.cos(theta)])
            else:
                negative_surf.points[i] = np.zeros(3)
                positive_surf.points[i]= np.array([ r * np.sin(theta) * np.cos(phi), 
                                                    r * np.sin(theta) * np.sin(phi),
                                                    r * np.cos(theta)])
            Geo_Cps.append([positive_surf,negative_surf])

    return Geo_Cps

Geo_Cps = Chemicalpressure(coeffvalues)
#print(Geo_Cps)
b = []
for i in Geo_Cps:
    #print(i[0].bounds)
    if i[0].bounds[1]+i[0].bounds[3]+i[0].bounds[5]  >0 or i[1].bounds[1]+i[1].bounds[3]+i[1].bounds[5] >0 :
        b.append(i)
#a = np.where(Geo_Cps[:][0].>0, True, False)

#print(len(b))
#print(b)

    #p=pv.Plotter()
    #p.add_mesh(psurf,style='surface', color='blue')
    #p.add_mesh(nsurf,style='surface', color='red')
    #p.show()

#Graphing the spheres just the atom location
'''
for i, coords in enumerate(xyz):
    Sphere = pv.Sphere()
    Sphere=Sphere.translate(coords, inplace=True)
    if numberofatoms[i] == 0:
        plotter.add_mesh(Sphere, color='blue')
    else:
        plotter.add_mesh(Sphere, color='red')
'''

#Appling the chemical presure lobe
print(len(Geo_Cps))
def apply_CP(geovalues, cellvalues, Geo_Cps, xyz, scale=10,):
    natoms_geo = len(geovalues)
    natoms_coords = len(xyz)
    lobes = []
    counter = 0
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
                        dist = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5
                        if(dist < 0.1): 
                        #if True:
                            Ellipsoid1=Geo_Cps[k1][1]
                            Ellipsoid2=Geo_Cps[k1][0]
                            #print(Geo_Cps[k1][1].bounds[1])
                            Ellipsoid1 = pv.ParametricEllipsoid(xradius=Geo_Cps[k1][1].bounds[1], yradius=Geo_Cps[k1][1].bounds[3], zradius=Geo_Cps[k1][1].bounds[5])  #negative
                            Ellipsoid2 = pv.ParametricEllipsoid(xradius=Geo_Cps[k1][0].bounds[1], yradius=Geo_Cps[k1][0].bounds[3], zradius=Geo_Cps[k1][0].bounds[5])  #positive
                            #NEED TO FIGURE OUT THE CORRECT WAY TO SCALE THE XYZ. This is the problem ^^^^^^^^
                            for i in range(0,1):
                                
                                print(Geo_Cps[k1][1].bounds)
                                print(Geo_Cps[k1][0].bounds)
                                xyztemp=[x2,y2,z2]
                                #print(xyztemp)
                                Ellipsoid1=Ellipsoid1.translate(xyztemp, inplace=True)
                                Ellipsoid2=Ellipsoid2.translate(xyztemp, inplace=True)
                                print(Ellipsoid1)
                                Ellipsoid1.points*=scale
                                Ellipsoid2.points*=scale

                                plotter.add_mesh(Ellipsoid1, style='surface', color='purple')
                                plotter.add_mesh(Ellipsoid2, style='surface', color='green')
                                #plotter.show()
                                #raise(SystemExit)
                                counter=counter+1
                                

CreateCPs = apply_CP(geovalues, cellvalues, Geo_Cps, xyz, scale=10)


plotter.show()
#Currently previous code need to rewrite into python.
'''
function apply_CP(geo_CPs, geo, cell, template) 
  var natoms_geo = geo.length
  var natoms_template = template.length
  var lobes = []
  var lobes_temp = []
  var counter = 0
  var j1,j2,j3,k1,k2,x1,y1,z1,x2,y2,z2,dist
  for( j1=-4; j1 < 5; j1++) {
    for( j2=-4; j2 < 5; j2++) {
      for( j3=-4; j3 < 5; j3++) {
        for(k1 = 0; k1 < natoms_geo; k1++) {
          for(k2 = 0; k2 < natoms_template; k2++) {
            x1=geo[k1][1]+cell[0][0]*j1+cell[1][0]*j2+cell[2][0]*j3
            y1=geo[k1][2]+cell[0][1]*j1+cell[1][1]*j2+cell[2][1]*j3
            z1=geo[k1][3]+cell[0][2]*j1+cell[1][2]*j2+cell[2][2]*j3
            x2=template[k2][1]
            y2=template[k2][2]
            z2=template[k2][3]
            dist = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5
            if(dist < 0.1) {
              lobes_temp[0] = new THREE.Mesh(geo_CPs[k1][0], plus)
              lobes_temp[1] = new THREE.Mesh(geo_CPs[k1][1], minus)
              for (let i=0; i<2; i++) {
                lobes_temp[i].translateX(x2)
                lobes_temp[i].translateY(y2)
                lobes_temp[i].translateZ(z2)
                lobes_temp[i].geometry.verticesNeedUpdate = true
                lobes_temp[i].geometry.facesNeedUpdate = true
                lobes_temp[i].geometry.elementsNeedUpdate = true
                lobes_temp[i].geometry.computeFaceNormals()
                lobes_temp[i].geometry.computeVertexNormals()
                scene.add(lobes_temp[i])
                lobes[counter] = lobes_temp[i]
                counter
              }
            }
          }
        }
      }
    }
  }
  return(lobes)
'''
'''
this code below is also from the html script and it creates the bonds between the atoms
function drawbonds_cylinder(atoms,atom1,color1,atom2,color2,d_min,d_max,radius) {
  var dist = 0.0;
  var d_x,d_y,d_z,theta,phi;
  var new_cylinder = [];
  var counter = 0;
  var nbonds=0;
  for(let j = 0; j < atoms.length ; j++) {
    if(atom1 === atom2) {
      for(let k = 0; k < atoms.length ; k++) {
        dist = ((atoms[j][1]-atoms[k][1])**2.0 + (atoms[j][2]-atoms[k][2])**2.0 + (atoms[j][3]-atoms[k][3])**2.0)**0.5;
        if((dist >= d_min) && (dist <= d_max) && (atoms[j][0] === atom1) && (atoms[k][0] === atom2)) {
          nbonds++;
	        d_x = (atoms[k][1]-atoms[j][1]);
	        d_y = (atoms[k][2]-atoms[j][2]);
	        d_z = (atoms[k][3]-atoms[j][3]);
	        theta=Math.acos(d_z/dist);
	        phi =0.0;
	        if((theta > 0) && (theta < Math.PI)) {
            var acin = Math.round( d_x / Math.sin(theta) / dist * 1000000 + Number.EPSILON ) / 1000000;
            phi = Math.acos(acin);
	        }
	        if(d_y < 0.0) {
	          phi=-phi;
	        }

	        var geometry2 = new THREE.CylinderBufferGeometry(radius,radius,dist,20); //IMAGE QUALITY
	        geometry2.rotateX(Math.PI/2);
	        geometry2.rotateY(theta);
	        geometry2.rotateZ(phi);
	        geometry2.translate(atoms[j][1]+d_x*0.5,atoms[j][2]+d_y*0.5,atoms[j][3]+d_z*0.5);
	        geometry2.computeFaceNormals();
	        geometry2.normalsNeedUpate = true;
	        new_cylinder[counter] = new THREE.Mesh( geometry2,color1);
	        scene.add(new_cylinder[counter]);
	        counter++;
        }
      }
    }
    else {
      for(let k = 0; k < atoms.length ; k++) {
        dist = ((atoms[j][1]-atoms[k][1])**2.0 + (atoms[j][2]-atoms[k][2])**2.0 + (atoms[j][3]-atoms[k][3])**2.0)**0.5;
        if((dist >= d_min) && (dist <= d_max) && (atoms[j][0] === atom1) && (atoms[k][0] === atom2)) {
          d_x = (atoms[k][1]-atoms[j][1]);
          d_y = (atoms[k][2]-atoms[j][2]);
          d_z = (atoms[k][3]-atoms[j][3]);
          theta=Math.acos(d_z/dist);
          phi =0.0;
          if((theta > 0) && (theta < Math.PI)) {
            var acin = Math.round( d_x / Math.sin(theta) / dist * 1000000 + Number.EPSILON ) / 1000000;
            phi = Math.acos(acin);
          }
          if(d_y < 0.0) {
            phi=-phi;
          }
          var geometry2 = new THREE.CylinderBufferGeometry(radius,radius,dist/2.0,20); //IMAGE QUALITY
          geometry2.rotateX(Math.PI/2);
          geometry2.rotateY(theta);
          geometry2.rotateZ(phi);
          geometry2.translate(atoms[j][1]+d_x*0.25,atoms[j][2]+d_y*0.25,atoms[j][3]+d_z*0.25);
          geometry2.computeFaceNormals();
          geometry2.normalsNeedUpate = true;
          new_cylinder[counter] = new THREE.Mesh( geometry2,color1);
          scene.add(new_cylinder[counter]);
          counter++;
       	  geometry2 = new THREE.CylinderBufferGeometry(radius,radius,dist/2.0,20); //IMAGE QUALITY
          geometry2.rotateX(Math.PI/2);
          geometry2.rotateY(theta);
          geometry2.rotateZ(phi);
          geometry2.translate(atoms[j][1]+d_x*0.75,atoms[j][2]+d_y*0.75,atoms[j][3]+d_z*0.75);
          geometry2.computeFaceNormals();
          geometry2.normalsNeedUpate = true;
          new_cylinder[counter] = new THREE.Mesh( geometry2,color2);
          scene.add(new_cylinder[counter]);
          counter++;
        }
      }
    }
  }
  return(new_cylinder);
}
'''

