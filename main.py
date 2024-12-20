from PIL import Image as im  # importing PIL to create png's

import math

import numpy as np 
import open3d as o3d

import json

# determins if the signed distance field is being imported or being calculated from a pcd
LOADING_DISTANCE_FIELD_SAVE = True
FIELD_FILE = "objectSaveFibSphere.json"
SAVE_FILE = "objectSave.json"

# the level at which the surface is defined as solid instead of void
ISO_CONTOUR_LEVEL = 4  # 25

# as long as the chunk size equals this any points that can make the surface solid will be within the neighboring 27 cells
CHUNK_SIZE = 25#math.ceil(ISO_CONTOUR_LEVEL * 1.5)

# the size and positioning of the calculated area
SAMPLING_SPACE_SIZE = [101, 101, 101]
SAMPLING_SPACE_OFFSET = [0, 0, 0]

# the maximum fill depth (to avoid any infinite loops in the case of an error in the code, part, or input settings)
MAX_FILL_DEPTH = 1000000

# the number of iterations of surface tension to smooth out the surface (the more the better the surface)
TENSION_ITTERATIONS = 225#150

# the delta time and spacing for various components of the simulation (space and time)
DELTA_X = 1
DELTA_Y = 1
DELTA_Z = 1

DT = 0.005



print("Packages Imported and Setup Complete")


# gets the hash for the chunk hash map
def GetChunk(point: tuple) -> tuple:
    return [int(point[0] // CHUNK_SIZE), int(point[1] // CHUNK_SIZE), int(point[2] // CHUNK_SIZE)]


# returns all the surrounding 27 chunks' positions (I'm too lazy to write out all of this everytime it's needed and it's faster to have it hard coded)
def GetChunks(point: tuple) -> list:  # includes the inputed point in this selection
    chunk = GetChunk(point)
    return [
        [chunk[0] - 1, chunk[1] - 1, chunk[2] - 1], [chunk[0], chunk[1] - 1, chunk[2] - 1], [chunk[0] + 1, chunk[1] - 1, chunk[2] - 1],
        [chunk[0] - 1, chunk[1]    , chunk[2] - 1], [chunk[0], chunk[1]    , chunk[2] - 1], [chunk[0] + 1, chunk[1]    , chunk[2] - 1],
        [chunk[0] - 1, chunk[1] + 1, chunk[2] - 1], [chunk[0], chunk[1] + 1, chunk[2] - 1], [chunk[0] + 1, chunk[1] + 1, chunk[2] - 1],
        
        [chunk[0] - 1, chunk[1] - 1, chunk[2]], [chunk[0], chunk[1] - 1, chunk[2]], [chunk[0] + 1, chunk[1] - 1, chunk[2]],
        [chunk[0] - 1, chunk[1]    , chunk[2]], [chunk[0], chunk[1]    , chunk[2]], [chunk[0] + 1, chunk[1]    , chunk[2]],
        [chunk[0] - 1, chunk[1] + 1, chunk[2]], [chunk[0], chunk[1] + 1, chunk[2]], [chunk[0] + 1, chunk[1] + 1, chunk[2]],

        [chunk[0] - 1, chunk[1] - 1, chunk[2] + 1], [chunk[0], chunk[1] - 1, chunk[2] + 1], [chunk[0] + 1, chunk[1] - 1, chunk[2] + 1],
        [chunk[0] - 1, chunk[1]    , chunk[2] + 1], [chunk[0], chunk[1]    , chunk[2] + 1], [chunk[0] + 1, chunk[1]    , chunk[2] + 1],
        [chunk[0] - 1, chunk[1] + 1, chunk[2] + 1], [chunk[0], chunk[1] + 1, chunk[2] + 1], [chunk[0] + 1, chunk[1] + 1, chunk[2] + 1]
    ]


# gets the neighbors to a point
def GetNeighbors(point: tuple) -> list:  # excludes the inputed point in this selection
    return [
        [point[0] - 1, point[1] - 1, point[2] - 1], [point[0] - 1, point[1], point[2] - 1], [point[0] - 1, point[1] + 1, point[2] - 1],
        [point[0]    , point[1] - 1, point[2] - 1], [point[0]    , point[1], point[2] - 1], [point[0]    , point[1] + 1, point[2] - 1],
        [point[0] + 1, point[1] - 1, point[2] - 1], [point[0] + 1, point[1], point[2] - 1], [point[0] + 1, point[1] + 1, point[2] - 1],

        [point[0] - 1, point[1] - 1, point[2]], [point[0] - 1, point[1], point[2]], [point[0] - 1, point[1] + 1, point[2]],
        [point[0]    , point[1] - 1, point[2]]                                    , [point[0]    , point[1] + 1, point[2]],
        [point[0] + 1, point[1] - 1, point[2]], [point[0] + 1, point[1], point[2]], [point[0] + 1, point[1] + 1, point[2]],

        [point[0] - 1, point[1] - 1, point[2] + 1], [point[0] - 1, point[1], point[2] + 1], [point[0] - 1, point[1] + 1, point[2] + 1],
        [point[0]    , point[1] - 1, point[2] + 1], [point[0]    , point[1], point[2] + 1], [point[0]    , point[1] + 1, point[2] + 1],
        [point[0] + 1, point[1] - 1, point[2] + 1], [point[0] + 1, point[1], point[2] + 1], [point[0] + 1, point[1] + 1, point[2] + 1]
    ]


# uses a flood fill algerithm to systematically flip the signs of an internal void/solid in a structure
def FloodFill(point: tuple, distanceField: list, signedGrid: list, sign: int) -> None:
    # getting the current sign as a reference for where to search
    signedGrid[point[0]][point[1]][point[2]] = sign
    neighbors = GetNeighbors(point)

    # looping through repeatidly going through all the neighbors
    for itteration in range(MAX_FILL_DEPTH):
        # a list of the new neighbors from the new searches
        newNeighbors = []

        # going through all neighboring points
        for neighboringPoint in neighbors:
            # stopping the flood fill from sampling points outside the grid
            try:
            
                # making sure the point hasn't been sampled yet (making it equal to 0/false/null)
                if not signedGrid[neighboringPoint[0]][neighboringPoint[1]][neighboringPoint[2]]:
                    # checking if the point is outside the boundary of a surface
                    if distanceField[neighboringPoint[0]][neighboringPoint[1]][neighboringPoint[2]] > ISO_CONTOUR_LEVEL:
                        
                        # updating the sign and getting all neighbors to search in the next itteration
                        signedGrid[neighboringPoint[0]][neighboringPoint[1]][neighboringPoint[2]] = sign
                        newNeighbors += GetNeighbors(neighboringPoint)
            
            # an exception for when a point is sampled outside of the boundary of the distance field
            except IndexError:
                pass
        
        # updating the neighoring points to the new selection generated
        neighbors = newNeighbors
        
        # ending the search once all points have been searched
        if not neighbors: return


# loading a pcd file
def LoadPCD(file: str) -> list:
    pcd = o3d.io.read_point_cloud(file)
    return np.asarray(pcd.points)


# saves a json file
def SaveJson(dictionary: dict, fileName: str) -> None:
    # generating the json file
    jsonFile = json.dumps(dictionary, indent=4)

    # saving the file to memory
    with open(fileName, "w") as outfile:
        outfile.write(jsonFile)


# loads a json file
def LoadJson(fileName: str) -> dict:
    return json.load(open(fileName))  # loading and return the file


# generates a png image
def ImageFromArray(list: list, name: str) -> im.fromarray:  # converts a 2d array of colors into a png image
    img_w, img_h = [len(list), len(list[0])]
    data = np.zeros((img_h, img_w, 3), dtype=np.uint8)
    data[100, 100] = [255, 0, 0]
    
    try:
        list[0][0].rgb
        for x in range(img_w):
            for y in range(img_h):
                data[y][x] = list[x][y].rgb
    except AttributeError:
        for x in range(img_w):
            for y in range(img_h):
                data[y][x] = list[x][y]
    
    img = im.fromarray(data, 'RGB')
    img.save(name)
    
    return img


# creates a fibonacci sphere
def fibonacci_sphere(samples: int):

    points = []
    phi = math.pi * (math.sqrt(5.) - 1.)  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((x, y, z))

    return points


# loading or calculating the distance field
if LOADING_DISTANCE_FIELD_SAVE:
    signedDistanceField = LoadJson(FIELD_FILE)["Signed Distance Field"]
else:
    # a list of all the unorder points of the point cloud
    lpointCloud = LoadPCD("test2.pcd")
    pointCloud = []
    for i in range(len(pointCloud)):
        pointCloud[i] = [pointCloud[i][0] / DELTA_X, pointCloud[i][1] / DELTA_Y, pointCloud[2] / DELTA_Z]

    """
    for p in lpointCloud:
        pointCloud.append([
            (p[0] * 0.5 + 0.25) * 90 + 50,
            (p[1] * 0.5 + 0.25) * 90 + 50,
            (p[2] * 0.5 + 0.25) * 90 + 50
        ])#"""


    #"""
    newPoints = fibonacci_sphere(575)
    for point in newPoints:
        pointCloud.append([point[0] * 34 + 50, point[1] * 34 + 50, point[2] * 34 + 50])
        #pointCloud.append([point[0]*25+34, point[1]*25+34, point[2]*25+34])
    newPoints = fibonacci_sphere(250)
    for point in newPoints:
        pointCloud.append([point[0] * 15 + 50, point[1] * 15 + 50, point[2] * 15 + 50])
        #pointCloud.append([point[0]*11+34, point[1]*11+34, point[2]*11+34])
    #"""

    """
    newPoints = fibonacci_sphere(10)
    for point in newPoints:
        pointCloud.append([point[0]*14+50, point[1]*14+50, point[2]*14+50])
    #"""

    """
    # generating a test cube for the point cloud
    for x in range(8):
        for y in range(8):
            for z in range(8):
                if 7 in [x, y, z] or 0 in [x, y, z]:
                    if [x, y, z] not in pointCloud: pointCloud.append([x*4+20, y*4+20, z*4+20])
    #"""



    # the hash map for sorting the points into chunks (for faster data calculations)
    chunkedPoints = {}

    # moving through the points and sorting them into chunks
    for point in pointCloud:
        chunk = str(GetChunk(point))
        if chunk in chunkedPoints:  # adding the point to an exsisting chunk
            chunkedPoints[chunk].append(point)
        else:
            chunkedPoints[chunk] = [point]  # creating a new chunk


    print("Chunking Complete")


    # a 3d structured space/field containing the calculed signed distance to the nearest point
    signedDistanceField = []

    # calculating the minimum unsigned distance for every point in the structure space/field
    print("Distance Field 0% |                    |", end='\r')
    for x_ in range(SAMPLING_SPACE_SIZE[0]):
        print(f"Distance Field {round(x_/SAMPLING_SPACE_SIZE[0]*100)}% |{'='*round(x_/SAMPLING_SPACE_SIZE[0]*20)}{' '*(20-round(x_/SAMPLING_SPACE_SIZE[0]*20))}|", end='\r')
        yLayer = []
        for y_ in range(SAMPLING_SPACE_SIZE[1]):
            zLayer = []
            for z_ in range(SAMPLING_SPACE_SIZE[2]):
                # getting the global position of the point
                x, y, z = x_ + SAMPLING_SPACE_OFFSET[0], y_ + SAMPLING_SPACE_OFFSET[1], z_ + SAMPLING_SPACE_OFFSET[2]
                
                neighboringChunks = GetChunks([x, y, z])
                nearest = ISO_CONTOUR_LEVEL + 9999  # chosing a value out of range of the surface to not conflict with it
                
                # searching through the neighboring points to find the nearest one
                for chunk_ in neighboringChunks:
                    chunk = str(chunk_)
                    if chunk in chunkedPoints:
                        for point in chunkedPoints[chunk]:
                            newPoint = [x-point[0], y-point[1], z-point[2]]
                            unsignedDistance = math.sqrt(newPoint[0]*newPoint[0] + newPoint[1]*newPoint[1] + newPoint[2]*newPoint[2])
                            if unsignedDistance < nearest: nearest = unsignedDistance  # updating the nearest point

                # adding the signed distance to the nearest point
                zLayer.append(nearest)
            
            # stacking the layers to create a 3d space
            yLayer.append(zLayer)
        signedDistanceField.append(yLayer)


    print("\x1b[2KDistance Field Calculated")


    # saving the json of the output
    SaveJson({
        "Signed Distance Field": signedDistanceField
    }, SAVE_FILE)
    print("Json File Saved")




# temporary
for x in range(SAMPLING_SPACE_SIZE[0]):
    for y in range(SAMPLING_SPACE_SIZE[1]):
        for z in range(SAMPLING_SPACE_SIZE[2]):
            signedDistanceField[x][y][z] = abs(signedDistanceField[x][y][z])





print("Surface Tension Calculated 0% |                    |", end='\r')

# generating surface normals
normalsGrid = [  # an array of 0's the same size as the signedDistanceField
    [       [0 for z in range(SAMPLING_SPACE_SIZE[2]-2)]
        for y in range(SAMPLING_SPACE_SIZE[1]-2)]
    for x in range(SAMPLING_SPACE_SIZE[0]-2)]

# apply surface tension repeadidly to smooth the surface and remove the lumps
for itteration in range(TENSION_ITTERATIONS):
    print(f"Surface Tension Calculated {round(itteration/TENSION_ITTERATIONS*100)}% |{'='*round(itteration/TENSION_ITTERATIONS*20)}{' '*(20-round(itteration/TENSION_ITTERATIONS*20))}|", end='\r')

    
    # looping through the entire sample space minus the bounds
    for x_ in range(SAMPLING_SPACE_SIZE[0] - 2):
        for y_ in range(SAMPLING_SPACE_SIZE[1] - 2):
            for z_ in range(SAMPLING_SPACE_SIZE[2] - 2):
                x, y, z = x_ + 1, y_ + 1, z_ + 1  # the grid coordinates for the signedDistanceField

                # calculating the slope across various points
                difX = (signedDistanceField[x + 1][y][z] - signedDistanceField[x - 1][y][z]) / (2 * DELTA_X)
                difY = (signedDistanceField[x][y + 1][z] - signedDistanceField[x][y - 1][z]) / (2 * DELTA_Y)
                difZ = (signedDistanceField[x][y][z + 1] - signedDistanceField[x][y][z - 1]) / (2 * DELTA_Z)

                # calculating the normalized normal
                l = math.sqrt(difX*difX + difY*difY + difZ*difZ)
                if l: normal = [difX / l, difY / l, difZ / l]
                else: normal = [0, 0, 0]

                # filling in the value
                normalsGrid[x_][y_][z_] = normal


    # copying the entirety of the signed distance field
    signedDistanceFieldCopy = [
    [       [signedDistanceField[x][y][z] for z in range(SAMPLING_SPACE_SIZE[2])]
        for y in range(SAMPLING_SPACE_SIZE[1])]
    for x in range(SAMPLING_SPACE_SIZE[0])]

    # looping through the sample space except the bounds
    for x_ in range(SAMPLING_SPACE_SIZE[0] - 5):
        for y_ in range(SAMPLING_SPACE_SIZE[1] - 5):
            for z_ in range(SAMPLING_SPACE_SIZE[2] - 5):
                x, y, z = x_ + 3, y_ + 3, z_ + 3  # the grid coordinates for the signedDistanceField
                x2, y2, z2 = x - 1, y - 1, z - 1

                # calculating the surface tension
                kappa = ((normalsGrid[x2][y2][z2][0] - normalsGrid[x2 - 2][y2][z2][0]) / (2 * DELTA_X) + (normalsGrid[x2][y2][z2][1] - normalsGrid[x2][y2 - 2][z2][1]) / (2 * DELTA_Y) + (normalsGrid[x2][y2][z2][2] - normalsGrid[x2][y2][z2 - 2][2]) / (2 * DELTA_Z))
                
                dMinusX = (signedDistanceField[x][y][z] - signedDistanceField[x - 1][y][z]) / DELTA_X
                dMinusY = (signedDistanceField[x][y][z] - signedDistanceField[x][y - 1][z]) / DELTA_Y
                dMinusZ = (signedDistanceField[x][y][z] - signedDistanceField[x][y][z - 1]) / DELTA_Z

                dPlusX = (signedDistanceField[x + 1][y][z] - signedDistanceField[x][y][z]) / DELTA_X
                dPlusY = (signedDistanceField[x][y + 1][z] - signedDistanceField[x][y][z]) / DELTA_Y
                dPlusZ = (signedDistanceField[x][y][z + 1] - signedDistanceField[x][y][z]) / DELTA_Z
                
                gradPlus = math.sqrt(
                    pow(max(dMinusX, 0.0), 2.0) + pow(min(dPlusX, 0.0), 2.0) +
                    pow(max(dMinusY, 0.0), 2.0) + pow(min(dPlusY, 0.0), 2.0) +
                    pow(max(dMinusZ, 0.0), 2.0) + pow(min(dPlusZ, 0.0), 2.0)
                )
                gradMinus = math.sqrt(
                    pow(min(dMinusX, 0.0), 2.0) + pow(max(dPlusX, 0.0), 2.0) +
                    pow(min(dMinusY, 0.0), 2.0) + pow(max(dPlusY, 0.0), 2.0) +
                    pow(min(dMinusZ, 0.0), 2.0) + pow(max(dPlusZ, 0.0), 2.0)
                )
                
                # applying the surface tension
                signedDistanceFieldCopy[x][y][z] += DT * (max(kappa, 0) * gradPlus + min(kappa, 0) * gradMinus)

    # replacing the old signed distance field with the new one
    signedDistanceField = signedDistanceFieldCopy


print("\x1b[2KSurface Tension Calculated")



# generating a new grid to keep track of the sign changes
signedGrid = [  # an array of 0's the same size as the signedDistanceField
    [       [0 for z in range(SAMPLING_SPACE_SIZE[2])]
        for y in range(SAMPLING_SPACE_SIZE[1])]
    for x in range(SAMPLING_SPACE_SIZE[0])]

# filling the inital void (the entire perimeter should be void if the program's setup was correctly done and set to encompass the entire point cloud)
print("Signs Calculated 0% |                    |", end='\r')
FloodFill([0, 0, 0], signedDistanceField, signedGrid, 1)  # this is the bottle neck

# starting a march across the grid to fill in all reigons
for x in range(SAMPLING_SPACE_SIZE[0]):
    print(f"Signs Calculated {round(x/SAMPLING_SPACE_SIZE[0]*100)}% |{'='*round(x/SAMPLING_SPACE_SIZE[0]*20)}{' '*(20-round(x/SAMPLING_SPACE_SIZE[0]*20))}|", end='\r')
    for y in range(SAMPLING_SPACE_SIZE[1]):
        # beginning the march
        sign = 1  # the entire perimeter should be void so it's safe to assume it as such
        for z in range(SAMPLING_SPACE_SIZE[2]):
            # checking if the point has already been filled
            if signedGrid[x][y][z]:
                # setting the sign to the already calculated value and moving on
                sign = signedGrid[x][y][z]
            else:
                # checking if the point's outside an object
                if signedDistanceField[x][y][z] > ISO_CONTOUR_LEVEL:
                    # flood filling the area
                    sign *= -1  # flipping the sign (this should only run once it's left the bounds of the object's surface)
                    FloodFill([x, y, z], signedDistanceField, signedGrid, sign)


# adjusting the signedDistanceField for the generated signs
for x in range(SAMPLING_SPACE_SIZE[0]):
    for y in range(SAMPLING_SPACE_SIZE[1]):
        for z in range(SAMPLING_SPACE_SIZE[2]):
            # adjusting the sign
            if signedGrid[x][y][z]:
                signedDistanceField[x][y][z] *= signedGrid[x][y][z]


print("\x1b[2KSigns Calculated")



# rendering a slice of the signed distance field
selection = []
for x in range(SAMPLING_SPACE_SIZE[0]):
    layer = []
    for y in range(SAMPLING_SPACE_SIZE[1]):
        val = signedDistanceField[x][y][SAMPLING_SPACE_SIZE[2]//2]
        v = min(max(   val * 10 + 60   , 0), 255)
        v = [v, v, v]
        if val < ISO_CONTOUR_LEVEL: v = [255, 0, v[0]]
        layer.append(v)
    selection.append(layer)

ImageFromArray(selection, "output.png")



# visulizing the outputed data

import skimage
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Create a figure and 3D axes
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Extract the isosurface
verts, faces, normals, values = skimage.measure.marching_cubes(np.array(signedDistanceField), ISO_CONTOUR_LEVEL)

# Plot the isosurface
ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2], cmap='jet', alpha=0.5)

plt.show()

