import math

# the level at which the surface is defined as solid instead of void
ISO_CONTOUR_LEVEL = 5  # 25

# as long as the chunk size equals this any points that can make the surface solid will be within the neighboring 27 cells
CHUNK_SIZE = ISO_CONTOUR_LEVEL

# the size and positioning of the calculated area
SAMPLING_SPACE_SIZE = [100, 100, 100]
SAMPLING_SPACE_OFFSET = [0, 0, 0]

# the maximum fill depth (to avoid any infinite loops in the case of an error in the code, part, or input settings)
MAX_FILL_DEPTH = 1000000

# the number of iterations of surface tension to smooth out the surface (the more the better the surface)
TENSION_ITTERATIONS = 25


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



# a list of all the unorder points of the point cloud
pointCloud = []



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


#"""
newPoints = fibonacci_sphere(175)
for point in newPoints:
    pointCloud.append([point[0]*25+34, point[1]*25+34, point[2]*25+34])
newPoints = fibonacci_sphere(75)
for point in newPoints:
    pointCloud.append([point[0]*11+34, point[1]*11+34, point[2]*11+34])
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


# generating surface normals
normalsGrid = [  # an array of 0's the same size as the signedDistanceField
    [       [0 for z in range(SAMPLING_SPACE_SIZE[2]-2)]
        for y in range(SAMPLING_SPACE_SIZE[1]-2)]
    for x in range(SAMPLING_SPACE_SIZE[0]-2)]

# looping through the entire sample space minus the bounds
for x_ in range(SAMPLING_SPACE_SIZE[0] - 2):
    for y_ in range(SAMPLING_SPACE_SIZE[1] - 2):
        for z_ in range(SAMPLING_SPACE_SIZE[2] - 2):
            x, y, z = x_ + 1, y_ + 1, z_ + 1  # the grid coordinates for the signedDistanceField

            # calculating the slope across various points
            difX = (signedDistanceField[x + 1][y][z] - signedDistanceField[x - 1][y][z]) / 2
            difY = (signedDistanceField[x][y + 1][z] - signedDistanceField[x][y - 1][z]) / 2
            difZ = (signedDistanceField[x][y][z + 1] - signedDistanceField[x][y][z - 1]) / 2

            # calculating the normalized normal
            l = math.sqrt(difX*difX + difY*difY + difZ*difZ)
            if l: normal = [difX / l, difY / l, difZ / l]
            else: normal = [0, 0, 0]

            # filling in the value
            normalsGrid[x_][y_][z_] = normal


print("Normals Calculated")


print("Surface Tension Calculated 0% |                    |", end='\r')

# apply surface tension repeadidly to smooth the surface and remove the lumps
for itteration in range(TENSION_ITTERATIONS):
    print(f"Surface Tension Calculated {round(itteration/TENSION_ITTERATIONS*100)}% |{'='*round(itteration/TENSION_ITTERATIONS*20)}{' '*(20-round(itteration/TENSION_ITTERATIONS*20))}|", end='\r')

    # looping through the sample space except the bounds
    for x_ in range(SAMPLING_SPACE_SIZE[0] - 4):
        for y_ in range(SAMPLING_SPACE_SIZE[1] - 4):
            for z_ in range(SAMPLING_SPACE_SIZE[2] - 4):
                x, y, z = x_ + 2, y_ + 2, z_ + 2  # the grid coordinates for the signedDistanceField

                # calculating the new signed distance based on the surface normal

                # improve this with the other missing terms so it doesn't grow unstable and just create noise
                surfaceTension = (normalsGrid[x][y][z][0] - normalsGrid[x - 2][y][z][0]) / 2 + (normalsGrid[x][y][z][1] - normalsGrid[x][y - 2][z][1]) / 2 + (normalsGrid[x][y][z][2] - normalsGrid[x][y][z - 2][2]) / 2
                signedDistanceField[x][y][z] = signedDistanceField[x][y][z] + surfaceTension*0.025

    # looping through the entire sample space minus the bounds
    for x_ in range(SAMPLING_SPACE_SIZE[0] - 2):
        for y_ in range(SAMPLING_SPACE_SIZE[1] - 2):
            for z_ in range(SAMPLING_SPACE_SIZE[2] - 2):
                x, y, z = x_ + 1, y_ + 1, z_ + 1  # the grid coordinates for the signedDistanceField

                # calculating the slope across various points
                difX = (signedDistanceField[x + 1][y][z] - signedDistanceField[x - 1][y][z]) / 2
                difY = (signedDistanceField[x][y + 1][z] - signedDistanceField[x][y - 1][z]) / 2
                difZ = (signedDistanceField[x][y][z + 1] - signedDistanceField[x][y][z - 1]) / 2

                # calculating the normalized normal
                l = math.sqrt(difX*difX + difY*difY + difZ*difZ)
                if l: normal = [difX / l, difY / l, difZ / l]
                else: normal = [0, 0, 0]

                # filling in the value
                normalsGrid[x_][y_][z_] = normal


print("\x1b[2KSurface Tension Calculated")




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

