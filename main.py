from PIL import Image as im  # importing PIL to create png's

import math

import numpy as np 
import open3d as o3d

import json

# determins if the signed distance field is being imported or being calculated from a pcd
LOADING_DISTANCE_FIELD_SAVE = False
GENERATE_NEW_VERSION_SDF = True
USING_OLD_SDF = False
FIELD_FILE = "tensioned_objectSaveShelledGood.json"  # "objectSaveFibSphere.json"  # tensioned_      "tensioned_objectSaveNew.json"#"objectSaveNew.json"    #"tensioned_objectSaveShelled.json"#""
SAVE_FILE = "objectSaveShelledTesting2.json"

TENSION_IMAGES_SAVE_PATH = "ShrinkWrap/"
SAVE_TENSION_IMAGES = True

# the level at which the surface is defined as solid instead of void
ISO_CONTOUR_LEVEL = 4

# as long as the chunk size equals this any points that can make the surface solid will be within the neighboring 27 cells
CHUNK_SIZE = ISO_CONTOUR_LEVEL#math.ceil(ISO_CONTOUR_LEVEL * 1.5)
SHELL_CHUNK_SIZE = 8#25  # the chunk size for the points generated to represent the bubbly shell around the approximated objects

# the size and positioning of the calculated area
SAMPLING_SPACE_SIZE = [101, 101, 101]
SAMPLING_SPACE_OFFSET = [0, 0, 0]

# the maximum fill depth (to avoid any infinite loops in the case of an error in the code, part, or input settings)
MAX_FILL_DEPTH = 1000000

# the number of iterations of surface tension to smooth out the surface (the more the better the surface)
SIMULATING_SURFACE_TENSION = False
TENSION_ITTERATIONS = 150#50
VELOCITY = -2.5  # constant velocity with respect to the gradient

# the delta time and spacing for various components of the simulation (space and time)
DELTA_X = 1
DELTA_Y = 1
DELTA_Z = 1

DT = 0.0055



print("Packages Imported and Setup Complete")


# gets the hash for the chunk hash map
def GetChunk(point: tuple, cs: float) -> tuple:
    return [int(point[0] // cs), int(point[1] // cs), int(point[2] // cs)]


# returns all the surrounding 27 chunks' positions (I'm too lazy to write out all of this everytime it's needed and it's faster to have it hard coded)
def GetChunks(point: tuple, cs: float) -> list:  # includes the inputed point in this selection
    chunk = GetChunk(point, cs)
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
                    if distanceField[neighboringPoint[0]][neighboringPoint[1]][neighboringPoint[2]] >= ISO_CONTOUR_LEVEL:
                        
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


# calculates the internal/external sections of parts
def CalculateSigns(signedGridInput: list) -> None:
    # filling the inital void (the entire perimeter should be void if the program's setup was correctly done and set to encompass the entire point cloud)
    print("Signs Calculated 0% |                    |", end='\r')
    FloodFill([0, 0, 0], signedDistanceField, signedGridInput, 1)  # this is the bottle neck

    # starting a march across the grid to fill in all reigons
    for x in range(SAMPLING_SPACE_SIZE[0]):
        print(f"Signs Calculated {round(x/SAMPLING_SPACE_SIZE[0]*100)}% |{'='*round(x/SAMPLING_SPACE_SIZE[0]*20)}{' '*(20-round(x/SAMPLING_SPACE_SIZE[0]*20))}|", end='\r')
        for y in range(SAMPLING_SPACE_SIZE[1]):
            # beginning the march
            sign = 1  # the entire perimeter should be void so it's safe to assume it as such
            for z in range(SAMPLING_SPACE_SIZE[2]):
                # checking if the point has already been filled
                if signedGridInput[x][y][z]:
                    # setting the sign to the already calculated value and moving on
                    sign = signedGridInput[x][y][z]
                else:
                    # checking if the point's outside an object
                    if signedDistanceField[x][y][z] >= ISO_CONTOUR_LEVEL:
                        # flood filling the area
                        sign *= -1  # flipping the sign (this should only run once it's left the bounds of the object's surface)
                        FloodFill([x, y, z], signedDistanceField, signedGridInput, sign)

    print("\x1b[2KSigns Calculated")


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
        chunk = str(GetChunk(point, CHUNK_SIZE))
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
                
                neighboringChunks = GetChunks([x, y, z], CHUNK_SIZE)
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



if USING_OLD_SDF:
    # ensuring the signs of any sdf's are corrected to not mess up any of the following calculations
    for x in range(SAMPLING_SPACE_SIZE[0]):
        for y in range(SAMPLING_SPACE_SIZE[1]):
            for z in range(SAMPLING_SPACE_SIZE[2]):
                signedDistanceField[x][y][z] = abs(signedDistanceField[x][y][z])




# stores and sorts a set of points into a grid of chunks
class ChunkGrid:

    # initializing the chunk grid
    def __init__(self, chunks: list, smallBound: list, largeBound: list, ChunkKey) -> None:
        self.chunks = chunks
        self.smallBound = smallBound
        self.largeBound = largeBound
        self.ChunkKey = ChunkKey

        self.chunkBounds = [len(self.chunks), len(self.chunks[0]), len(self.chunks[0][0])]
    
    # finds the nearest point given a starting position
    def FindNearestPoint(self, samplePosition: list) -> list:  # returns the nearest point and the distance to it ([pointX, pointY, pointZ], dst)
        # bringing the point to the start of the grid (saves time searching through empty space)
        correctedPosition = [
            max(min(samplePosition[0], self.largeBound[0]), self.smallBound[0]),
            max(min(samplePosition[1], self.largeBound[1]), self.smallBound[1]),
            max(min(samplePosition[2], self.largeBound[2]), self.smallBound[2])
        ]

        # getting the starting chunk position
        startingChunkPosition = self.ChunkKey(correctedPosition)
        occupied, finalRadius = self.ExpandingSearch(startingChunkPosition, 0)
        if not occupied: return [(999999, 999999, 999999), 999999]  # ending the search if no points were found
        
        # finding the nearest point within the current search
        nearestPoint = [0, 0, 0]
        nearestDepth = 999999

        # going through the points
        for point in occupied:
            # calculating the distance
            dx, dy, dz = samplePosition[0] - point[0], samplePosition[1] - point[1], samplePosition[2] - point[2]
            dst = dx*dx + dy*dy + dz*dz

            # checking if this is the smallest distance found so far
            if dst < nearestDepth:
                # updating the information on the nearest point
                nearestDepth = dst
                nearestPoint = point
        
        # checking if the radius is grater than the minimum size of the box search
        if finalRadius < nearestDepth:
            # doing an extended search in case there is a closer point within this range
            maxExtendedDepth = math.ceil(finalRadius - nearestDepth)  # the max possible radius out

            # doing the extended search
            newOccupied, newFinalRadius = self.ExpandingSearch(startingChunkPosition, finalRadius + 1, maxSearchDepth=maxExtendedDepth)

            # checking if any of these new points are closer (if any)
            for point in newOccupied:
                # calculating the distance
                dx, dy, dz = samplePosition[0] - point[0], samplePosition[1] - point[1], samplePosition[2] - point[2]
                dst = dx*dx + dy*dy + dz*dz

                # checking if this is the smallest distance found so far
                if dst < nearestDepth:
                    # updating the information on the nearest point
                    nearestDepth = dst
                    nearestPoint = point
        
        # returning the conclusions of the search
        return (nearestPoint[0], nearestDepth)
    
    # expands a search until a chunk with points is found
    def ExpandingSearch(self, startingPos: list, radius: int, maxSearchDepth: int = 9999) -> list:  # (returns a list with all the chunk positions that have points, returns the final radius)
        # searching for any occupied chunks within the search area
        occupiedChunks = []
        for x in range(-radius, radius + 1):
            # making sure the position is valid
            if startingPos[0] + x < 0 or startingPos[0] + x >= self.chunkBounds[0]: continue
            
            # searching along the y axis
            for y in range(-radius, radius + 1):
                # making sure the position is valid
                if startingPos[1] + y < 0 or startingPos[1] + y >= self.chunkBounds[1]: continue
                
                # checking each chunk for occupancy
                if startingPos[2] - radius >= 0:  # startingPos[2] - radius < self.chunkBounds[2] and 
                    if self.chunks[startingPos[0] + x][startingPos[1] + y][startingPos[2] - radius]: occupiedChunks += self.chunks[startingPos[0] + x][startingPos[1] + y][startingPos[2] - radius]
                if startingPos[2] + radius < self.chunkBounds[2]:  # and startingPos[2] + radius >= 0
                    if self.chunks[startingPos[0] + x][startingPos[1] + y][startingPos[2] + radius]: occupiedChunks += self.chunks[startingPos[0] + x][startingPos[1] + y][startingPos[2] + radius]

            # searching along the z axis:
            for z in range(-radius, radius + 1):
                # making sure the position is valid
                if startingPos[2] + z < 0 or startingPos[2] + z >= self.chunkBounds[2]: continue

                # checking each chunk for occupancy
                if startingPos[1] - radius >= 0:  # startingPos[1] - radius < self.chunkBounds[1] and 
                    if self.chunks[startingPos[0] + x][startingPos[1] - radius][startingPos[2] + z]: occupiedChunks += self.chunks[startingPos[0] + x][startingPos[1] - radius][startingPos[2] + z]
                if startingPos[1] + radius < self.chunkBounds[1]:  # and startingPos[1] + radius >= 0
                    if self.chunks[startingPos[0] + x][startingPos[1] + radius][startingPos[2] + z]: occupiedChunks += self.chunks[startingPos[0] + x][startingPos[1] + radius][startingPos[2] + z]

        # searching along the y and z axis
        for y in range(-radius, radius + 1):
            # making sure the position is valid
            if startingPos[1] + y < 0 or startingPos[1] + y >= self.chunkBounds[1]: continue

            # searching along the z axis
            for z in range(-radius, radius + 1):
                # making sure the position is valid
                if startingPos[2] + z < 0 or startingPos[2] + z >= self.chunkBounds[2]: continue

                # checking each chunk for occupancy
                if startingPos[0] - radius >= 0:  # startingPos[0] - radius < self.chunkBounds[0] and 
                    if self.chunks[startingPos[0] - radius][startingPos[1] + y][startingPos[2] + z]: occupiedChunks += self.chunks[startingPos[0] - radius][startingPos[1] + y][startingPos[2] + z]
                if startingPos[0] + radius < self.chunkBounds[0]:  # and startingPos[0] + radius >= 0
                    if self.chunks[startingPos[0] + radius][startingPos[1] + y][startingPos[2] + z]: occupiedChunks += self.chunks[startingPos[0] + radius][startingPos[1] + y][startingPos[2] + z]

        # returning the found occupied chunks    or return the empty list if the search depth limit has been reached
        if occupiedChunks or radius == maxSearchDepth: return (occupiedChunks, radius)

        # continung the search for occupied chunks
        return self.ExpandingSearch(startingPos, radius + 1)  # expanding the search recursively


# generates a sorted chunk grid from a list of points
def GenerateChunks(points: list, chunkSize: float) -> ChunkGrid:
    smallest = [999999, 999999, 999999]
    largest = [-999999, -999999, -999999]

    # finding the smallest and largest points to fit the bounding box perfectly around them
    for point in points:
        smallest = [min(smallest[0], point[0]), min(smallest[1], point[1]), min(smallest[2], point[2])]
        largest = [max(largest[0], point[0]), max(largest[1], point[1]), max(largest[2], point[2])]
    
    # a function to generate the chunk key
    ChunkKey = lambda pos: [
        int((pos[0] - smallest[0]) // chunkSize),  # translating the point than scaling it to fit the local space of the chunk grid
        int((pos[1] - smallest[1]) // chunkSize),
        int((pos[2] - smallest[2]) // chunkSize)
    ]

    # finding the total size of the chunk grid
    gridSize = [
        math.ceil(abs(largest[0] - smallest[0]) / chunkSize),
        math.ceil(abs(largest[1] - smallest[1]) / chunkSize),
        math.ceil(abs(largest[2] - smallest[2]) / chunkSize)
    ]

    # creating the grid for the chunks
    chunkGrid = [
        [       [[] for z in range(gridSize[2])]
            for y in range(gridSize[1])]
        for x in range(gridSize[0])]
    
    # adding the points to the grid
    for point in points:
        # finding the chunk key and adding the point
        key = ChunkKey(point)
        chunkGrid[key[0]][key[1]][key[2]].append([point[0], point[1], point[2]])
    
    # returning the results
    return ChunkGrid(
        chunkGrid,
        smallest,
        largest,
        ChunkKey
    )



# generating an updated sdf
if GENERATE_NEW_VERSION_SDF:
    # generating a shell around the objects

    # generating a new grid to keep track of the sign changes
    signedGrid = [  # an array of 0's the same size as the signedDistanceField
        [       [0 for z in range(SAMPLING_SPACE_SIZE[2])]
            for y in range(SAMPLING_SPACE_SIZE[1])]
        for x in range(SAMPLING_SPACE_SIZE[0])]

    CalculateSigns(signedGrid)

    # all the surface points which will be sorted later
    surfacePoints = []

    # creating the shell around all objects surfaces (still bubbly)
    print("Calculating Object Shell 0% |                    |", end='\r')
    for x in range(1, SAMPLING_SPACE_SIZE[0] - 1):
        print(f"Calculating Object Shell {round(x/SAMPLING_SPACE_SIZE[0]*100)}% |{'='*round(x/SAMPLING_SPACE_SIZE[0]*20)}{' '*(20-round(x/SAMPLING_SPACE_SIZE[0]*20))}|", end='\r')
        for y in range(1, SAMPLING_SPACE_SIZE[1] - 1):
            for z in range(1, SAMPLING_SPACE_SIZE[2] - 1):
                # calculating the new distance from the shell of the object
                if max([
                        max(signedGrid[x - 1][y][z], signedGrid[x + 1][y][z]),
                        max(signedGrid[x][y - 1][z], signedGrid[x][y + 1][z]),
                        max(signedGrid[x][y][z - 1], signedGrid[x][y][z + 1])
                    ]) > 0:
                    if (signedGrid[x][y][z] < 0) or (signedDistanceField[x][y][z] < ISO_CONTOUR_LEVEL):
                        # finding the chunk
                        surfacePoints.append([x, y, z])  # adding the point

    print("\x1b[2KCalculated Object Shell")


    print("Calculating New Distance Field 0% |                    |", end='\r')

    # creating a new signed distance field that uses the distance to the calculated approximate surface
    newSignedDistanceField = [  # an array of 0's the same size as the signedDistanceField
        [       [0 for z in range(SAMPLING_SPACE_SIZE[2])]
            for y in range(SAMPLING_SPACE_SIZE[1])]
        for x in range(SAMPLING_SPACE_SIZE[0])]
    
    # generating the chunk grid and sorting the points
    chunkGrid = GenerateChunks(surfacePoints, SHELL_CHUNK_SIZE)

    # finding all the new distances
    for x in range(SAMPLING_SPACE_SIZE[0]):
        print(f"Calculating New Distance Field {round(x/SAMPLING_SPACE_SIZE[0]*100)}% |{'='*round(x/SAMPLING_SPACE_SIZE[0]*20)}{' '*(20-round(x/SAMPLING_SPACE_SIZE[0]*20))}|", end='\r')
        for y in range(SAMPLING_SPACE_SIZE[1]):
            for z in range(SAMPLING_SPACE_SIZE[2]):
                # finding the nearest point
                nearestPoint, distance = chunkGrid.FindNearestPoint([x, y, z])

                # finding the correct sign for the current position
                sign = signedGrid[x][y][z]
                if signedDistanceField[x][y][z] < ISO_CONTOUR_LEVEL and sign > -1: sign = -1
                if not sign: sign = 1

                # updating the new signed distance field with the new signed-distance
                newSignedDistanceField[x][y][z] = math.sqrt(distance) * sign + ISO_CONTOUR_LEVEL

    print("\x1b[2KCalculated New Distance Field")


    # uploading a new output image
    selection = []
    for x in range(SAMPLING_SPACE_SIZE[0]):
        layer = []
        for y in range(SAMPLING_SPACE_SIZE[1]):
            val = min(max(newSignedDistanceField[x][y][SAMPLING_SPACE_SIZE[2]//2] * 100 + 50, 0), 255)
            
            layer.append(val)
        selection.append(layer)

    ImageFromArray(selection, "new_output.png")

    print("New Signed Distance Field Created")
    signedDistanceField = newSignedDistanceField  # temporary till I update the rest of the code
    #"""



# simulating surface tension
if SIMULATING_SURFACE_TENSION:
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
        for x_ in range(SAMPLING_SPACE_SIZE[0] - 4):
            for y_ in range(SAMPLING_SPACE_SIZE[1] - 4):
                for z_ in range(SAMPLING_SPACE_SIZE[2] - 4):
                    x, y, z = x_ + 2, y_ + 2, z_ + 2  # the grid coordinates for the signedDistanceField
                    x2, y2, z2 = x - 1, y - 1, z - 1  # the grid coordinates for the normals

                    # calculating the surface tension
                    
                    kappa = (
                        (normalsGrid[x2 + 1][y2][z2][0] - normalsGrid[x2 - 1][y2][z2][0]) / (2 * DELTA_X) +
                        (normalsGrid[x2][y2 + 1][z2][1] - normalsGrid[x2][y2 - 1][z2][1]) / (2 * DELTA_Y) +
                        (normalsGrid[x2][y2][z2 + 1][2] - normalsGrid[x2][y2][z2 - 1][2]) / (2 * DELTA_Z)
                    )
                    
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
                    signedDistanceFieldCopy[x][y][z] += DT * (max(kappa + VELOCITY, 0) * gradPlus + min(kappa + VELOCITY, 0) * gradMinus)

        # replacing the old signed distance field with the new one
        signedDistanceField = signedDistanceFieldCopy

        if SAVE_TENSION_IMAGES:
            # rendering a slice of the signed distance field
            selection = []
            for x in range(SAMPLING_SPACE_SIZE[0]):
                layer = []
                for y in range(SAMPLING_SPACE_SIZE[1]):
                    val = signedDistanceField[x][y][SAMPLING_SPACE_SIZE[2]//2]
                    v = min(max(   val * 10 + 60   , 0), 255)
                    v = [v, v, v]
                    #if val < ISO_CONTOUR_LEVEL: v = [255, 0, v[0]]
                    layer.append(v)
                selection.append(layer)

            ImageFromArray(selection, f"{TENSION_IMAGES_SAVE_PATH}output{itteration}.png")
    

    print("\x1b[2KSurface Tension Calculated")



# updating the signs for older-type sdf's
if not GENERATE_NEW_VERSION_SDF and USING_OLD_SDF:
    # generating a new grid to keep track of the sign changes
    signedGrid = [  # an array of 0's the same size as the signedDistanceField
        [       [0 for z in range(SAMPLING_SPACE_SIZE[2])]
            for y in range(SAMPLING_SPACE_SIZE[1])]
        for x in range(SAMPLING_SPACE_SIZE[0])]

    CalculateSigns(signedGrid)

    # adjusting the signedDistanceField for the generated signs
    for x in range(SAMPLING_SPACE_SIZE[0]):
        for y in range(SAMPLING_SPACE_SIZE[1]):
            for z in range(SAMPLING_SPACE_SIZE[2]):
                # adjusting the sign
                if signedGrid[x][y][z]:
                    signedDistanceField[x][y][z] *= signedGrid[x][y][z]



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


# saving the json of the output
SaveJson({
    "Signed Distance Field": signedDistanceField
}, f"tensioned_{SAVE_FILE}")
print("Json File Saved")


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

