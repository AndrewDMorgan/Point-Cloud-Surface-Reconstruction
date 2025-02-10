#include <iostream>
#include "matar.h"

using namespace mtr;

const bool LOADING_DISTANCE_FIELD_SAVE = false;
const bool GENERATE_NEW_VERSION_SDF = true;
const bool USING_OLD_SDF = true;

const std::string FIELD_FILE = "tensioned_StandFordBunny_untensioned_workingFixed!.json";  // "tensioned_StandFordBunny_untensioned_test2.json"
const std::string SAVE_FILE = "StanfordBunnyTestShrink.json";
const bool SAVE_TENSIONED_JSON = false;

const std::string OBJ_SAVE_FILE = "StanfordBunnyTestShrink.obj";

const bool SHOW_PLOT = false;


// bunny obj file (I transfered/reformatted it to a .pcd file)      35k points (so a fair bit considering it's in python and not paralyzed)
// https://github.com/alecjacobson/common-3d-test-models/blob/master/data/stanford-bunny.obj


// tensioning image series file information
const std::string TENSION_IMAGES_SAVE_PATH = "ShrinkWrap4/";
const bool SAVE_TENSION_IMAGES = false;

// loading a pcd file
const bool LOAD_PCD = false;
const std::string PCD_FILE = "standfordBunny.pcd";  // "test2.pcd";  // "TestFile.pcd"  # 
const double POINTS_LOADED_PERCENT = 100.0 / (100.0);

// test generators for certain shapes for testing
const bool GENERATE_HALLOW_SPHERE = true;  // fib sphere should be 30/68
const bool GENERATE_SOLID_SPHERE = false;
const bool GENERATE_CUBE = false;

// https://people.math.sc.edu/Burkardt/data/pcd/p213.pcd: TestFile.pcd

// the level at which the surface is defined as solid instead of void
double ISO_CONTOUR_LEVEL = 4.0;

// as long as the chunk size equals this any points that can make the surface solid will be within the neighboring 27 cells
double CHUNK_SIZE = ISO_CONTOUR_LEVEL;
const double SHELL_CHUNK_SIZE = 10;  // the chunk size for the points generated to represent the bubbly shell around the approximated objects
const int CHUNK_BUFFER_SIZE = 250;
const int POINT_CLOUD_SIZE = 50000;
const int SURFACE_POINTS_BUFFER_SIZE = 250000;

// the size and positioning of the calculated area
const int SAMPLING_SPACE_SIZE[3] = {101, 101, 101};  //[101, 101, 101]  # [500, 500, 500]#
double SAMPLING_SPACE_OFFSET[3] = {0, 0, 0};

// the maximum fill depth (to avoid any infinite loops in the case of an error in the code, part, or input settings)
const int MAX_FILL_DEPTH = 1000000000;  // 1 billion. Hopefully that's not too little, the grids can get fairly large

// the number of iterations of surface tension to smooth out the surface (the more the better the surface)
const bool SIMULATING_SURFACE_TENSION = false;
const int TENSION_ITTERATIONS = 2500;
const double VELOCITY = 0.75;  // constant velocity with respect to the gradient

// shifts all verticies inward at the end based on the surface normals of the scalar field (signed distance field)
const bool SCALE_VERTS_INWARD = false;
const double VERT_MOVEMENT = 0;

// the delta time and spacing for various components of the simulation (space and time)
double INVERSE_DELTA_X = 1 / (1);
double INVERSE_DELTA_Y = 1 / (1);  // these are the inverses to save computation speed sense it's almost always divided by
double INVERSE_DELTA_Z = 1 / (1);

const double DT = 0.0055 * 0.1;

// auto gets certain perameters like the grid size and offset (returns the results and sets them for the run)
const bool AUTO_SET = true;  // changes the offset and dt to fit the object into the specified grid size

// some random numbers to be saved
const double PI = 3.14159;



// utility functionsd
double Min (double v1, double v2)
{
    if (v1 < v2) return v1;
    return v2;
}

double Max (double v1, double v2)
{
    if (v1 > v2) return v1;
    return v2;
}



// stores and sorts a set of points into a grid of chunks
class ChunkGrid
{

private:

    // the grid for the chunks
    CArray <double> chunks;
    CArray <int> chunkStackedSize;
    // chunks = CArray <double> (x,y,z,w) ;

    // the bounds of the chunk
    int gridSizeX;
    int gridSizeY;
    int gridSizeZ;

    double smallBoundX;
    double smallBoundY;
    double smallBoundZ;
    double largeBoundX;
    double largeBoundY;
    double largeBoundZ;
    double chunkSize;
    int bufferSize;

    // initializing the grid
    public: ChunkGrid (int _gridSizeX, int _gridSizeY, int _gridSizeZ, int _smallBoundX, int _smallBoundY, int _smallBoundZ, int _largeBoundX, int _largeBoundY, int _largeBoundZ, double _chunkSize, int _bufferSize)
    {
        // creating the chunks
        chunks = CArray <double> (_gridSizeX, _gridSizeY, _gridSizeZ, _bufferSize, 3);
        chunkStackedSize = CArray <int> (_gridSizeX, _gridSizeY, _gridSizeZ);
        chunkStackedSize.set_values(0);
        bufferSize = _bufferSize;

        chunkSize = _chunkSize;

        // defining the bounds of the chunk grid
        gridSizeX = _gridSizeX;
        gridSizeY = _gridSizeY;
        gridSizeZ = _gridSizeZ;

        // defining the bounds of the data sets
        smallBoundX = _smallBoundX;
        smallBoundY = _smallBoundY;
        smallBoundZ = _smallBoundZ;
        largeBoundX = _largeBoundX;
        largeBoundY = _largeBoundY;
        largeBoundZ = _largeBoundZ;
    }

    // finds the nearest point to a given position
    public: double FindNearestPoint (double startingPositionX, double startingPositionY, double startingPositionZ)
    {
        // brining the point to the start of the grid to save time from searching empty space
        double correctedPositionX = Max(Min(startingPositionX, largeBoundX), smallBoundX);
        double correctedPositionY = Max(Min(startingPositionY, largeBoundY), smallBoundY);
        double correctedPositionZ = Max(Min(startingPositionZ, largeBoundZ), smallBoundZ);

        // finding the starting chunk position
        double startingChunkPositionX = std::floor((correctedPositionX - smallBoundX) / chunkSize);
        double startingChunkPositionY = std::floor((correctedPositionY - smallBoundY) / chunkSize);
        double startingChunkPositionZ = std::floor((correctedPositionZ - smallBoundZ) / chunkSize);

        // getting the distance to the nearest point to the position
        int searchRadius = 0;
        double finalDistance = ExpandingSearch(startingChunkPositionX, startingChunkPositionY, startingChunkPositionZ, startingPositionX, startingPositionY, startingPositionZ, searchRadius, 500);

        // checking if the point is further away than the nearest edge of the search (cubes aren't circular)
        if ((double) searchRadius < finalDistance)
        {
            // finding the max depth the search needs to be extended to
            int maxExtendedDepth = (int) std::ceil((double) searchRadius - finalDistance);

            // doing an extended search
            searchRadius++;  // searching beyond the current search
            double expandedDistance = ExpandingSearch(startingChunkPositionX, startingChunkPositionY, startingChunkPositionZ, startingPositionX, startingPositionY, startingPositionZ, searchRadius, maxExtendedDepth);

            // checking if the new search returned a closer point
            if (expandedDistance < finalDistance) finalDistance = expandedDistance;
        }

        // returning the distance found
        return finalDistance;
    }

    // recursively expands a search until chunks with points are found; then returns the nearest distance
    private: double ExpandingSearch (double chunkPositionX, double chunkPositionY, double chunkPositionZ, double positionX, double positionY, double positionZ, int &radius, int maxSearchDepth)
    {
        bool valid = false;  // for if no valid positions were found (not searching beyond bounds)
        double minDistance = 99999999.0;  // the minimum distance found

        // searching through all chunks in the location
        for (int x = -radius; x <= radius; x++)
        {
            // making sure the position is valid
            if (chunkPositionX + x >= 0 && chunkPositionX + x < gridSizeX)
            {
                

                // searching along the y axis
                for (int y = -radius; y <= radius; y++)
                {
                    // making sure the position is valid
                    if (chunkPositionY + y >= 0 && chunkPositionY + y < gridSizeY)
                    {
                        // checking each chunk for occupancy
                        if (chunkPositionZ - radius >= 0) {
                            valid = true;

                            // finding the number of points in the chunk
                            int numPointsInChunk = chunkStackedSize(chunkPositionX + x, chunkPositionY + y, chunkPositionZ - radius);
                            if (numPointsInChunk) {
                                // looping through all the points in the chunk
                                for (int i = 0; i < numPointsInChunk; i++)
                                {
                                    // checking the distance and checking if it's the new minimum
                                    double pointX = chunks(chunkPositionX + x, chunkPositionY + y, chunkPositionZ - radius, i, 0);
                                    double pointY = chunks(chunkPositionX + x, chunkPositionY + y, chunkPositionZ - radius, i, 1);
                                    double pointZ = chunks(chunkPositionX + x, chunkPositionY + y, chunkPositionZ - radius, i, 2);
                                    double distance = GetDistance(pointX, pointY, pointZ, positionX, positionY, positionZ);

                                    // checking if the point is the new closest
                                    if (distance < minDistance) minDistance = distance;
                                }
                            }
                        }
                        if (chunkPositionZ + radius < gridSizeZ) {
                            valid = true;

                            // finding the number of points in the chunk
                            int numPointsInChunk = chunkStackedSize(chunkPositionX + x, chunkPositionY + y, chunkPositionZ + radius);
                            if (numPointsInChunk) {
                                // looping through all the points in the chunk
                                for (int i = 0; i < numPointsInChunk; i++)
                                {
                                    // checking the distance and checking if it's the new minimum
                                    double pointX = chunks(chunkPositionX + x, chunkPositionY + y, chunkPositionZ + radius, i, 0);
                                    double pointY = chunks(chunkPositionX + x, chunkPositionY + y, chunkPositionZ + radius, i, 1);
                                    double pointZ = chunks(chunkPositionX + x, chunkPositionY + y, chunkPositionZ + radius, i, 2);
                                    double distance = GetDistance(pointX, pointY, pointZ, positionX, positionY, positionZ);

                                    // checking if the point is the new closest
                                    if (distance < minDistance) minDistance = distance;
                                }
                            }
                        }
                    }
                }


                // searching along the z axis
                for (int z = -radius; z <= radius; z++)
                {
                    // making sure the position is valid
                    if (chunkPositionZ + z >= 0 && chunkPositionZ + z < gridSizeZ)
                    {
                        // checking each chunk for occupancy
                        if (chunkPositionY - radius >= 0) {
                            valid = true;

                            // finding the number of points in the chunk
                            int numPointsInChunk = chunkStackedSize(chunkPositionX + x, chunkPositionY - radius, chunkPositionZ + z);
                            if (numPointsInChunk) {
                                // looping through all the points in the chunk
                                for (int i = 0; i < numPointsInChunk; i++)
                                {
                                    // checking the distance and checking if it's the new minimum
                                    double pointX = chunks(chunkPositionX + x, chunkPositionY - radius, chunkPositionZ + z, i, 0);
                                    double pointY = chunks(chunkPositionX + x, chunkPositionY - radius, chunkPositionZ + z, i, 1);
                                    double pointZ = chunks(chunkPositionX + x, chunkPositionY - radius, chunkPositionZ + z, i, 2);
                                    double distance = GetDistance(pointX, pointY, pointZ, positionX, positionY, positionZ);

                                    // checking if the point is the new closest
                                    if (distance < minDistance) minDistance = distance;
                                }
                            }
                        }
                        if (chunkPositionY + radius < gridSizeY) {
                            valid = true;

                            // finding the number of points in the chunk
                            int numPointsInChunk = chunkStackedSize(chunkPositionX + x, chunkPositionY + radius, chunkPositionZ + z);
                            if (numPointsInChunk) {
                                // looping through all the points in the chunk
                                for (int i = 0; i < numPointsInChunk; i++)
                                {
                                    // checking the distance and checking if it's the new minimum
                                    double pointX = chunks(chunkPositionX + x, chunkPositionY + radius, chunkPositionZ + z, i, 0);
                                    double pointY = chunks(chunkPositionX + x, chunkPositionY + radius, chunkPositionZ + z, i, 1);
                                    double pointZ = chunks(chunkPositionX + x, chunkPositionY + radius, chunkPositionZ + z, i, 2);
                                    double distance = GetDistance(pointX, pointY, pointZ, positionX, positionY, positionZ);

                                    // checking if the point is the new closest
                                    if (distance < minDistance) minDistance = distance;
                                }
                            }
                        }
                    }
                }
            }
        }

        // searching along the y and z axis
        for (int y = -radius; y <= radius; y++)
        {
            // making sure the position is valid
            if (chunkPositionY + y >= 0 && chunkPositionY + y < gridSizeY)
            {
                // searching through the z axis
                for (int z = -radius; z <= radius; z++)
                {
                    // making sure the position is valid
                    if (chunkPositionZ + z >= 0 && chunkPositionZ + z < gridSizeZ)
                    {
                        // checking each chunk for occupancy
                        if (chunkPositionX - radius >= 0) {
                            valid = true;

                            // finding the number of points in the chunk
                            int numPointsInChunk = chunkStackedSize(chunkPositionX - radius, chunkPositionY + y, chunkPositionZ + z);
                            if (numPointsInChunk) {
                                // looping through all the points in the chunk
                                for (int i = 0; i < numPointsInChunk; i++)
                                {
                                    // checking the distance and checking if it's the new minimum
                                    double pointX = chunks(chunkPositionX - radius, chunkPositionY + y, chunkPositionZ + z, i, 0);
                                    double pointY = chunks(chunkPositionX - radius, chunkPositionY + y, chunkPositionZ + z, i, 1);
                                    double pointZ = chunks(chunkPositionX - radius, chunkPositionY + y, chunkPositionZ + z, i, 2);
                                    double distance = GetDistance(pointX, pointY, pointZ, positionX, positionY, positionZ);

                                    // checking if the point is the new closest
                                    if (distance < minDistance) minDistance = distance;
                                }
                            }
                        }
                        if (chunkPositionX + radius < gridSizeZ) {
                            valid = true;

                            // finding the number of points in the chunk
                            int numPointsInChunk = chunkStackedSize(chunkPositionX + radius, chunkPositionY + y, chunkPositionZ + z);
                            if (numPointsInChunk) {
                                // looping through all the points in the chunk
                                for (int i = 0; i < numPointsInChunk; i++)
                                {
                                    // checking the distance and checking if it's the new minimum
                                    double pointX = chunks(chunkPositionX + radius, chunkPositionY + y, chunkPositionZ + z, i, 0);
                                    double pointY = chunks(chunkPositionX + radius, chunkPositionY + y, chunkPositionZ + z, i, 1);
                                    double pointZ = chunks(chunkPositionX + radius, chunkPositionY + y, chunkPositionZ + z, i, 2);
                                    double distance = GetDistance(pointX, pointY, pointZ, positionX, positionY, positionZ);

                                    // checking if the point is the new closest
                                    if (distance < minDistance) minDistance = distance;
                                }
                            }
                        }
                    }
                }
            }
        }


        // returning the distance if points have been found
        if (minDistance < 99999999.0 || !valid || radius >= maxSearchDepth) return minDistance;

        radius++;  // increasing the search radius
        return ExpandingSearch(chunkPositionX, chunkPositionY, chunkPositionZ, positionX, positionY, positionZ, radius, maxSearchDepth);
    }

    // gets the distance to a point
    private: double GetDistance (double pointX, double pointY, double pointZ, double positionX, double positionY, double positionZ)
    {
        // finding the difference between the two points
        double dx = positionX - pointX;
        double dy = positionY - pointY;
        double dz = positionZ - pointZ;
        
        // gettting the distance
        double distance = dx*dx + dy*dy + dz*dz;

        // returning the distance
        return distance;
    }

    // adds a new point to the chunk grid
    public: void AddChunk (int indexX, int indexY, int indexZ, double pointX, double pointY, double pointZ)
    {
        // getting the number of stacked points in the chunk
        int numberStacked = chunkStackedSize(indexX, indexY, indexZ) + 1;
        chunkStackedSize(indexX, indexY, indexZ) += 1;
        
        // adding the point
        chunks(indexX, indexY, indexZ, numberStacked, 0) = pointX;
        chunks(indexX, indexY, indexZ, numberStacked, 1) = pointY;
        chunks(indexX, indexY, indexZ, numberStacked, 2) = pointZ;
        
        return;
    }
};


// generates a sorted chunk grid from a list of points
ChunkGrid GenerateChunks (CArray <double> points, double chunkSize, int numPoints)
{
    // finding the bounds of the data set
    double smallBoundX = 999999999.0;
    double smallBoundY = 999999999.0;
    double smallBoundZ = 999999999.0;

    double largeBoundX = -999999999.0;
    double largeBoundY = -999999999.0;
    double largeBoundZ = -999999999.0;

    // going through the points and finding the boudns of the data set
    for (int i = 0; i < numPoints; i++)
    {
        // finding the minimum values
        smallBoundX = Min(points(i, 0), smallBoundX);
        smallBoundY = Min(points(i, 1), smallBoundX);
        smallBoundZ = Min(points(i, 2), smallBoundX);

        // finding the(int)  maximum values
        largeBoundX = Max(points(i, 0), largeBoundX);
        largeBoundY = Max(points(i, 1), largeBoundX);
        largeBoundZ = Max(points(i, 2), largeBoundX);
    }

    // bringing the bounds to whole numbers
    smallBoundX = std::floor(smallBoundX);
    smallBoundY = std::floor(smallBoundY);
    smallBoundZ = std::floor(smallBoundZ);

    largeBoundX = std::ceil(largeBoundX);
    largeBoundY = std::ceil(largeBoundY);
    largeBoundZ = std::ceil(largeBoundZ);

    // finding the grid size
    int gridSizeX = (int) std::ceil(abs(largeBoundX - smallBoundX) / chunkSize) + 1;
    int gridSizeY = (int) std::ceil(abs(largeBoundY - smallBoundY) / chunkSize) + 1;
    int gridSizeZ = (int) std::ceil(abs(largeBoundZ - smallBoundZ) / chunkSize) + 1;

    // creating the grid
    ChunkGrid chunkGrid = ChunkGrid(gridSizeX, gridSizeY, gridSizeZ, smallBoundX, smallBoundY, smallBoundZ, largeBoundX, largeBoundY, largeBoundZ, chunkSize, CHUNK_BUFFER_SIZE);

    // going through all the points and adding them to the chunk grid
    for (int i = 0; i < numPoints; i++)
    {
        // finding the chunk position
        int chunkPosX = (int) ((points(i, 0) - smallBoundX) / chunkSize);
        int chunkPosY = (int) ((points(i, 1) - smallBoundY) / chunkSize);
        int chunkPosZ = (int) ((points(i, 2) - smallBoundZ) / chunkSize);

        // adding the point
        chunkGrid.AddChunk(chunkPosX, chunkPosY, chunkPosZ, points(i, 0), points(i, 1), points(i, 2));
    }

    return chunkGrid;
}



// does a flood fill (used for calculating the signs)
void FloodFill (int indexX, int indexY, int indexZ, CArray <double> &distanceField, CArray <short> &signedGrid, short sign, int &numFilled)
{
    // all the directions neighbors can be in
    CArray <double> changeX = CArray <double> (6);
    changeX.set_values(0);
    changeX(0) = 1; changeX(1) = -1;
    CArray <double> changeY = CArray <double> (6);
    changeY.set_values(0);
    changeY(2) = 1; changeY(3) = -1;
    CArray <double> changeZ = CArray <double> (6);
    changeZ.set_values(0);
    changeZ(4) = 1; changeZ(5) = -1;

    // getting the grid size
    int gridSizeX = signedGrid.dims(0);
    int gridSizeY = signedGrid.dims(1);
    int gridSizeZ = signedGrid.dims(2);

    // sizeX(1 + sizeY(1 + sizeZ))
    // an array for neighboring points (max size: sizeX + sizeY*sizeX + sizeX*sizeY*sizeZ)
    int neighborsArraySize = gridSizeX * (1 + gridSizeY * (1 + gridSizeZ));
    CArray <int> neighborsX = CArray <int> (neighborsArraySize);
    CArray <int> neighborsY = CArray <int> (neighborsArraySize);
    CArray <int> neighborsZ = CArray <int> (neighborsArraySize);

    CArray <int> newNeighborsX = CArray <int> (neighborsArraySize);
    CArray <int> newNeighborsY = CArray <int> (neighborsArraySize);
    CArray <int> newNeighborsZ = CArray <int> (neighborsArraySize);

    CArray <int> newNeighborsIndexes = CArray <int> (CHUNK_BUFFER_SIZE);

    // adding the initial point
    //int indexKey = indexX + (indexY + indexZ * gridSizeY) * gridSizeX;
    neighborsX(0) = indexX;
    neighborsY(0) = indexY;
    neighborsZ(0) = indexZ;
    int numNeighbors = 1;

    // itterating till all points are flood filled
    for (int itteration = 0; itteration < MAX_FILL_DEPTH; itteration++)
    {
        int numNewNeighbors = 0;  // resetting the number of new neighbors

        // going through all current neighbors
        for (int neighborIndex = 0; neighborIndex < numNeighbors; neighborIndex++)
        {
            // checking if the current position is yet to be filled
            int xIndex = neighborsX(neighborIndex);
            int yIndex = neighborsY(neighborIndex);
            int zIndex = neighborsZ(neighborIndex);
            if (!signedGrid(xIndex, yIndex, zIndex))
            {
                // making sure the position is in a hollow space
                if (distanceField(xIndex, yIndex, zIndex) >= ISO_CONTOUR_LEVEL)
                {
                    // adding the point to the sign grid
                    signedGrid(xIndex, yIndex, zIndex) = sign;
                    numFilled++;

                    // adding all possible new neighbors
                    for (int i = 0; i < 6; i++)
                    {
                        // getting the neighbors position
                        int neighborXPos = xIndex + changeX(i);
                        int neighborYPos = yIndex + changeY(i);
                        int neighborZPos = zIndex + changeZ(i);

                        // checking if the point is valid
                        if (
                            neighborXPos >= 0 && neighborXPos < gridSizeX &&
                            neighborYPos >= 0 && neighborYPos < gridSizeY &&
                            neighborZPos >= 0 && neighborZPos < gridSizeZ
                            ) {

                                // adding the point to the new neighbors
                                int indexKey = xIndex + (yIndex + zIndex * gridSizeY) * gridSizeX;
                                newNeighborsX(indexKey) = xIndex;
                                newNeighborsY(indexKey) = yIndex;
                                newNeighborsZ(indexKey) = zIndex;
                                newNeighborsIndexes(numNewNeighbors) = indexKey;
                                numNewNeighbors++;

                        }
                    }
                }
            }
        }

        // checking if the search/fill is complete
        if (!numNewNeighbors) return;

        // copying the new neighbors to the old ones
        numNeighbors = numNewNeighbors;

        // looping through the new neighbors and adding them to the old ones
        for (int i = 0; i < numNewNeighbors; i++)
        {
            // getting the index for the point
            int indexKey = newNeighborsIndexes(i);  // doesn't need to be reset because a variable tracks its length

            // checking that the point isn't a duplicate
            if (newNeighborsX(indexKey) >= 0)
            {
                // clearing the positions from the array
                neighborsX(i) = newNeighborsX(indexKey);
                neighborsY(i) = newNeighborsY(indexKey);
                neighborsZ(i) = newNeighborsZ(indexKey);

                neighborsX(indexKey) = -1;  // the rest don't need to be cleared
            }
        }
    }

    return ;  // exiting the funciton
}


// calculates the internal/external sections of the object(s) which is represented by signs
void CalculateSigns (CArray <double> &distanceField, CArray <short> &signedGrid)
{
    int numFilled = 0;  // the number of tiles filled (can be used to find the percentage calculated more accuretly)

    // looping through every point on the grid to fill all regions
    for (int x = 0; x < SAMPLING_SPACE_SIZE[0]; x++)
    {
        // going over the y axis
        for (int y = 0; y < SAMPLING_SPACE_OFFSET[1]; y++)
        {
            short sign = 1;  // the entire perimeter (sides of the bounding box) should be void/air

            // looping through the z axis
            for (int z = 0; z < SAMPLING_SPACE_OFFSET[2]; z++)
            {
                // checking if the sign is already known
                if (signedGrid(x, y, z)) sign = signedGrid(x, y, z);
                else {
                    // checking if the point is outside of the wall of an object
                    if (distanceField(x, y, z) >= ISO_CONTOUR_LEVEL)
                    {
                        // flood filling the area starting at the current point
                        sign *= -1;
                        FloodFill(x, y, z, distanceField, signedGrid, sign, numFilled);
                    }
                }

            }
        }

    }

    return;  // exiting the function
}



// data loading functions
// add them here (don't feel like figuring it out right now)



// creates a fibonacci sphere
void FibonacciSphere (int samples, int &pointCloudStartingIndex, CArray <double> &pointCloud, double scalingFactor, double offset)
{
    double phi = PI * (sqrt(0.5) - 1.0);

    // going through all samples and generating the points
    for (int i = 0; i < samples; i++)
    {
        double y = 1 - (i / (float) (samples - 1)) * 2.0;  // y goes from 1 to -1
        double radius = sqrt(1 - y*y);  // radius at y

        double theta = phi * i;  // golden angle increment

        // getting the rest of the position
        double x = cos(theta) * radius;
        double z = sin(theta) * radius;

        // adding the point
        pointCloudStartingIndex++;
        pointCloud(pointCloudStartingIndex, 0) = x * scalingFactor + offset;
        pointCloud(pointCloudStartingIndex, 1) = y * scalingFactor + offset;
        pointCloud(pointCloudStartingIndex, 2) = z * scalingFactor + offset;
    }
}



// the main script
int main()
{

    // the distance field
    CArray <double> distanceField = CArray <double> (SAMPLING_SPACE_SIZE[0], SAMPLING_SPACE_SIZE[1], SAMPLING_SPACE_SIZE[2]);


    // loading or calculatiung the distance field
    if (LOADING_DISTANCE_FIELD_SAVE)
    {
        // load file
    }
    else
    {

        // the point cloud for all the initial points to be stored in
        CArray <double> pointCloud = CArray <double> (POINT_CLOUD_SIZE, 3);
        int pointCloudIndex = 0;

        // loading a pcd file
        if (LOAD_PCD)
        {
            // add pcd file loading (too lazy, I'll do it later)
        }
        else
        {

            // temporary shape generators for testing
            if (GENERATE_HALLOW_SPHERE) {
                FibonacciSphere(575, pointCloudIndex, pointCloud, 34, 50);
                FibonacciSphere(250, pointCloudIndex, pointCloud, 15, 50);
            }
            if (GENERATE_SOLID_SPHERE) {
                FibonacciSphere(200, pointCloudIndex, pointCloud, 14, 50);
            }
            if (GENERATE_CUBE) {
                for (int x = 0; x < 8; x++)
                {
                    for (int y = 0; y < 8; y++)
                    {
                        for (int z = 0; z < 8; z++)
                        {
                            if( x==7 || y==7 || z==7 || !x*y*z)
                            {
                                // checkinh if the point is already in the point cloud
                                bool valid = true;
                                for (int i = 0; i < pointCloudIndex; i++)
                                {
                                    if (pointCloud(i, 0) == x && pointCloud(i, 1) == y && pointCloud(i, 2) == z) {valid = false; break;}
                                }

                                // adding the point
                                if (valid)
                                {
                                    pointCloud(pointCloudIndex, 0) = x;
                                    pointCloud(pointCloudIndex, 1) = y;
                                    pointCloud(pointCloudIndex, 2) = z;
                                    pointCloudIndex++;
                                }
                            }
                        }
                    }
                }
            }


            // auto setting some of the parameters
            if (AUTO_SET)
            {
                // to store the bounds of the point cloud
                int lowestX = 999999999;
                int lowestY = 999999999;
                int lowestZ = 999999999;

                int largestX = -999999999;
                int largestY = -999999999;
                int largestZ = -999999999;

                // going through the point cloud and finding the bounds
                for (int i = 0; i < pointCloudIndex; i++)
                {
                    // updating the lowest
                    lowestX = Min(lowestX, pointCloud(i, 0));
                    lowestY = Min(lowestY, pointCloud(i, 1));
                    lowestZ = Min(lowestZ, pointCloud(i, 2));

                    // updating the largest
                    largestX = Max(largestX, pointCloud(i, 0));
                    largestY = Max(largestY, pointCloud(i, 1));
                    largestZ = Max(largestZ, pointCloud(i, 2));
                }

                // using the bounds to adjust the parameters to fit the object
                double cloudSizeX = largestX - lowestX;
                double cloudSizeY = largestY - lowestY;
                double cloudSizeZ = largestZ - lowestZ;

                // adding space around the object to the point cloud size and calculating other parameters
                double newDeltaScaleX = cloudSizeX / SAMPLING_SPACE_SIZE[0];
                double newDeltaScaleY = cloudSizeY / SAMPLING_SPACE_SIZE[1];
                double newDeltaScaleZ = cloudSizeZ / SAMPLING_SPACE_SIZE[2];

                double newIsoContourLevel = (newDeltaScaleX + newDeltaScaleY + newDeltaScaleZ) / 3 * ISO_CONTOUR_LEVEL;
                
                // updating the cloud size
                cloudSizeX += newIsoContourLevel * 12;
                cloudSizeY += newIsoContourLevel * 12;
                cloudSizeZ += newIsoContourLevel * 12;

                // updating the grid offset and delta position scaling
                double newGridOffsetX = lowestX - newIsoContourLevel * 6;
                double newGridOffsetY = lowestY - newIsoContourLevel * 6;
                double newGridOffsetZ = lowestZ - newIsoContourLevel * 6;

                newDeltaScaleX = cloudSizeX / SAMPLING_SPACE_SIZE[0];
                newDeltaScaleY = cloudSizeY / SAMPLING_SPACE_SIZE[1];
                newDeltaScaleZ = cloudSizeZ / SAMPLING_SPACE_SIZE[2];

                // calculating the new inverse delta scaling factors
                INVERSE_DELTA_X = 1 / newDeltaScaleX;
                INVERSE_DELTA_Y = 1 / newDeltaScaleY;
                INVERSE_DELTA_Z = 1 / newDeltaScaleZ;

                // updating the Iso contour level and chunk size
                double dAvg = (newDeltaScaleX + newDeltaScaleY + newDeltaScaleZ) / 3;
                ISO_CONTOUR_LEVEL *= dAvg;

                CHUNK_SIZE *= dAvg;

                // updating the offset
                SAMPLING_SPACE_OFFSET[0] = newGridOffsetX;
                SAMPLING_SPACE_OFFSET[1] = newGridOffsetY;
                SAMPLING_SPACE_OFFSET[2] = newGridOffsetZ;
            }


            // generating a new distance function
            
            // putting all the points of the point cloud into chunks
            ChunkGrid chunkGrid = GenerateChunks(pointCloud, CHUNK_SIZE, pointCloudIndex);

            std::cout << "chunked point cloud" << std::endl;

            // looping through all points and finding the point that is nearest to the cell
            for (int x_ = 0; x_ < SAMPLING_SPACE_SIZE[0]; x_++)
            {
                for (int y_ = 0; y_ < SAMPLING_SPACE_SIZE[1]; y_++)
                {
                    for (int z_ = 0; z_ < SAMPLING_SPACE_SIZE[2]; z_++)
                    {
                        // getting the global position of the point
                        double x = (double) x_ / INVERSE_DELTA_X + SAMPLING_SPACE_OFFSET[0];
                        double y = (double) y_ / INVERSE_DELTA_Y + SAMPLING_SPACE_OFFSET[1];
                        double z = (double) z_ / INVERSE_DELTA_Z + SAMPLING_SPACE_OFFSET[2];

                        // getting the nearest position to the point
                        double distance = chunkGrid.FindNearestPoint(x, y, z);
                        distanceField(x_, y_, z_) = distance;
                    }
                }
            }

            std::cout << "found all distances" << std::endl;

            // save the distance field (too lazy, add later)
        }


        // correct old sdf's to work with the generation of the new ones
        if (USING_OLD_SDF)
        {
            // looping through all the points and taking the absolute value of it
            for (int x = 0; x < SAMPLING_SPACE_SIZE[0]; x++)
            {
                for (int y = 0; y < SAMPLING_SPACE_SIZE[1]; y++)
                {
                    for (int z = 0; z < SAMPLING_SPACE_SIZE[2]; z++)
                    {
                        distanceField(x, y, z) = abs(distanceField(x, y, z));
                    }
                }
            }
        }


        // generating an updated sdf
        if (GENERATE_NEW_VERSION_SDF)
        {
            // generating a shell around all objects

            // a grid for the signs
            CArray <short> signedGrid = CArray <short> (SAMPLING_SPACE_SIZE[0], SAMPLING_SPACE_SIZE[1], SAMPLING_SPACE_SIZE[2]);

            std::cout << "ready to calculate signs" << std::endl;

            // calculating the signs for the objects and distance field
            CalculateSigns(distanceField, signedGrid);

            std::cout << "signs calculated" << std::endl;

            // an array for all surface points on the grid
            CArray <double> surfacePoints = CArray <double> (SURFACE_POINTS_BUFFER_SIZE, 3);
            int surfaceIndex = 0;  // the current index for the surface points array

            // itterating over every point and calculating if it's a surface point or not
            for (int x = 1; x < SAMPLING_SPACE_SIZE[0] - 1; x++)
            {
                for (int y = 1; y < SAMPLING_SPACE_SIZE[1] - 1; y++)
                {
                    for (int z = 1; z < SAMPLING_SPACE_SIZE[2] - 1; z++)
                    {
                        // making sure there's a hollow point next to the cell
                        if (
                            signedGrid(x - 1, y, z) > 0 && signedGrid(x + 1, y, z) > 0 &&
                            signedGrid(x, y - 1, z) > 0 && signedGrid(x, y + 1, z) > 0 &&
                            signedGrid(x, y, z - 1) > 0 && signedGrid(x, y, z + 1) > 0
                            ) {
                            
                            // checking if the current position is inside an object / on the inside of the surface wall of an object
                            if (signedGrid(x, y, z) < 0 || distanceField(x, y, z) < ISO_CONTOUR_LEVEL)
                            {
                                // adding the surface point to the array of surface points
                                surfacePoints(surfaceIndex, 0) = x;
                                surfacePoints(surfaceIndex, 1) = y;
                                surfacePoints(surfaceIndex, 2) = z;
                                surfaceIndex++;
                            }

                        }
                    }
                }
            }

            std::cout << "surface points gathered" << std::endl;

            // sorting the points into chunks
            ChunkGrid chunkedSurfacePoints = GenerateChunks(surfacePoints, SHELL_CHUNK_SIZE, surfaceIndex);
            
            std::cout << "surface points chunked" << std::endl;
            
            // correcting the iso contour sense the grid spacing has changed for this next step
            double oldIsoContourLevel;  // how is this broken?????
            std::cout << "old sio defined" << std::endl;
            if (AUTO_SET && !LOADING_DISTANCE_FIELD_SAVE) {
                std::cout << "into if" << std::endl;
                oldIsoContourLevel = ISO_CONTOUR_LEVEL;
                std::cout << "old iso set" << std::endl;
                ISO_CONTOUR_LEVEL = ISO_CONTOUR_LEVEL * 3 / (1/INVERSE_DELTA_X + 1/INVERSE_DELTA_Y + 1/INVERSE_DELTA_Z);
                std::cout << "iso set" << std::endl;
            }
            std::cout << "auto calculated" << std::endl;

            // finding all the new distances
            for (int x = 0; x < SAMPLING_SPACE_SIZE[0]; x++)
            {
                for (int y = 0; y < SAMPLING_SPACE_SIZE[1]; y++)
                {
                    for (int z = 0; z < SAMPLING_SPACE_SIZE[2]; z++)
                    {
                        std::cout << "before distance check" << std::endl;
                        // finding the distance to the nearest point
                        double distance = chunkedSurfacePoints.FindNearestPoint(x, y, z);

                        std::cout << "found nearest point" << std::endl;
                        // finding the correct sign for the current position
                        short sign = signedGrid(x, y, z);
                        if (distanceField(x, y, z) < oldIsoContourLevel && sign > -1) sign = -1;
                        if (!sign) sign = 1;

                        std::cout << "found sign" << std::endl;

                        // updating the signed distance funciton with the distance and sign for the given position
                        distanceField(x, y, z) = sqrt(distance) * (double) sign + ISO_CONTOUR_LEVEL;
                        std::cout << "got signed distance" << std::endl;
                    }
                }
            }

            std::cout << "surface distance found" << std::endl;

        }


        // simulating surface tension
        if (SIMULATING_SURFACE_TENSION)
        {
            // adding simulation here (too lazy rn, add later)
        }


        // updating the signs for older-type sdf's
        if (!GENERATE_NEW_VERSION_SDF && USING_OLD_SDF)
        {
            // a grid for the signs
            CArray <short> signedGrid = CArray <short> (SAMPLING_SPACE_SIZE[0], SAMPLING_SPACE_SIZE[1], SAMPLING_SPACE_SIZE[2]);

            // calculating the signs for the objects and distance field
            CalculateSigns(distanceField, signedGrid);

            // creating a distance field with correct signs
            CArray <double> signedDistanceField = CArray <double> (SAMPLING_SPACE_SIZE[0], SAMPLING_SPACE_SIZE[1], SAMPLING_SPACE_SIZE[2]);

            // itterating through all points and correcting the signs
            for (int x = 0; x < SAMPLING_SPACE_SIZE[0]; x++)
            {
                for (int y = 0; y < SAMPLING_SPACE_SIZE[1]; y++)
                {
                    for (int z = 0; z < SAMPLING_SPACE_SIZE[2]; z++)
                    {
                        signedDistanceField(x, y, z) = distanceField(x, y, z) * signedGrid(x, y, z);
                    }
                }
            }
        }


        // other stuff after this

    }
    
    return 0;
}


