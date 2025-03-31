#include <iostream>
#include <fstream>
#include <string>
#include <ranges>
#include <vector>
#include <sstream>

#include <chrono>
#include "matar.h"

#include "marchingCubes.h"
#include "octree.hpp"

using namespace mtr;  // so I don't have to write mtr:: a thousand times per line


//=========================================================================================================
//                                          Program Settings
//=========================================================================================================


const bool LOADING_DISTANCE_FIELD_SAVE = false;
const bool GENERATE_NEW_VERSION_SDF = true;
const bool USING_OLD_SDF = true;

//const std::string FIELD_FILE = "tensioned_StandFordBunny_untensioned_workingFixed!.json";  // "tensioned_StandFordBunny_untensioned_test2.json"
//const std::string SAVE_FILE = "StanfordBunnyTestShrink.json";
//const bool SAVE_TENSIONED_JSON = false;

//const std::string OBJ_SAVE_FILE = "StanfordBunnyTestShrink.obj";

//const bool SHOW_PLOT = false;


// bunny obj file (I transfered/reformatted it to a .pcd file)      35k points (so a fair bit considering it's in python and not paralyzed)
// https://github.com/alecjacobson/common-3d-test-models/blob/master/data/stanford-bunny.obj


// tensioning image series file information
const std::string TENSION_IMAGES_SAVE_PATH = "ShrinkWrap4/";
const bool SAVE_TENSION_IMAGES = false;

// loading a pcd file
const bool LOAD_PCD = true;
const std::string PCD_FILE = "dragon.pcd.ply";//"standfordBunny.pcd";  // "test2.pcd";  // "TestFile.pcd"  # 
const double POINTS_LOADED_PERCENT = 100.0 / (100.0);

// test generators for certain shapes for testing
const bool GENERATE_HALLOW_SPHERE = false;  // fib sphere should be 30/68
const bool GENERATE_SOLID_SPHERE = false;
const bool GENERATE_CUBE = false;

// https://people.math.sc.edu/Burkardt/data/pcd/p213.pcd: TestFile.pcd

// the level at which the surface is defined as solid instead of void
double ISO_CONTOUR_LEVEL = 6.0;

// octree setup settings
const int OCTREE_POINT_BUFFER_SIZE = 25;
const int MAX_OCTREE_DEPTH = 10;

// as long as the chunk size equals this any points that can make the surface solid will be within the neighboring 27 cells
double CHUNK_SIZE = ISO_CONTOUR_LEVEL;
const double SHELL_CHUNK_SIZE = 10;  // the chunk size for the points generated to represent the bubbly shell around the approximated objects
//const int CHUNK_BUFFER_SIZE = 25;
const int POINT_CLOUD_SIZE = 50000;
const int SURFACE_POINTS_BUFFER_SIZE = 250000;

// the size and positioning of the calculated area
const int SAMPLING_SPACE_SIZE[3] = {250, 250, 250};  //[101, 101, 101]  # [500, 500, 500]#
double SAMPLING_SPACE_OFFSET[3] = {0, 0, 0};

// the maximum fill depth (to avoid any infinite loops in the case of an error in the code, part, or input settings)
const int MAX_FILL_DEPTH = 100000000;  // 100 million. Hopefully that's not too little, the grids can get fairly large

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



//=========================================================================================================
//                                        Basic Utility Functions
//=========================================================================================================


// utility functions
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



//=========================================================================================================
//                                       Distance Field Generation
//=========================================================================================================


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
        chunks = CArray <double> (_gridSizeX, _gridSizeY, _gridSizeZ, _bufferSize, 3);  // failing right here for memory
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
        int startingChunkPositionX = (int) ((correctedPositionX - smallBoundX) / chunkSize);
        int startingChunkPositionY = (int) ((correctedPositionY - smallBoundY) / chunkSize);
        int startingChunkPositionZ = (int) ((correctedPositionZ - smallBoundZ) / chunkSize);

        // getting the distance to the nearest point to the position
        int searchRadius = 0;
        double finalDistance = sqrt(ExpandingSearch(startingChunkPositionX, startingChunkPositionY, startingChunkPositionZ, startingPositionX, startingPositionY, startingPositionZ, searchRadius, 500));
        
        // checking if the point is further away than the nearest edge of the search (cubes aren't circular)
        if ((double) searchRadius < finalDistance)
        {
            // finding the max depth the search needs to be extended to
            int maxExtendedDepth = (int) std::ceil((double) searchRadius - finalDistance);

            // doing an extended search
            searchRadius++;  // searching beyond the current search
            double expandedDistance = sqrt(ExpandingSearch(startingChunkPositionX, startingChunkPositionY, startingChunkPositionZ, startingPositionX, startingPositionY, startingPositionZ, searchRadius, maxExtendedDepth));

            // checking if the new search returned a closer point
            if (expandedDistance < finalDistance) finalDistance = expandedDistance;
        }

        // returning the distance found
        return finalDistance;
    }

    // recursively expands a search until chunks with points are found; then returns the nearest distance
    private: double ExpandingSearch (int chunkPositionX, int chunkPositionY, int chunkPositionZ, double positionX, double positionY, double positionZ, int &radius, int maxSearchDepth)
    {
        bool valid = false;  // for if no valid positions were found (not searching beyond bounds)
        double minDistance = 99999999.0;  // the minimum distance found

        // searching through all chunks in the location
        for (int x = -radius; x <= radius; x++)
        {
            // making sure the position is valid
            if (chunkPositionX + x < 0 || chunkPositionX + x >= gridSizeX) continue;
            
            // searching along the y axis
            for (int y = -radius + 1; y <= radius; y++)
            {
                // making sure the position is valid
                if (chunkPositionY + y < 0 || chunkPositionY + y >= gridSizeY) continue;
                
                // checking each chunk for occupancy
                if (chunkPositionZ - radius >= 0) {
                    valid = true;
                    
                    // finding the number of points in the chunk
                    int numPointsInChunk = chunkStackedSize(chunkPositionX + x, chunkPositionY + y, chunkPositionZ - radius);
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
                if (chunkPositionZ + radius < gridSizeZ) {
                    valid = true;

                    // finding the number of points in the chunk
                    int numPointsInChunk = chunkStackedSize(chunkPositionX + x, chunkPositionY + y, chunkPositionZ + radius);
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


            // searching along the z axis
            for (int z = -radius; z <= radius; z++)
            {
                // making sure the position is valid
                if (chunkPositionZ + z < 0 || chunkPositionZ + z >= gridSizeZ) continue;

                // checking each chunk for occupancy
                if (chunkPositionY - radius >= 0) {
                    valid = true;

                    // finding the number of points in the chunk
                    int numPointsInChunk = chunkStackedSize(chunkPositionX + x, chunkPositionY - radius, chunkPositionZ + z);
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
                if (chunkPositionY + radius < gridSizeY) {
                    valid = true;

                    // finding the number of points in the chunk
                    int numPointsInChunk = chunkStackedSize(chunkPositionX + x, chunkPositionY + radius, chunkPositionZ + z);
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


        // searching along the y and z axis
        for (int y = -radius + 1; y <= radius; y++)
        {
            // making sure the position is valid
            if (chunkPositionY + y < 0 || chunkPositionY + y >= gridSizeY) continue;

            // searching through the z axis
            for (int z = -radius + 1; z <= radius; z++)
            {
                // making sure the position is valid
                if (chunkPositionZ + z < 0 || chunkPositionZ + z >= gridSizeZ) continue;

                // checking each chunk for occupancy
                if (chunkPositionX - radius >= 0) {
                    valid = true;
                    
                    // finding the number of points in the chunk
                    int numPointsInChunk = chunkStackedSize(chunkPositionX - radius, chunkPositionY + y, chunkPositionZ + z);
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
                if (chunkPositionX + radius < gridSizeX) {
                    valid = true;

                    // finding the number of points in the chunk
                    int numPointsInChunk = chunkStackedSize(chunkPositionX + radius, chunkPositionY + y, chunkPositionZ + z);
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


        // returning the distance if points have been found
        if (minDistance < 99999999.0 || radius >= maxSearchDepth) return minDistance;
        if (!valid) return 0;  // no valid positions are available (every chunk has been searched)       shouldn't really ever happen

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
        int numberStacked = chunkStackedSize(indexX, indexY, indexZ);
        chunkStackedSize(indexX, indexY, indexZ) += 1;
        
        // adding the point
        chunks(indexX, indexY, indexZ, numberStacked, 0) = pointX;
        chunks(indexX, indexY, indexZ, numberStacked, 1) = pointY;
        chunks(indexX, indexY, indexZ, numberStacked, 2) = pointZ;
        
        return;
    }
};


// generates a sorted chunk grid from a list of points
ChunkGrid GenerateChunks (CArray <double> &points, double chunkSize, int numPoints)
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
    int gridSizeX = (int) std::ceil(abs(largeBoundX - smallBoundX) / chunkSize) + 9;
    int gridSizeY = (int) std::ceil(abs(largeBoundY - smallBoundY) / chunkSize) + 9;
    int gridSizeZ = (int) std::ceil(abs(largeBoundZ - smallBoundZ) / chunkSize) + 9;

    // creating the grid
    int bufferSize = gridSizeX * gridSizeY * gridSizeZ;  // the maximum number of possible points per chunk
    ChunkGrid chunkGrid = ChunkGrid(gridSizeX, gridSizeY, gridSizeZ, smallBoundX, smallBoundY, smallBoundZ, largeBoundX, largeBoundY, largeBoundZ, chunkSize, bufferSize);

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



//=========================================================================================================
//                                          Sign Calculations
//=========================================================================================================


// does a flood fill (used for calculating the signs)
void FloodFill (int startX, int startY, int startZ, CArray <double> &distanceField, CArray <int8_t> &signedGrid, int8_t sign, int &numFilled)
{
    // all the directions neighbors can be in
    CArray <int> changeX = CArray <int> (6);
    changeX.set_values(0);
    changeX(0) = 1; changeX(1) = -1;
    CArray <int> changeY = CArray <int> (6);
    changeY.set_values(0);
    changeY(2) = 1; changeY(3) = -1;
    CArray <int> changeZ = CArray <int> (6);
    changeZ.set_values(0);
    changeZ(4) = 1; changeZ(5) = -1;

    // getting the grid size
    int gridSizeX = signedGrid.dims(0);
    int gridSizeY = signedGrid.dims(1);
    int gridSizeZ = signedGrid.dims(2);

    // sizeX(1 + sizeY(1 + sizeZ))
    // an array for neighboring points (max size: sizeX + sizeY*sizeX + sizeX*sizeY*sizeZ)
    int neighborsArraySize = gridSizeX * (1 + gridSizeY * (1 + gridSizeZ));     // change these to the addresses of them defined in the calculate signs function to save time allocating memory repeatedly
    CArray <int> neighborsX = CArray <int> (neighborsArraySize);
    CArray <int> neighborsY = CArray <int> (neighborsArraySize);
    CArray <int> neighborsZ = CArray <int> (neighborsArraySize);
    
    CArray <int> newNeighborsX = CArray <int> (neighborsArraySize);
    CArray <int> newNeighborsY = CArray <int> (neighborsArraySize);
    CArray <int> newNeighborsZ = CArray <int> (neighborsArraySize);  // move this allocation earlier on

    CArray <int> newNeighborsIndexes = CArray <int> (neighborsArraySize);

    // adding the initial point
    //int indexKey = indexX + (indexY + indexZ * gridSizeY) * gridSizeX;
    neighborsX(0) = startX;
    neighborsY(0) = startY;
    neighborsZ(0) = startZ;
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
                        //std::cout << "pre change" << std::endl;
                        int neighborXPos = xIndex + changeX(i);
                        int neighborYPos = yIndex + changeY(i);
                        int neighborZPos = zIndex + changeZ(i);
                        //std::cout << "post change" << std::endl;

                        // checking if the point is valid
                        if (
                            neighborXPos >= 0 && neighborXPos < gridSizeX &&
                            neighborYPos >= 0 && neighborYPos < gridSizeY &&
                            neighborZPos >= 0 && neighborZPos < gridSizeZ
                            ) {

                                // adding the point to the new neighbors
                                int indexKey = neighborXPos + (neighborYPos + neighborZPos * gridSizeY) * gridSizeX;
                                //std::cout << "ready to add points; " << neighborsArraySize << ";      " << neighborXPos << ", " << neighborYPos << ", " << neighborZPos << "        :" << indexKey << std::endl;
                                newNeighborsX(indexKey) = neighborXPos;
                                newNeighborsY(indexKey) = neighborYPos;
                                newNeighborsZ(indexKey) = neighborZPos;
                                //std::cout << "final adding: " << numNewNeighbors << "        : " <<  << std::endl;
                                newNeighborsIndexes(numNewNeighbors) = indexKey;
                                numNewNeighbors++;
                                //std::cout << "points added" << std::endl;

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

                newNeighborsX(indexKey) = -1;  // the rest don't need to be cleared
            }
        }
    }

    return ;  // exiting the funciton
}


// calculates the internal/external sections of the object(s) which is represented by signs
void CalculateSigns (CArray <double> &distanceField, CArray <int8_t> &signedGrid)
{
    int numFilled = 0;  // the number of tiles filled (can be used to find the percentage calculated more accuretly)

    // filling the initial space around the objects
    FloodFill(0, 0, 0, distanceField, signedGrid, 1, numFilled);  // this isn't actually filling anything. This needs fixing
    std::cout << "Sign at 0, 0, 0: " << (int) signedGrid(1, 1, 1) << std::endl;
    std::cout << "Distance at 0, 0, 0: " << distanceField(1, 1, 1) << std::endl;

    // looping through every point on the grid to fill all regions
    for (int x = 0; x < SAMPLING_SPACE_SIZE[0]; x++)
    {
        // going over the y axis
        for (int y = 0; y < SAMPLING_SPACE_SIZE[1]; y++)
        {
            int8_t sign = 1;  // the entire perimeter (sides of the bounding box) should be void/air

            // looping through the z axis
            for (int z = 0; z < SAMPLING_SPACE_SIZE[2]; z++)
            {
                // checking if the sign is already known
                if (signedGrid(x, y, z)) sign = signedGrid(x, y, z);
                else {
                    // checking if the point is outside of the wall of an object
                    if (distanceField(x, y, z) >= ISO_CONTOUR_LEVEL)
                    {
                        // flood filling the area starting at the current point
                        sign *= -1;
                        
                        std::cout << "(before) filled: " << numFilled << "      sign is: " << (int) sign << "\nPosition: " << x << ", " << y << ", " << z << std::endl;
                        FloodFill(x, y, z, distanceField, signedGrid, sign, numFilled);
                        std::cout << "(after) filled: " << numFilled << std::endl;
                    }
                }
            }
        }
    }

    return;  // exiting the function
}



//=========================================================================================================
//                                       Data Loading & Saving
//=========================================================================================================


// data loading functions
// add them here (don't feel like figuring it out right now)




//=========================================================================================================
//                                          Shape Generation
//=========================================================================================================


// creates a fibonacci sphere
void FibonacciSphere (int samples, int &pointCloudStartingIndex, CArray <double> &pointCloud, double scalingFactor, double offset)
{
    double phi = PI * (sqrt(0.5) - 1.0);

    // going through all samples and generating the points
    for (int i = 0; i < samples; i++)
    {
        double y = 1 - (i / (double) (samples - 1)) * 2.0;  // y goes from 1 to -1
        double radius = sqrt(1 - y*y);  // radius at y

        double theta = phi * (double) i;  // golden angle increment

        // getting the rest of the position
        double x = cos(theta) * radius;
        double z = sin(theta) * radius;

        // adding the point
        pointCloud(pointCloudStartingIndex, 0) = x * scalingFactor + offset;
        pointCloud(pointCloudStartingIndex, 1) = y * scalingFactor + offset;
        pointCloud(pointCloudStartingIndex, 2) = z * scalingFactor + offset;
        pointCloudStartingIndex++;
    }
}



//=========================================================================================================
//                                          Main Script
//=========================================================================================================


// the main script
int main()
{

    // the distance field
    CArray <double> distanceField = CArray <double> (SAMPLING_SPACE_SIZE[0], SAMPLING_SPACE_SIZE[1], SAMPLING_SPACE_SIZE[2]);


    // loading or calculatiung the distance field
    if (LOADING_DISTANCE_FIELD_SAVE)
    {
        //
    }
    else
    {

        // the point cloud for all the initial points to be stored in
        CArray <double> pointCloud = CArray <double> (POINT_CLOUD_SIZE, 3);
        int pointCloudIndex = 0;

        // loading a pcd file
        if (LOAD_PCD)
        {
            std::cout << "Read to load pcd" << std::endl;
            // add pcd file loading (too lazy, I'll do it later)
            // load file
            /*
            0:  # .PCD v.5 - Point Cloud Data file format
            1:  VERSION .5
            2:  FIELDS x y z
            3:  SIZE 4 4 4
            4:  TYPE F F F
            5:  COUNT 1 1 1
            6:  WIDTH 35947
            7:  HEIGHT 1
            8:  POINTS 35947
            9:  DATA ascii
            10: -0.037830 0.127940 0.004475
            */
            std::string line;
            std::ifstream myfile (PCD_FILE);
            if (myfile.is_open())
            {
                //std::cout << "Valid file" << std::endl;
                int index = 0;
                while ( getline (myfile, line, '\n') )
                {
                    //std::cout << "Loading a line" << std::endl;
                    if (index > 9) {
                        int xyzIndex = 0;
                        std::string numStr;
                        std::stringstream lineStream (line);
                        while(getline(lineStream, numStr, ' '))
                        {
                            //std::cout << "loading number: " << numStr << ": (" << index << ", " << xyzIndex << ") -> " << std::stod(numStr) << std::endl;
                            pointCloud(pointCloudIndex, xyzIndex) = std::stod(numStr) * 250.;
                            xyzIndex += 1;
                        }
                        //std::cout << "Number loaded: " << pointCloud(pointCloudIndex, 0) << ", " <<  pointCloud(pointCloudIndex, 1) << ", " <<  pointCloud(pointCloudIndex, 2) << std::endl;
                        pointCloudIndex ++;
                    }
                    index += 1;
                }
                myfile.close();
            } else {
                std::cout << "NOOOOOO" << std::endl;
            }
            std::cout << "Num points after loading: " << pointCloudIndex << std::endl;
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
        }

        std::cout << "Num points before auto: " << pointCloudIndex << std::endl;

        // auto setting some of the parameters
        double samplingSpaceX = SAMPLING_SPACE_SIZE[0];
        double samplingSpaceY = SAMPLING_SPACE_SIZE[1];
        double samplingSpaceZ = SAMPLING_SPACE_SIZE[2];
        if (AUTO_SET)
        {
            // to store the bounds of the point cloud
            double lowestX = 999999999.0;
            double lowestY = 999999999.0;
            double lowestZ = 999999999.0;

            double largestX = -999999999.0;
            double largestY = -999999999.0;
            double largestZ = -999999999.0;

            // going through the point cloud and finding the bounds
            for (int i = 0; i < pointCloudIndex; i++)
            {
                // updating the lowest
                lowestX = std::min <double> (lowestX, pointCloud(i, 0));
                lowestY = std::min <double> (lowestY, pointCloud(i, 1));
                lowestZ = std::min <double> (lowestZ, pointCloud(i, 2));

                // updating the largest
                largestX = std::max <double> (largestX, pointCloud(i, 0));
                largestY = std::max <double> (largestY, pointCloud(i, 1));
                largestZ = std::max <double> (largestZ, pointCloud(i, 2));
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
            cloudSizeX += newIsoContourLevel * 4;
            cloudSizeY += newIsoContourLevel * 4;
            cloudSizeZ += newIsoContourLevel * 4;

            // updating the grid offset and delta position scaling
            double newGridOffsetX = lowestX - newIsoContourLevel * 2;
            double newGridOffsetY = lowestY - newIsoContourLevel * 2;
            double newGridOffsetZ = lowestZ - newIsoContourLevel * 2;

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

            // updating the bounding box size information
            samplingSpaceX = cloudSizeX;
            samplingSpaceY = cloudSizeY;
            samplingSpaceZ = cloudSizeZ;

            std::cout << "sampling offset: " << SAMPLING_SPACE_OFFSET[0] << ", " << SAMPLING_SPACE_OFFSET[1] << ", " << SAMPLING_SPACE_OFFSET[2] << "\nIso: " << ISO_CONTOUR_LEVEL << "     deltas: " << 1/INVERSE_DELTA_X << ", " << 1/INVERSE_DELTA_Y << ", " << 1/INVERSE_DELTA_Z << "\nCloud size: " << cloudSizeX << ", " << cloudSizeY << ", " << cloudSizeZ << std::endl;
            std::cout << "lowest: " << lowestX << ", " << lowestY << ", " << lowestZ << "\nHighest: " << largestX << ", " << largestY << ", " << largestZ << std::endl;
            std::cout << "Chunk size: " << CHUNK_SIZE << std::endl;



            //Octree octree = Octree(pointCloud, pointCloudIndex, cloudSizeX, cloudSizeY, cloudSizeZ, SAMPLING_SPACE_OFFSET[0], SAMPLING_SPACE_OFFSET[1], SAMPLING_SPACE_OFFSET[2], 10, 25);
            //octree.SubDivide();  // somehow this is all working


            /*    test code
            double testX = pointCloud(7, 0);
            double testY = pointCloud(7, 1);
            double testZ = pointCloud(7, 2);
            std::cout << "test pos: " << testX << ", " << testY << ", " << testZ << std::endl;
            int node = octree.GetLeafIndex(testX, testY, testZ);  // not sure if it's finding the points alright or not?
            std::cout << "node: " << node << std::endl;
            int pointIndex = octree.positionIndexesPlusOne(node, 0);
            std::cout << "point index: " << pointIndex << "        total points: " << pointCloudIndex << std::endl;
            std::cout << "num poses: " << octree.numPositionIndexs(node) << std::endl;
            octree.points(pointIndex, 0);
            std::cout << "leaf node index: " << node << std::endl;
            double nodePointX = octree.GetNodePointX(node, 0);
            double nodePointY = octree.GetNodePointY(node, 0);
            double nodePointZ = octree.GetNodePointZ(node, 0);

            std::cout << "node point: " << nodePointX << ", " << nodePointY << ", " << nodePointZ << std::endl;*/

        }


        // generating a new distance function
        
        // putting all the points of the point cloud into chunks
        ChunkGrid chunkGrid = GenerateChunks(pointCloud, CHUNK_SIZE, pointCloudIndex);

        // creating the octree for the data
        /*Octree::Octree pointCloudOctree = Octree::Octree (pointCloud, pointCloudIndex, samplingSpaceX, samplingSpaceY, samplingSpaceZ, SAMPLING_SPACE_OFFSET[0], SAMPLING_SPACE_OFFSET[1], SAMPLING_SPACE_OFFSET[2], MAX_OCTREE_DEPTH, OCTREE_POINT_BUFFER_SIZE, 0.05);
        pointCloudOctree.SubDivide();  // sorting all the point cloud points into the octree and generating the various nodes

        // generating the traversal table
        std::cout << "creating traversal table cache" << std::endl;

        pointCloudOctree.GenerateNodeReferenceChache();  // IT SEEMS TO WORK!!!!!!!

        std::cout << "created traversal table cache" << std::endl;

        CArray <double> distanceFieldOctree = CArray <double> (pointCloudOctree.GetSize());
        
        const CArray <double> pos = pointCloudOctree.GetLeafNodePosition(pointCloudOctree.GetLeafIndex(23.456, 75.346, 14.3456));//25.3467, 28.0, 36.2356));

        // make the code return the min distance if the queue runs out
        //std::cout << pointCloudOctree.NearestNeighborSearch           (7.970978545547169, 31.27784731879688, 7.970978545547169) << std::endl;  // pos(0), pos(1), pos(2)
        std::cout << pointCloudOctree.ApproximateNearestNeighborSearch(7.970978545547169, 31.27784731879688, 7.970978545547169) << std::endl;  // pos(0), pos(1), pos(2)  // seems to work fine for node positions, so it might work fine for this application, although the error should be figured out
        std::cout << chunkGrid.       FindNearestPoint                (7.970978545547169, 31.27784731879688, 7.970978545547169) << std::endl;  // pos(0), pos(1), pos(2)
        */
        // something weird is happening. It works fine above, but once run multiple times below it crashes by the queue running dry. What could cause this discrepancy?
        // cache is almost correct I think, it works at the bounds, although it is still a bit off in those areas, not sure why

        /*
        auto t1 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 1000000; i++) {pointCloudOctree.NearestNeighborSearch(pos(0), pos(1), pos(2));}
        auto t2 = std::chrono::high_resolution_clock::now();
        auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::chrono::duration <double, std::milli> ms_double = t2 - t1;
        std::cout << ms_double.count() << "ms\n";

        t1 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 1000000; i++) {pointCloudOctree.ApproximateNearestNeighborSearch(pos(0), pos(1), pos(2));}
        t2 = std::chrono::high_resolution_clock::now();
        ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        ms_double = t2 - t1;
        std::cout << ms_double.count() << "ms\n";

        t1 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 1000000; i++) {chunkGrid.FindNearestPoint(pos(0), pos(1), pos(2));}
        t2 = std::chrono::high_resolution_clock::now();
        ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        ms_double = t2 - t1;
        std::cout << ms_double.count() << "ms\n";*/


        // the cached neighbors might be broken at the very edges of the bounds of the octree
        // it may need to ascend all the way to the root node, allow it to go heigher
        // the cache seems to be fixed after going to depth = 0 although the edges still crash and have an issue


        /*double px = -100;//pointCloud(25, 0) - 2;
        double py = -100;//pointCloud(25, 1) + 3;
        double pz = -100;//pointCloud(25, 2) + 7;

        std::cout << "nearest point cloud search: " << pointCloudOctree.NearestNeighborSearch(px, py, pz) << std::endl;
        std::cout << "nearest spatial partition search: " << chunkGrid.FindNearestPoint(px, py, pz) << std::endl;*/

        /*std::time_t start, end;
        std::time(&start);
        for (int i = 0; i < 1000000; i++) pointCloudOctree.NearestNeighborSearch(px, py, pz);
        std::time(&end);
        std::cout <<  double(end-start) << "       " << pointCloudOctree.NearestNeighborSearch(px, py, pz) << std::endl;*/

        /*std::time(&start);
        for (int i = 0; i < 1000000; i++) chunkGrid.FindNearestPoint(px, py, pz);
        std::time(&end);
        std::cout <<  double(end-start) << "       " << chunkGrid.FindNearestPoint(px, py, pz) << std::endl;*/

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
                    double distance = chunkGrid.FindNearestPoint(x, y, z);  // pointCloudOctree.NearestNeighborSearch(x, y, z);//

                    //std::cout << "ready to start neighbor search..." << std::endl;
                    //double pointDst = pointCloudOctree.ApproximateNearestNeighborSearch(x, y, z);  // ???????
                    //std::cout << "searched neighbor" << std::endl;
                    /*if (chunkGrid.FindNearestPoint(x, y, z) != pointDst)
                    {
                        //std::cout << "octree nns failed: \n    x, y, z: " << x << ", " << y << ", " << z << "\n    Octree Distance: " << pointCloudOctree.ApproximateNearestNeighborSearch(x, y, z) << "\n    Spatial Distance: " << distance << std::endl;
                    }*/

                    distanceField(x_, y_, z_) = distance;
                }
            }
        }
        
        std::cout << "found all distances" << std::endl;

        // save the distance field (too lazy, add later)

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
            CArray <int8_t> signedGrid = CArray <int8_t> (SAMPLING_SPACE_SIZE[0], SAMPLING_SPACE_SIZE[1], SAMPLING_SPACE_SIZE[2]);
            //signedGrid.set_values(0);

            std::cout << "ready to calculate signs" << std::endl;

            // calculating the signs for the objects and distance field
            CalculateSigns(distanceField, signedGrid);

            std::cout << "signs calculated" << std::endl;

            // an array for all surface points on the grid
            CArray <double> surfacePoints = CArray <double> (SURFACE_POINTS_BUFFER_SIZE, 3);
            int surfaceIndex = 0;  // the current index for the surface points array

            std::cout << "calculating surface points" << std::endl;

            // itterating over every point and calculating if it's a surface point or not
            for (int x = 1; x < SAMPLING_SPACE_SIZE[0] - 1; x++)
            {
                for (int y = 1; y < SAMPLING_SPACE_SIZE[1] - 1; y++)
                {
                    for (int z = 1; z < SAMPLING_SPACE_SIZE[2] - 1; z++)
                    {
                        // making sure there's a hollow point next to the cell
                        if (
                            signedGrid(x - 1, y, z) > 0 || signedGrid(x + 1, y, z) > 0 ||
                            signedGrid(x, y - 1, z) > 0 || signedGrid(x, y + 1, z) > 0 ||
                            signedGrid(x, y, z - 1) > 0 || signedGrid(x, y, z + 1) > 0
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

            std::cout << "generating chunks; Number of points: " << surfaceIndex << std::endl;

            // sorting the points into chunks
            ChunkGrid chunkedSurfacePoints = GenerateChunks(surfacePoints, SHELL_CHUNK_SIZE, surfaceIndex);

            // creating the octree for the data
            //Octree::Octree surfaceOctree = Octree::Octree (surfacePoints, surfaceIndex, SAMPLING_SPACE_SIZE[0], SAMPLING_SPACE_SIZE[1], SAMPLING_SPACE_SIZE[2], SAMPLING_SPACE_OFFSET[0], SAMPLING_SPACE_OFFSET[1], SAMPLING_SPACE_OFFSET[2], MAX_OCTREE_DEPTH, OCTREE_POINT_BUFFER_SIZE, 0.2);
            //surfaceOctree.SubDivide();  // sorting all the point cloud points into the octree and generating the various cells

            std::cout << "chunked the surface points" << std::endl;

            // correcting the iso contour sense the grid spacing has changed for this next step
            double oldIsoContourLevel;
            if (AUTO_SET && !LOADING_DISTANCE_FIELD_SAVE) {
                oldIsoContourLevel = ISO_CONTOUR_LEVEL;
                ISO_CONTOUR_LEVEL = ISO_CONTOUR_LEVEL * 3 / (1/INVERSE_DELTA_X + 1/INVERSE_DELTA_Y + 1/INVERSE_DELTA_Z);
            }

            // finding all the new distances
            for (int x = 0; x < SAMPLING_SPACE_SIZE[0]; x++)
            {
                for (int y = 0; y < SAMPLING_SPACE_SIZE[1]; y++)
                {
                    for (int z = 0; z < SAMPLING_SPACE_SIZE[2]; z++)
                    {
                        // finding the distance to the nearest point
                        double distance = chunkedSurfacePoints.FindNearestPoint(x, y, z);  // surfaceOctree.NearestNeighborSearch(x, y, z);//

                        // finding the correct sign for the current position
                        short sign = signedGrid(x, y, z);
                        if (distanceField(x, y, z) < oldIsoContourLevel && sign > -1) sign = -1;
                        if (!sign) sign = 1;

                        // updating the signed distance funciton with the distance and sign for the given position
                        distanceField(x, y, z) = distance * (double) sign + ISO_CONTOUR_LEVEL;
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
            CArray <int8_t> signedGrid = CArray <int8_t> (SAMPLING_SPACE_SIZE[0], SAMPLING_SPACE_SIZE[1], SAMPLING_SPACE_SIZE[2]);

            // calculating the signs for the objects and distance field
            CalculateSigns(distanceField, signedGrid);

            // itterating through all points and correcting the signs
            for (int x = 0; x < SAMPLING_SPACE_SIZE[0]; x++)
            {
                for (int y = 0; y < SAMPLING_SPACE_SIZE[1]; y++)
                {
                    for (int z = 0; z < SAMPLING_SPACE_SIZE[2]; z++)
                    {
                        distanceField(x, y, z) = distanceField(x, y, z) * (double) signedGrid(x, y, z);
                    }
                }
            }
        }

    }


    // running the code through marching cubes & returning it as an stl file
    MarchingCubes::MarchingCubesToSTL(distanceField, ISO_CONTOUR_LEVEL, 1.0, 1.0, 1.0);  // 1.0/INVERSE_DELTA_X, 1.0/INVERSE_DELTA_Y, 1.0/INVERSE_DELTA_Z (the distance field at this point has a scale represent 1x1x1 cells)


    return 0;
}


