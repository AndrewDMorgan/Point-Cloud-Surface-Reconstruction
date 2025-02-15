#include <stdio.h>
#include "matar.h"

using namespace mtr;


// from https://www.geeksforgeeks.org/program-for-factorial-of-a-number/
int Factorial(int n)
{
    int res = 1;
    for (int i = 2; i <= n; i++)
        res *= i;
    return res;
}



// manages and stores the otree data structure
class Octree
{

private:

    CArray <int> childPointReferences;  // each cell has 8 values each being an index to the same array of the child
    CArray <int> positionIndexesPlusOne;  // (subtract 1 from the index to get the array index) the index of a point within a 1d array of the points being stored in the octree
    CArray <int> numPositionIndexs;  // the number of positions in each given grid cell (esentially a stack)

    int numberChildReferences;

    CArray <double> points;  // the array of points being sorted into the grid
    int numPoints;

    // the position and size of the bounding box of the octree
    double offsetX, offsetY, offsetZ;
    double sizeX, sizeY, sizeZ;

    int rootNodeReferenceIndex;
    int maxDepth;  // the max subdivision depth

    CArray <double> depthSizeScalars;

    // the offset orders for the children
    int childOffsetOrderX[8] = {1, 1, 0, 0, 1, 1, 0, 0};
    int childOffsetOrderY[8] = {0, 1, 1, 0, 0, 1, 1, 0};
    int childOffsetOrderZ[8] = {0, 0, 0, 1, 1, 1, 1, 0};

    CArray <int> offsetOrderChildIndex;

    // the constructor
    public: Octree (CArray <double> &_points, int _numPoints, double _sizeX, double _sizeY, double _sizeZ, double _offsetX, double _offsetY, double _offsetZ, int _maxDepth, int bufferSize)
    {
        
        // is the maximum possible cells 4^depth?   or is it some sort of summation?

        // creating the point reference array that stores the octree
        int max1DDepth = pow(4.0, (double) _maxDepth);//Factorial(_maxDepth);
        childPointReferences = CArray <int> (max1DDepth, 8);

        // creating the point reference array
        positionIndexesPlusOne = CArray <int> (max1DDepth, bufferSize);
        numPositionIndexs = CArray <int> (max1DDepth);
        numPositionIndexs.set_values(0);  // making sure the default is 

        maxDepth = _maxDepth;
        points = _points;  // storing the address of the points inside the class
        numPoints = _numPoints;

        // setting the bounds of the octree
        offsetX = _offsetX;
        offsetY = _offsetY;
        offsetZ = _offsetZ;

        sizeX = _sizeX;
        sizeY = _sizeY;
        sizeZ = _sizeZ;

        numberChildReferences = 0;

        // pre-calculating all the scalars for the different depth levels (widths of the cells)
        depthSizeScalars = CArray <double> (maxDepth);
        for (int depth = 0; depth < maxDepth; depth++)
        {
            depthSizeScalars(depth) = pow(2, (double) depth);
        }

        // getting the reverse child offset order (for quick reverse look ups)
        offsetOrderChildIndex = CArray <int> (8, 8, 8);
        offsetOrderChildIndex(1, 0, 0) = 0;
        offsetOrderChildIndex(1, 1, 0) = 1;
        offsetOrderChildIndex(0, 1, 0) = 2;
        offsetOrderChildIndex(0, 0, 1) = 3;
        offsetOrderChildIndex(1, 0, 1) = 4;
        offsetOrderChildIndex(1, 1, 1) = 5;
        offsetOrderChildIndex(0, 1, 1) = 6;
        offsetOrderChildIndex(0, 0, 0) = 7;

    }

    // subdivides the octree repeatidly
    public: void SubDivide()
    {
        // starting the subdivision and saving the final index to act as the inital root node
        rootNodeReferenceIndex = SubDivideRecursive(0, 0, 0, offsetX, offsetY, offsetZ, 1);
        std::cout << "num children ->>>>  " << numberChildReferences << std::endl;
    }

    // subdivides the octree to fit the points      returns the reference index
    private: int SubDivideRecursive (int shiftX, int shiftY, int shiftZ, double tilePositionX, double tilePositionY, double tilePositionZ, int depth)
    {

        // getting the current size
        double scaledDepth = depthSizeScalars(depth - 1);    // 2^x is the scaling fator for size    it goes /1 /2 /4 /8 /16...
        double tileSizeX = sizeX / scaledDepth;
        double tileSizeY = sizeY / scaledDepth;
        double tileSizeZ = sizeZ / scaledDepth;

        // getting the current position
        tilePositionX += tileSizeX * shiftX;
        tilePositionY += tileSizeY * shiftY;
        tilePositionZ += tileSizeZ * shiftZ;

        // checking if the max depth was reached
        if (depth == maxDepth)
        {
            childPointReferences(numberChildReferences, 0) = -1;
            childPointReferences(numberChildReferences, 1) = -1;
            childPointReferences(numberChildReferences, 2) = -1;
            childPointReferences(numberChildReferences, 3) = -1;
            childPointReferences(numberChildReferences, 4) = -1;  // no references available
            childPointReferences(numberChildReferences, 5) = -1;
            childPointReferences(numberChildReferences, 6) = -1;
            childPointReferences(numberChildReferences, 7) = -1;

            // adding any possible remaining points
            for (int i = 0; i < points.dims(0); i++)
            {
                // checking if the point lays within the bounding box of the node
                if (
                    points(i, 0) >= tilePositionX && points(i, 0) < tilePositionX + tileSizeX &&
                    points(i, 1) >= tilePositionY && points(i, 1) < tilePositionY + tileSizeY &&
                    points(i, 2) >= tilePositionZ && points(i, 2) < tilePositionZ + tileSizeZ
                ) {
                    
                    // adding the point
                    positionIndexesPlusOne(numberChildReferences, numPositionIndexs(numberChildReferences)) = i + 1;
                    
                    numPositionIndexs(numberChildReferences)++;
                }
            }

            // returning the index
            numberChildReferences++;
            return numberChildReferences - 1;
        }
        // starting at the root node, and dividing each cell until it no longer has a value or has only one value
        
        // checking how many points are within the bounds of the cell
        int lastIndex = -1;
        int numberBoundingPoints = 0;
        for (int i = 0; i < points.dims(0); i++)
        {
            // checking if the point lays within the bounding box of the node
            if (
                points(i, 0) >= tilePositionX && points(i, 0) < tilePositionX + tileSizeX &&
                points(i, 1) >= tilePositionY && points(i, 1) < tilePositionY + tileSizeY &&
                points(i, 2) >= tilePositionZ && points(i, 2) < tilePositionZ + tileSizeZ
            ) {
                
                lastIndex = i;  // for if there's only 1 point this is chosen
                numberBoundingPoints++;
            }
        }

        // checking if there's no points in the cell
        if (!numberBoundingPoints)
        {
            childPointReferences(numberChildReferences, 0) = -1;
            childPointReferences(numberChildReferences, 1) = -1;
            childPointReferences(numberChildReferences, 2) = -1;
            childPointReferences(numberChildReferences, 3) = -1;
            childPointReferences(numberChildReferences, 4) = -1;  // no references available
            childPointReferences(numberChildReferences, 5) = -1;
            childPointReferences(numberChildReferences, 6) = -1;
            childPointReferences(numberChildReferences, 7) = -1;

            // returning the index
            numberChildReferences++;
            return numberChildReferences - 1;
        }

        // checking if there's a single point in the cell
        if (numberBoundingPoints == 1)
        {
            childPointReferences(numberChildReferences, 0) = -1;
            childPointReferences(numberChildReferences, 1) = -1;
            childPointReferences(numberChildReferences, 2) = -1;
            childPointReferences(numberChildReferences, 3) = -1;
            childPointReferences(numberChildReferences, 4) = -1;  // no references available
            childPointReferences(numberChildReferences, 5) = -1;
            childPointReferences(numberChildReferences, 6) = -1;
            childPointReferences(numberChildReferences, 7) = -1;

            // adding the point
            positionIndexesPlusOne(numberChildReferences, numPositionIndexs(numberChildReferences)) = lastIndex + 1;
            numPositionIndexs(numberChildReferences)++;

            // returning the index
            numberChildReferences++;
            return numberChildReferences - 1;
        }

        // continuing to sub divide into more cells

        // getting the references
        
        int childIndex1 = SubDivideRecursive(childOffsetOrderX[0], childOffsetOrderY[0], childOffsetOrderZ[0], tilePositionX, tilePositionY, tilePositionZ, depth + 1);
        int childIndex2 = SubDivideRecursive(childOffsetOrderX[1], childOffsetOrderY[1], childOffsetOrderZ[1], tilePositionX, tilePositionY, tilePositionZ, depth + 1);
        int childIndex3 = SubDivideRecursive(childOffsetOrderX[2], childOffsetOrderY[2], childOffsetOrderZ[2], tilePositionX, tilePositionY, tilePositionZ, depth + 1);
        int childIndex4 = SubDivideRecursive(childOffsetOrderX[3], childOffsetOrderY[3], childOffsetOrderZ[3], tilePositionX, tilePositionY, tilePositionZ, depth + 1);
        int childIndex5 = SubDivideRecursive(childOffsetOrderX[4], childOffsetOrderY[4], childOffsetOrderZ[4], tilePositionX, tilePositionY, tilePositionZ, depth + 1);
        int childIndex6 = SubDivideRecursive(childOffsetOrderX[5], childOffsetOrderY[5], childOffsetOrderZ[5], tilePositionX, tilePositionY, tilePositionZ, depth + 1);
        int childIndex7 = SubDivideRecursive(childOffsetOrderX[6], childOffsetOrderY[6], childOffsetOrderZ[6], tilePositionX, tilePositionY, tilePositionZ, depth + 1);
        int childIndex8 = SubDivideRecursive(childOffsetOrderX[7], childOffsetOrderY[7], childOffsetOrderZ[7], tilePositionX, tilePositionY, tilePositionZ, depth + 1);

        // saving the points
        childPointReferences(numberChildReferences, 0) = childIndex1;
        childPointReferences(numberChildReferences, 1) = childIndex2;
        childPointReferences(numberChildReferences, 2) = childIndex3;
        childPointReferences(numberChildReferences, 3) = childIndex4;
        childPointReferences(numberChildReferences, 4) = childIndex5;
        childPointReferences(numberChildReferences, 5) = childIndex6;
        childPointReferences(numberChildReferences, 6) = childIndex7;
        childPointReferences(numberChildReferences, 7) = childIndex8;

        // incrementing the number of references
        numberChildReferences++;

        return numberChildReferences - 1;  // returning the index of the child node ( - 1 because the index was incremented just before)
    
    }


    // grabs the leaf cell for a given position
    public: int GetLeafIndex (double posX, double posY, double posZ)
    {
        int nodeIndex = rootNodeReferenceIndex;  // getting the root node
        int lastNodeIndex;

        // tracking the current position
        double baseX = offsetX;
        double baseY = offsetY;
        double baseZ = offsetZ;

        std::cout << "position: " << baseX << ", " << baseY << ", " << baseZ << std::endl;

        double widthScalar, cellSizeX, cellSizeY, cellSizeZ;

        // itterating through the levels of the tree till the leaf cell is found
        for (int i = 1; i < maxDepth; i++)  // should leave once it's hit the depth limit
        {
            // getting the cell dimensions
            widthScalar = depthSizeScalars(i);
            cellSizeX = sizeX / widthScalar;
            cellSizeY = sizeY / widthScalar;
            cellSizeZ = sizeZ / widthScalar;

            // getting the offset for the child
            int childOffsetX = (int) ((posX - baseX) / cellSizeX);
            int childOffsetY = (int) ((posY - baseY) / cellSizeY);
            int childOffsetZ = (int) ((posZ - baseZ) / cellSizeZ);

            std::cout << childOffsetX << ", " << childOffsetY << ", " << childOffsetZ << std::endl;

            int childIndex = offsetOrderChildIndex(childOffsetX, childOffsetY, childOffsetZ);  // the child's index
            
            lastNodeIndex = nodeIndex;
            nodeIndex = childPointReferences(nodeIndex, childIndex);  // the new node index
            if (nodeIndex < 0) return lastNodeIndex;  // making sure the node won't go off of an empty cell

            // adjusting the base position for the corner of the cell
            baseX += cellSizeX * (double) childOffsetX;
            baseY += cellSizeY * (double) childOffsetY;
            baseZ += cellSizeZ * (double) childOffsetZ;

            std::cout << "position: " << baseX << ", " << baseY << ", " << baseZ << std::endl;

            // checking if the search has concluded
            if (numPositionIndexs(nodeIndex)) return nodeIndex;  // checking if a value has been added -- only leaf nodes contain values

        }

        return nodeIndex;  // returning the found node
    }


    // gets the point at a given node index
    public: double GetNodePointX (int nodeIndex, int bufferIndex) {return points(positionIndexesPlusOne(nodeIndex, bufferIndex) - 1, 0);}
    public: double GetNodePointY (int nodeIndex, int bufferIndex) {return points(positionIndexesPlusOne(nodeIndex, bufferIndex) - 1, 1);}
    public: double GetNodePointZ (int nodeIndex, int bufferIndex) {return points(positionIndexesPlusOne(nodeIndex, bufferIndex) - 1, 2);}

};


