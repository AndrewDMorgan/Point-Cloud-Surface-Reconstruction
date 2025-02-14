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

    }

    // subdivides the octree repeatidly
    public: void SubDivide()
    {
        // starting the subdivision and saving the final index to act as the inital root node
        rootNodeReferenceIndex = SubDivideRecursive(0, 0, 0, offsetX, offsetY, offsetZ, 1);
    }

    // subdivides the octree to fit the points      returns the reference index
    private: int SubDivideRecursive (int shiftX, int shiftY, int shiftZ, double tilePositionX, double tilePositionY, double tilePositionZ, int depth)
    {

        // getting the current size
        double scaledDepth = pow(2, depth - 1);    // 2^x is the scaling fator for size    it goes /1 /2 /4 /8 /16...
        double tileSizeX = sizeX / scaledDepth;
        double tileSizeY = sizeY / scaledDepth;
        double tileSizeZ = sizeZ / scaledDepth;

        // getting the current position
        tilePositionX += tileSizeX * shiftX;
        tilePositionY += tileSizeY * shiftY;
        tilePositionZ += tileSizeZ * shiftZ;

        //std::cout << "entering recursive itteration   0:1" << std::endl;

        // checking if the max depth was reached
        if (depth == maxDepth)
        {
            std::cout << "max depth, adding child   1:1" << std::endl;
            childPointReferences(numberChildReferences, 0) = -1;
            childPointReferences(numberChildReferences, 1) = -1;
            childPointReferences(numberChildReferences, 2) = -1;
            childPointReferences(numberChildReferences, 3) = -1;
            childPointReferences(numberChildReferences, 4) = -1;  // no references available
            childPointReferences(numberChildReferences, 5) = -1;
            childPointReferences(numberChildReferences, 6) = -1;
            childPointReferences(numberChildReferences, 7) = -1;
            std::cout << "added child references   1:2" << std::endl;

            // adding any possible remaining points
            for (int i = 0; i < points.dims(0); i++)
            {
                // checking if the point lays within the bounding box of the node
                if (
                    points(i, 0) >= tilePositionX && points(i, 0) < tilePositionX + tileSizeX &&
                    points(i, 1) >= tilePositionY && points(i, 1) < tilePositionY + tileSizeY &&
                    points(i, 2) >= tilePositionZ && points(i, 2) < tilePositionZ + tileSizeZ
                ) {
                    
                    std::cout << "adding the point to buffer   1:3" << std::endl;
                    std::cout << numberChildReferences << ", " << numPositionIndexs(numberChildReferences) << std::endl;
                    // adding the point
                    positionIndexesPlusOne(numberChildReferences, numPositionIndexs(numberChildReferences)) = i + 1;
                    numPositionIndexs(numberChildReferences)++;
                }
            }
            std::cout << "added points   1:4       " << numPositionIndexs(numberChildReferences) << std::endl;

            // returning the index
            numberChildReferences++;
            return numberChildReferences - 1;
        }
        // starting at the root node, and dividing each cell until it no longer has a value or has only one value
        
        //std::cout << "searching through points   0:2" << std::endl;

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
                
                lastIndex = i;  // if there's only 1 point this is chosen
                numberBoundingPoints++;
            }
        }

        //std::cout << "checking if there aren't any points   0:3" << std::endl;
        
        // checking if there's no points in the cell
        if (!numberBoundingPoints)
        {
            std::cout << "None left: " << numberChildReferences << "     Depth: " << depth << std::endl;

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

        //std::cout << "checking if there is only a single point   0:4" << std::endl;
        
        // checking if there's a single point in the cell
        if (numberBoundingPoints == 1)
        {
            std::cout << "One left: " << numberChildReferences << "     Depth: " << depth << std::endl;

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

        //std::cout << "getting new references   0:5" << std::endl;

        // continuing to sub divide into more cells

        // getting the references
        int childIndex1 = SubDivideRecursive(1, 0, 0, tilePositionX, tilePositionY, tilePositionZ, depth + 1);
        int childIndex2 = SubDivideRecursive(1, 1, 0, tilePositionX, tilePositionY, tilePositionZ, depth + 1);
        int childIndex3 = SubDivideRecursive(0, 1, 0, tilePositionX, tilePositionY, tilePositionZ, depth + 1);
        int childIndex4 = SubDivideRecursive(0, 0, 1, tilePositionX, tilePositionY, tilePositionZ, depth + 1);
        int childIndex5 = SubDivideRecursive(1, 0, 1, tilePositionX, tilePositionY, tilePositionZ, depth + 1);
        int childIndex6 = SubDivideRecursive(1, 1, 1, tilePositionX, tilePositionY, tilePositionZ, depth + 1);
        int childIndex7 = SubDivideRecursive(0, 1, 1, tilePositionX, tilePositionY, tilePositionZ, depth + 1);
        int childIndex8 = SubDivideRecursive(0, 0, 0, tilePositionX, tilePositionY, tilePositionZ, depth + 1);

        std::cout << "adding child references: " << numberChildReferences << std::endl;

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

};


