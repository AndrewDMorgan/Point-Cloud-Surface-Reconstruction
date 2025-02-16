#include <stdio.h>
#include "matar.h"


namespace Octree
{

    using namespace mtr;


    // manages and stores an otree data structure

    template<typename T>  // the type for the internal data for each node to store/calculate; this allows a more efficent representation/calculation of specific functions within a specified domain
    class Octree
    {

    private:

        CArray <int> childPointReferences;  // each cell has 8 values each being an index to the same array of the child
        CArray <int> positionIndexesPlusOne;  // (subtract 1 from the index to get the array index) the index of a point within a 1d array of the points being stored in the octree
        CArray <int> numPositionIndexs;  // the number of positions in each given grid cell (esentially a stack)

        CArray <T> gridValues;  // the value for each node; calculated for a storage efficent data storage structure

        int numberChildReferences;

        CArray <double> points;  // the array of points being sorted into the grid
        int numPoints;

        // the position and size of the bounding box of the octree
        double offsetX, offsetY, offsetZ;
        double sizeX, sizeY, sizeZ;

        int rootNodeReferenceIndex;
        int maxDepth;  // the max subdivision depth

        CArray <double> depthSizeScalars;
        CArray <int> depthIndexBufferSearch;

        // the offset orders for the children
        int childOffsetOrderX[8] = {1, 1, 0, 0, 1, 1, 0, 0};
        int childOffsetOrderY[8] = {0, 1, 1, 0, 0, 1, 1, 0};
        int childOffsetOrderZ[8] = {0, 0, 0, 1, 1, 1, 1, 0};

        CArray <int> offsetOrderChildIndex;
        int lastLeafDepth;

        // the constructor
        public: Octree (CArray <double> &_points, int _numPoints, double _sizeX, double _sizeY, double _sizeZ, double _offsetX, double _offsetY, double _offsetZ, int _maxDepth, int bufferSize)
        {
            
            // is the maximum possible cells 4^depth?   or is it some sort of summation?

            // creating the point reference array that stores the octree
            int max1DDepth = pow(4.0, (double) _maxDepth);
            childPointReferences = CArray <int> (max1DDepth, 8);
            gridValues = CArray <T> (max1DDepth);

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

            depthIndexBufferSearch = CArray <int> (maxDepth);

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

            // pre-defining/allocating some variables
            int childOffsetX;
            int childOffsetY;
            int childOffsetZ;
            int childIndex;

            double widthScalar, cellSizeX, cellSizeY, cellSizeZ;

            depthIndexBufferSearch(0) = nodeIndex;

            // itterating through the levels of the tree till the leaf cell is found
            for (int i = 1; i < maxDepth; i++)  // should leave once it's hit the depth limit
            {
                lastLeafDepth = i;
                depthIndexBufferSearch(i) = nodeIndex;

                // getting the cell dimensions
                widthScalar = depthSizeScalars(i);
                cellSizeX = sizeX / widthScalar;
                cellSizeY = sizeY / widthScalar;
                cellSizeZ = sizeZ / widthScalar;

                // getting the offset for the child
                childOffsetX = (int) ((posX - baseX) / cellSizeX);
                childOffsetY = (int) ((posY - baseY) / cellSizeY);
                childOffsetZ = (int) ((posZ - baseZ) / cellSizeZ);

                childIndex = offsetOrderChildIndex(childOffsetX, childOffsetY, childOffsetZ);  // the child's index
                
                lastNodeIndex = nodeIndex;
                nodeIndex = childPointReferences(nodeIndex, childIndex);  // the new node index
                if (nodeIndex < 0) return lastNodeIndex;  // making sure the node won't go off of an empty cell

                // adjusting the base position for the corner of the cell
                baseX += cellSizeX * (double) childOffsetX;
                baseY += cellSizeY * (double) childOffsetY;
                baseZ += cellSizeZ * (double) childOffsetZ;

                // checking if the search has concluded
                if (numPositionIndexs(nodeIndex)) return nodeIndex;  // checking if a value has been added -- only leaf nodes contain values

            }

            return nodeIndex;  // returning the found node
        }


        // gets the nearest neighbor but also gets the base leaf node
        public: double NearestNeighborSearch (double samplePositionX, double samplePositionY, double samplePositionZ)
        {
            // getting the leaf node to begin the search
            int leafNodeIndex = GetLeafIndex(samplePositionX, samplePositionY, samplePositionZ);

            return NearestNeighborSearch(leafNodeIndex, samplePositionX, samplePositionY, samplePositionZ);  // returning the distance
        }

        // preforms a nearest neighbor search based on a given point
        public: double NearestNeighborSearch (int leafNodeIndex, double samplePositionX, double samplePositionY, double samplePositionZ)
        {
            double distance = 99999999.0;  // the minimum distance found
            
            // checking the initial point
            BranchSearch(depthIndexBufferSearch(lastLeafDepth), -1, distance, samplePositionX, samplePositionY, samplePositionZ);
            
            // itterating over the neighboring cells
            int currentDepth = lastLeafDepth;
            for (int i = 0; i < lastLeafDepth; i++)
            {
                if (distance < 99999999.0) return sqrt(distance); // if any point is found the distance is returned

                // starting a new search
                currentDepth--;
                BranchSearch(depthIndexBufferSearch(currentDepth), depthIndexBufferSearch(currentDepth + 1), distance, samplePositionX, samplePositionY, samplePositionZ);
            }

            // returning the distance
            return sqrt(distance);
        }

        // searching down the entirety of a branch
        private: void BranchSearch (int startingIndex, int voidedIndex, double &distance, double &samplePositionX, double &samplePositionY, double &samplePositionZ)
        {
            // going to each child and casting a new branch
            for (int childIndex = 0; childIndex < 8; childIndex++)
            {
                int newChildIndex = childPointReferences(startingIndex, childIndex);
                if (newChildIndex == voidedIndex) continue;  // making sure the point hasn't already been searched

                // checking if the point is an end point
                if (newChildIndex == -1)
                {
                    // itterating through all points, getting the distances, and ending the cycle
                    for (int pointId = 0; pointId < numPositionIndexs(startingIndex); pointId++)
                    {
                        // calculating the distance
                        double dx = samplePositionX - points(positionIndexesPlusOne(startingIndex, pointId) - 1, 0);
                        double dy = samplePositionY - points(positionIndexesPlusOne(startingIndex, pointId) - 1, 1);
                        double dz = samplePositionZ - points(positionIndexesPlusOne(startingIndex, pointId) - 1, 2);
                        double newDistance = dx*dx + dy*dy + dz*dz;

                        // checking if the distance is a new minimum
                        if (newDistance < distance) distance = newDistance;
                    }

                    return;  // ending this branch
                }

                // starting a further search
                BranchSearch(newChildIndex, -99999, distance, samplePositionX, samplePositionY, samplePositionZ);
            }

            return;
        }


        // gets a specific grid values adress
        public: T GetGridValue (int index) {return &gridValues(index);}

        // gets the point at a given node index
        public: double GetNodePointX (int nodeIndex, int bufferIndex) {return points(positionIndexesPlusOne(nodeIndex, bufferIndex) - 1, 0);}
        public: double GetNodePointY (int nodeIndex, int bufferIndex) {return points(positionIndexesPlusOne(nodeIndex, bufferIndex) - 1, 1);}
        public: double GetNodePointZ (int nodeIndex, int bufferIndex) {return points(positionIndexesPlusOne(nodeIndex, bufferIndex) - 1, 2);}

    };
}

