#include <iostream>
#include "matar.h"

// the maximum depth a binary search can go to
const int MAX_BINARY_SEARCH_ITTERATIONS = 250;

/*

create a hashed-grid based on the leaf positions relative to the x/y/z coords for that portion of the octree
use this to search as it's space is still reduced but also easier to acsess

how do you find the x/y/z coords of a global offset of positions?

info would need to be stored about the intial offsets of each row/collumn/slice and how large each cell is so that it can be undone in order to back trace to get to a neighboring cell


use some form of hashing with the leaf nodes to find a better/faster nearest neighbor search
https://link.springer.com/article/10.1007/s00138-017-0889-4








create an array with pointers to all neighboring cells

in search function:
    use another array that store whenever a cell is visisted during a nn search.
    traverse the array counting the distance across each node (figure out that too) to find the current distance to stop at the limit when expanding (whenever a point is found, the distance can simply be shrunk until no remaining cells exsist. This can all be done recursively)


to find a neighbor node, could you just find the coord of the center of the face of the urrent node and shift it by some small value and solve for a leaf node?
    How do you deal with finding multiple neighbors? do you find the larger cell of the same depth than add all leaf nodes bellow?
    
    Go to the corner of the cell plus a tiny offset (small enough it can't jump more than a cell)
    go down until equal to the current depth
    find all leaf nodes falling under (that still touch the face on at least one side) that node and add them

*/


namespace Octree
{

    using namespace mtr;


    // a binary search algerithm (returns the index)
    template <typename T>  // has to be a type with greater than and less than opperators
    int BinarySearch (const CArray <T> &points, const int arraySize, const T searchValue)
    {
        // dividing the space in two until the position is found
        int currentIndex = 0;
        int dividedSize = arraySize;

        // looping till the value is found (or max itterations)
        for (int i = 0; i < MAX_BINARY_SEARCH_ITTERATIONS; i++)
        {

            // dividing the space in two
            int halfWidth = (int) (dividedSize * 0.5);
            dividedSize -= halfWidth;

            T middleValue = points(currentIndex + halfWidth);

            int currentIndexOld = currentIndex;
            // checking if the value is smaller or greater than the input value
            if (middleValue == searchValue) return currentIndex + halfWidth;  // concluding the search once the value is found
            if (middleValue < searchValue) currentIndex += halfWidth;  // checking which half the value is in the move the index; allowing for odd and even sized arrays

        }

        std::cout << "Binary Search Failed--Reached Max Depth: Octree.hpp  line 75 (Octree::BinarySearch)\n    Try increasing max search itterations ( MAX_BINARY_SEARCH_ITTERATIONS ) at top of header file" << std::endl;
        return currentIndex;  // returning the last index (the search ran out of search depth likely)
    }


    // manages and stores an otree data structure
    class Octree
    {

    private:

        CArray <int> childPointReferences;  // each cell has 8 values each being an index to the same array of the child
        CArray <int> positionIndexesPlusOne;  // (subtract 1 from the index to get the array index) the index of a point within a 1d array of the points being stored in the octree
        CArray <int> numPositionIndexs;  // the number of positions in each given grid cell (esentially a stack)

        CArray <double> leafNodePositions;  // the global positions of the leaf nodes (different from the indexes)
        CArray <int> leafNodeIndexes;

        int numberChildReferences;
        int numberOfLeafNodes;

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
        const int childOffsetOrderX[8] = {1, 1, 0, 0, 1, 1, 0, 0};
        const int childOffsetOrderY[8] = {0, 1, 1, 0, 0, 1, 1, 0};
        const int childOffsetOrderZ[8] = {0, 0, 0, 1, 1, 1, 1, 0};

        CArray <int> offsetOrderChildIndex;
        int lastLeafDepth;

        CArray <int> cellNeighborReferences;
        CArray <int> cellNeighborBufferSize;

        
        // the constructors
        public: Octree (const CArray <double> &_points, int _numPoints, double _sizeX, double _sizeY, double _sizeZ, double _offsetX, double _offsetY, double _offsetZ, int _maxDepth, int bufferSize)
        {
            int totalPoints = (int) pow(4.0, (double) _maxDepth);
            OctreeConstructor(_points, _numPoints, _sizeX, _sizeY, _sizeZ, _offsetX, _offsetY, _offsetZ, _maxDepth, bufferSize, totalPoints);
        }

        // scales the memory for better performance (risky, it may not allocate enough, but that's unlikely)
        public: Octree (const CArray <double> &_points, int _numPoints, double _sizeX, double _sizeY, double _sizeZ, double _offsetX, double _offsetY, double _offsetZ, int _maxDepth, int bufferSize, double memoryScalar)
        {
            int totalPoints = (int) (pow(4.0, (double) _maxDepth) * memoryScalar);
            OctreeConstructor(_points, _numPoints, _sizeX, _sizeY, _sizeZ, _offsetX, _offsetY, _offsetZ, _maxDepth, bufferSize, totalPoints);
        }

        // finishes constructing the object/instance
        public: void OctreeConstructor (const CArray <double> &_points, int _numPoints, double _sizeX, double _sizeY, double _sizeZ, double _offsetX, double _offsetY, double _offsetZ, int _maxDepth, int bufferSize, int maxNumberOfPositions)
        {
            
            // creating the point reference array that stores the octree
            int max1DDepth = maxNumberOfPositions;
            childPointReferences = CArray <int> (max1DDepth, 8);

            // creating the point reference array
            positionIndexesPlusOne = CArray <int> (max1DDepth, bufferSize);
            numPositionIndexs = CArray <int> (max1DDepth);
            numPositionIndexs.set_values(0);  // making sure the default is 
            leafNodePositions = CArray <double> (max1DDepth, 3);
            leafNodeIndexes = CArray <int> (max1DDepth);

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
                depthSizeScalars(depth) = 1.0 / pow(2, (double) depth);  // the inverse to reduce the number of division operations
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
            numberOfLeafNodes = 0;  // resetting the number of leaf notes

            // starting the subdivision and saving the final index to act as the inital root node
            rootNodeReferenceIndex = SubDivideRecursive(0, 0, 0, offsetX, offsetY, offsetZ, 1);
        }

        // subdivides the octree to fit the points      returns the reference index
        private: int SubDivideRecursive (int shiftX, int shiftY, int shiftZ, double tilePositionX, double tilePositionY, double tilePositionZ, int depth)
        {
            // getting the current size
            double scaledDepth = depthSizeScalars(depth - 1);    // 2^x is the scaling fator for size    it goes /1 /2 /4 /8 /16...
            double tileSizeX = sizeX * scaledDepth;
            double tileSizeY = sizeY * scaledDepth;
            double tileSizeZ = sizeZ * scaledDepth;

            // getting the current position
            tilePositionX += tileSizeX * shiftX;
            tilePositionY += tileSizeY * shiftY;
            tilePositionZ += tileSizeZ * shiftZ;

            // checking if the max depth was reached
            if (depth == maxDepth)
            {
                // adding the position of the leaf node
                leafNodePositions(numberOfLeafNodes, 0) = tilePositionX + tileSizeX * 0.5;
                leafNodePositions(numberOfLeafNodes, 1) = tilePositionY + tileSizeY * 0.5;  // the position of the center of the current cell
                leafNodePositions(numberOfLeafNodes, 2) = tilePositionZ + tileSizeZ * 0.5;
                leafNodeIndexes(numberOfLeafNodes) = numberChildReferences;
                numberOfLeafNodes++;

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
                    if (numberBoundingPoints > 1) break;  // if there's more than 2 points, the node needs to be divided so no further searching is needed
                }
            }

            // checking if there's no points in the cell
            if (!numberBoundingPoints)
            {
                // adding the position of the leaf node
                leafNodePositions(numberOfLeafNodes, 0) = tilePositionX + tileSizeX * 0.5;
                leafNodePositions(numberOfLeafNodes, 1) = tilePositionY + tileSizeY * 0.5;  // the position of the center of the current cell
                leafNodePositions(numberOfLeafNodes, 2) = tilePositionZ + tileSizeZ * 0.5;
                leafNodeIndexes(numberOfLeafNodes) = numberChildReferences;
                numberOfLeafNodes++;

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
                // adding the position of the leaf node
                leafNodePositions(numberOfLeafNodes, 0) = tilePositionX + tileSizeX * 0.5;
                leafNodePositions(numberOfLeafNodes, 1) = tilePositionY + tileSizeY * 0.5;  // the position of the center of the current cell
                leafNodePositions(numberOfLeafNodes, 2) = tilePositionZ + tileSizeZ * 0.5;
                leafNodeIndexes(numberOfLeafNodes) = numberChildReferences;
                numberOfLeafNodes++;

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


        // generates and chaches pointers (in the form of an index) to all neighboring nodes from each node
        public: void GenerateNodeReferenceChache ()
        {
            // creating the reference array
            int nodeRefferenceIndex = 0;
            cellNeighborReferences = CArray <int> (numberChildReferences, maxDepth*maxDepth*3 + 3);  // pretty largeish array unfortunately
            cellNeighborBufferSize = CArray <int> (numberChildReferences);
            cellNeighborBufferSize.set_values(0);  // the intial buffer position is at 0 for all elements

            // finding all leaf nodes and then finding their neighbors
            NodeReferenceCacheRecursiveAscent(1, rootNodeReferenceIndex, offsetX, offsetY, offsetZ, nodeRefferenceIndex);

            return;  // ending the function
        }

        // recursively jumps down all branches and starts the generation of the leaf nodes' neighbors
        private: void NodeReferenceCacheRecursiveAscent (const int depth, const int childIndex, double posX, double posY, double posZ, int &nodeRefferenceIndex)
        {
            int newPosX, newPosY, newPosZ;

            // going through all children for this node
            for (int i = 0; i < 8; i++)
            {
                // checking if the child is a leaf node
                if (childPointReferences(childIndex, i) == -1) {
                    GenerateLeafNeighbors(depth, childIndex, posX, posY, posZ, nodeRefferenceIndex);  // generating the cells neighbors
                    return;  // the leaf node had all neighbors generated, no need to continue
                }
                
                // finding the position of the child node
                double scaledDepth = depthSizeScalars(depth - 1);
                double nodeSizeX = sizeX * scaledDepth;
                double nodeSizeY = sizeY * scaledDepth;  // I'm guessing tracking the position is faster than a binary search to lookup the positional information
                double nodeSizeZ = sizeZ * scaledDepth;

                // getting the current position
                newPosX = posX + nodeSizeX * (double) childOffsetOrderX[i];
                newPosY = posY + nodeSizeY * (double) childOffsetOrderY[i];
                newPosZ = posZ + nodeSizeZ * (double) childOffsetOrderZ[i];

                // continuning the recursive ascent
                NodeReferenceCacheRecursiveAscent(depth + 1, childPointReferences(childIndex, i), newPosX, newPosY, newPosZ, nodeRefferenceIndex);
            }

            return;  // ending the function
        }

        // generates the neighbors for a given leaf node
        private: void GenerateLeafNeighbors (const int depth, const int childIndex, const double posX, const double posY, const double posZ, int &nodeRefferenceIndex)
        {
            // getting the size of the node
            double scaledDepth = depthSizeScalars(depth - 1);
            double nodeSizeX = sizeX * scaledDepth;
            double nodeSizeY = sizeY * scaledDepth;
            double nodeSizeZ = sizeZ * scaledDepth;

            // finding the center of the leaf node
            double searchPosX = posX + nodeSizeX * 0.5;
            double searchPosY = posY + nodeSizeY * 0.5;  // i'm guessing it's  faster to track
            double searchPosZ = posZ + nodeSizeZ * 0.5;


            // getting the path to this node (so it can be back-traced) along with the shifts to find neighboring sections
            int ascentNodeIndex = rootNodeReferenceIndex;
            int previousAscentNodeIndex;

            double baseX = offsetX;
            double baseY = offsetY;
            double baseZ = offsetZ;

            int childOffsetX, childOffsetY, childOffsetZ;
            double widthScalar, cellSizeX, cellSizeY, cellSizeZ;
            int ascentChildIndex;

            depthIndexBufferSearch(0) = ascentNodeIndex;

            // looping until the inital leaf-node has been found again (to generate the path)
            for (int i = 1; i < depth; i++)  // by using the current depth, no if checks are needed to check for a leaf node sense it's already known this is a leaf node
            {
                lastLeafDepth = i;

                // finding the size of the node
                widthScalar = depthSizeScalars(i);
                cellSizeX = sizeX * widthScalar;
                cellSizeY = sizeY * widthScalar;
                cellSizeZ = sizeZ * widthScalar;

                // getting the current offset
                childOffsetX = (int) ((searchPosX - baseX) / cellSizeX);
                childOffsetY = (int) ((searchPosY - baseY) / cellSizeY);
                childOffsetZ = (int) ((searchPosZ - baseZ) / cellSizeZ);

                ascentChildIndex = offsetOrderChildIndex(childOffsetX, childOffsetY, childOffsetZ);  // the index of the child within a given parent node
                ascentNodeIndex = childPointReferences(ascentNodeIndex, ascentChildIndex);

                // adjusting the base position
                baseX += cellSizeX * (double) childOffsetX;
                baseY += cellSizeY * (double) childOffsetY;
                baseZ += cellSizeZ * (double) childOffsetZ;


                // adding the node to the path
                depthIndexBufferSearch(i) = ascentNodeIndex;
            }


            // finding the neighbors within the child's parent node
                // find the nodes offset in each axis, and do 1 - offset to get the neighboring node
            
            // getting the current offset to find the opposing neighbors
            childOffsetX = (int) ((searchPosX - baseX) / nodeSizeX);
            childOffsetY = (int) ((searchPosY - baseY) / nodeSizeY);
            childOffsetZ = (int) ((searchPosZ - baseZ) / nodeSizeZ);
            
            // adding the nodes based on 1 - the offset for each axis
            cellNeighborReferences(childIndex, 0) = childPointReferences(
                    depthIndexBufferSearch(depth - 1),
                    offsetOrderChildIndex(1 - childOffsetX, childOffsetY, childOffsetZ)
                );  // oposing x axis neighbor ^
            cellNeighborReferences(childIndex, 1) = childPointReferences(
                    depthIndexBufferSearch(depth - 1),
                    offsetOrderChildIndex(childOffsetX, 1 - childOffsetY, childOffsetZ)
                );  // oposing y axis neighbor ^
            cellNeighborReferences(childIndex, 2) = childPointReferences(
                    depthIndexBufferSearch(depth - 1),
                    offsetOrderChildIndex(childOffsetX, childOffsetY, 1 - childOffsetZ)
                );  // oposing z axis neighbor ^
            cellNeighborBufferSize(childIndex) = 3;  // three nodes were added; it's not possible for anything to already be in the buffer


            // finding the other 3 side neighbors
            //

            /*
            what happens if the other cell is smaller????
            */

            // finding the other children of the parten of this node
            /*int node;
            double nodePosX, nodePosY, nodePosZ;
            for (int i = 0; i < 8; i++)  // looping through all 8 children
            {
                node = childPointReferences(depthIndexBufferSearch(depth - 1), i);  // finding the neighboring child
                if (node == childIndex) continue;  // skipping the leaf node that's currently being looked at

                nodePosX = posX + (double) childOffsetOrderX[i] * depthSizeScalars(depth - 2);
                nodePosY = posY + (double) childOffsetOrderY[i] * depthSizeScalars(depth - 2);  // getting the position of the child
                nodePosZ = posZ + (double) childOffsetOrderZ[i] * depthSizeScalars(depth - 2);

                // checking if the point borders this node or not
                if (
                    false
                ) {}
            }*/
            
            // finding the other 3 cells (1 from the larger cube above, 1 to the side, and 1 to the front)
            // how do you back trace until it shifts in the correcct direction? Have a way to see shifts and just look at that?

            // going back down the path and branching out to the neighbors

            /*
            go back down the tree, branch out based on the shift coordinate looking for a valid child until one is found and then...
            decending while staying in line with the face and branching out to any children that fit width wise until at the bottom
            add those leaf nodes to the final buffer

            edge case:
                
                how do you handle the corners/walls of the bounding box for the entire octree?
                - There wouldn't be any neighbors to one side
            
            */


            return;  // ending the function
        }


        // grabs the leaf cell for a given position
        public: int GetLeafIndex (double posX_, double posY_, double posZ_)
        {
            // correcting the position so it fits within the octree (otherwise it won't return a valid position and will throw errors)
            double posX = std::min <double> (std::max <double> (posX_, offsetX), sizeX + offsetX);
            double posY = std::min <double> (std::max <double> (posY_, offsetY), sizeY + offsetY);
            double posZ = std::min <double> (std::max <double> (posZ_, offsetZ), sizeZ + offsetZ);

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
                lastLeafDepth = i;  // tracking the depth
                
                // getting the cell dimensions
                widthScalar = depthSizeScalars(i);
                cellSizeX = sizeX * widthScalar;
                cellSizeY = sizeY * widthScalar;
                cellSizeZ = sizeZ * widthScalar;

                // getting the offset for the child
                childOffsetX = (int) ((posX - baseX) / cellSizeX);
                childOffsetY = (int) ((posY - baseY) / cellSizeY);
                childOffsetZ = (int) ((posZ - baseZ) / cellSizeZ);

                childIndex = offsetOrderChildIndex(childOffsetX, childOffsetY, childOffsetZ);  // the child's index
                
                lastNodeIndex = nodeIndex;
                nodeIndex = childPointReferences(nodeIndex, childIndex);  // the new node index
                if (nodeIndex < 0) {  // making sure the node won't go off of an empty cell
                    nodeIndex = lastNodeIndex;
                    break;
                }

                // adjusting the base position for the corner of the cell
                baseX += cellSizeX * (double) childOffsetX;
                baseY += cellSizeY * (double) childOffsetY;
                baseZ += cellSizeZ * (double) childOffsetZ;

                // checking if the search has concluded
                if (numPositionIndexs(nodeIndex)) break;  // checking if a value has been added -- only leaf nodes contain values

                // adding the node to the path
                depthIndexBufferSearch(i) = nodeIndex;

            }

            // adding the node to the path
            depthIndexBufferSearch(lastLeafDepth) = nodeIndex;

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
            double distanceTraveledX = 0.0;
            double distanceTraveledY = 0.0;
            double distanceTraveledZ = 0.0;
            for (int i = 0; i < lastLeafDepth; i++)
            {
                distanceTraveledX += sizeX * depthSizeScalars(currentDepth - 1);  // tracking the distance traveled between cells (so an accuret extended search can be performed)
                distanceTraveledY += sizeY * depthSizeScalars(currentDepth - 1);
                distanceTraveledZ += sizeZ * depthSizeScalars(currentDepth - 1);

                if (distance < 99999999.0) break;  // if any point is found the distance is returned

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


        // gets a leaf node
        public: int GetLeafNodeFromLeafIndex (int leafIndex) {return leafNodeIndexes(leafIndex);}

        // gets the number of leaf nodes
        public: int GetNumberOfLeafNodes () {return numberOfLeafNodes;}

        // gets the position of a given node in global space (from a node index)
        public: CArray <double> GetLeafNodePosition (int leafNodeIndex)
        {
            CArray <double> outputPosition = CArray <double> (3);  // stores the final position that's calculated

            // doing a reverse lookup to find the position  (based on the inherent method used to construct the indexes, they should be in assending order; thus a binary search is being used)
            int leafIndex = BinarySearch <int> (leafNodeIndexes, numberOfLeafNodes, leafNodeIndex);

            // getting the position
            outputPosition(0) = leafNodePositions(leafIndex, 0);
            outputPosition(1) = leafNodePositions(leafIndex, 1);
            outputPosition(2) = leafNodePositions(leafIndex, 2);

            return outputPosition;  // returning the position
        }

        // gets the leaf index from a general index
        public: int GetLeafIndexFromIndex (int leafNodeIndex) {return BinarySearch <int> (leafNodeIndexes, numberOfLeafNodes, leafNodeIndex);}

        // returns the size of the octree
        public: int GetSize () {return numberChildReferences;}


        // gets the point at a given node index
        public: double GetNodePointX (int nodeIndex, int bufferIndex) {return points(positionIndexesPlusOne(nodeIndex, bufferIndex) - 1, 0);}
        public: double GetNodePointY (int nodeIndex, int bufferIndex) {return points(positionIndexesPlusOne(nodeIndex, bufferIndex) - 1, 1);}
        public: double GetNodePointZ (int nodeIndex, int bufferIndex) {return points(positionIndexesPlusOne(nodeIndex, bufferIndex) - 1, 2);}

    };
}

