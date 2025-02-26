#include <iostream>
#include "matar.h"

// the maximum depth a binary search can go to
const int MAX_BINARY_SEARCH_ITTERATIONS = 250;



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
        int halfWidth, currentIndexOld;
        for (int i = 0; i < MAX_BINARY_SEARCH_ITTERATIONS; i++)
        {
            // dividing the space in two
            halfWidth = (int) (dividedSize * 0.5);
            dividedSize -= halfWidth;

            T middleValue = points(currentIndex + halfWidth);

            currentIndexOld = currentIndex;
            // checking if the value is smaller or greater than the input value
            if (middleValue == searchValue) return currentIndex + halfWidth;  // concluding the search once the value is found
            if (middleValue < searchValue) currentIndex += halfWidth;  // checking which half the value is in the move the index; allowing for odd and even sized arrays
        }

        std::cout << "Binary Search Failed--Reached Max Depth: Octree.hpp; Octree::BinarySearch\n    Try increasing max search itterations ( MAX_BINARY_SEARCH_ITTERATIONS ) at top of header file" << std::endl;
        return -1;  // returning the last index (the search ran out of search depth likely)
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
            NodeReferenceCacheRecursiveAscent(1, rootNodeReferenceIndex, nodeRefferenceIndex);

            return;  // ending the function
        }

        // recursively jumps down all branches and starts the generation of the leaf nodes' neighbors
        private: void NodeReferenceCacheRecursiveAscent (const int depth, const int childIndex, int &nodeRefferenceIndex)
        {
            // going through all children for this node
            for (int i = 0; i < 8; i++)
            {
                // checking if the child is a leaf node
                if (childPointReferences(childIndex, i) == -1) {
                    GenerateLeafNeighbors(depth, childIndex, nodeRefferenceIndex);  // generating the cells neighbors
                    return;  // the leaf node had all neighbors generated, no need to continue
                }

                // continuning the recursive ascent
                NodeReferenceCacheRecursiveAscent(depth + 1, childPointReferences(childIndex, i), nodeRefferenceIndex);
            }

            return;  // ending the function
        }

        // generates the neighbors for a given leaf node
        private: void GenerateLeafNeighbors (const int depth, const int childIndex, int &nodeRefferenceIndex)
        {
            // getting the size of the node
            double scaledDepth = depthSizeScalars(depth);
            double nodeSizeX = sizeX * scaledDepth;
            double nodeSizeY = sizeY * scaledDepth;
            double nodeSizeZ = sizeZ * scaledDepth;

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

            // stores the offset values
            CArray <int> shiftBuffer = CArray <int> (depth + 1, 3);

            CArray <double> trueNodePosition = GetLeafNodePosition(childIndex);  // probably a lot slower, but I couldn't figure out how to get the position tracking to work for some reason and the whole generation algerithm is running fast enough considering it's a one time operation
            double searchPosX = trueNodePosition(0);
            double searchPosY = trueNodePosition(1);
            double searchPosZ = trueNodePosition(2);

            // looping until the inital leaf-node has been found again (to generate the path)
            for (int i = 1; i < depth; i++)  // by using the current depth, no if checks are needed to check for a leaf node sense it's already known this is a leaf node
            {
                // finding the size of the node
                widthScalar = depthSizeScalars(i);
                cellSizeX = sizeX * widthScalar;
                cellSizeY = sizeY * widthScalar;
                cellSizeZ = sizeZ * widthScalar;

                // getting the current offset
                childOffsetX = static_cast <int> ((searchPosX - baseX) / cellSizeX);
                childOffsetY = static_cast <int> ((searchPosY - baseY) / cellSizeY);
                childOffsetZ = static_cast <int> ((searchPosZ - baseZ) / cellSizeZ);

                // cashing the shift
                shiftBuffer(i, 0) = childOffsetX;
                shiftBuffer(i, 1) = childOffsetY;
                shiftBuffer(i, 2) = childOffsetZ;

                ascentChildIndex = offsetOrderChildIndex(childOffsetX, childOffsetY, childOffsetZ);  // the index of the child within a given parent node
                ascentNodeIndex = childPointReferences(ascentNodeIndex, ascentChildIndex);

                // adjusting the base position
                baseX += cellSizeX * static_cast <double> (childOffsetX);
                baseY += cellSizeY * static_cast <double> (childOffsetY);
                baseZ += cellSizeZ * static_cast <double> (childOffsetZ);

                // adding the node to the path
                depthIndexBufferSearch(i) = ascentNodeIndex;
            }

            depthIndexBufferSearch(depth) = childIndex;

            // getting the current offset to find the opposing neighbors
            childOffsetX = static_cast <int> ((searchPosX - baseX) / nodeSizeX);
            childOffsetY = static_cast <int> ((searchPosY - baseY) / nodeSizeY);
            childOffsetZ = static_cast <int> ((searchPosZ - baseZ) / nodeSizeZ);
            
            // cashing the shift
            shiftBuffer(depth, 0) = childOffsetX;
            shiftBuffer(depth, 1) = childOffsetY;
            shiftBuffer(depth, 2) = childOffsetZ;

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

            // getting all nodes on the oposing faces of other parent nodes (not under the main parent node that harbors the first 3)
            FindChildFace(depth, childIndex, nodeRefferenceIndex, shiftBuffer, 0);
            FindChildFace(depth, childIndex, nodeRefferenceIndex, shiftBuffer, 1);
            FindChildFace(depth, childIndex, nodeRefferenceIndex, shiftBuffer, 2);

            return;  // ending the function
        }

        // finds a given oposing face based on the axis; x, y, or z)
        private: void FindChildFace (const int depth, const int childIndex, int &nodeRefferenceIndex, const CArray <int> &shiftBuffer, const int axis)
        {
            // going up the tree until the correct shift value has been found
            int lookingForShift = 1 - shiftBuffer(depth, axis);

            // going up the tree until the shift is found
            int newChildIndex;
            for (int i = depth - 1; i > 0; i--)
            {
                // checking if it's the correct shift
                if (shiftBuffer(i, axis) == lookingForShift) {
                    // finding the offset order
                    int offsetIndex, indexI, indexJ;
                    if (!axis) {
                        offsetIndex = offsetOrderChildIndex(shiftBuffer(depth, axis), 0, 0);
                        indexI = 1; indexJ = 2;
                    }
                    else if (axis == 1) {
                        offsetIndex = offsetOrderChildIndex(0, shiftBuffer(depth, axis), 0);  // finding the correct order based on the current axis
                        indexI = 0; indexJ = 2;
                    }
                    else {
                        offsetIndex = offsetOrderChildIndex(0, 0, shiftBuffer(depth, axis));
                        indexI = 0; indexJ = 1;
                    }

                    // getting the face nodes
                    newChildIndex = childPointReferences(depthIndexBufferSearch(i), offsetIndex);
                    if (newChildIndex == -1) return;  // ending early, there is no valid face on this side?
                    GetChildFaceNodesRecursive(
                            depth,
                            newChildIndex,
                            nodeRefferenceIndex, shiftBuffer, axis,
                            indexI, indexJ
                        );

                    return;  // ending the search sense the face was found
                }
            }

            return;  // ending the function, no valid faces opose this side (likely the node was on the very very edge)
        }

        // gets all the nodes on a given face
        private: void GetChildFaceNodesRecursive (const int depth, const int childIndex, int &nodeRefferenceIndex, const CArray <int> &shiftBuffer, const int axis, const int indexI, const int indexJ)
        {
            // finding the offsets
            CArray <int> offsetKey = CArray <int> (3);  // to manage the offset orders
            offsetKey(axis) = shiftBuffer(depth, axis);  // adding the initial shift

            // going through i and j and continuing the search for each node
            int newChildIndex;
            for (int i = 0; i <= 1; i++) {
                for (int j = 0; j <= 1; j++) {
                    // constructing the index
                    offsetKey(indexI) = i;
                    offsetKey(indexJ) = j;

                    // checking if the child node is empty (in which case this is a leaf node and it'll be added to the traversal cashe)
                    newChildIndex = childPointReferences(childIndex, offsetOrderChildIndex(offsetKey(0), offsetKey(1), offsetKey(2)));
                    if (newChildIndex == -1)
                    {
                        // getting the face node index and adding it to the traversal table buffer
                        cellNeighborReferences(depthIndexBufferSearch(depth), cellNeighborBufferSize(depthIndexBufferSearch(depth))) = childPointReferences(
                            childIndex,
                            offsetOrderChildIndex(offsetKey(0), offsetKey(1), offsetKey(2))
                        );
                        cellNeighborBufferSize(depthIndexBufferSearch(depth))++;
                    } else {  // continuing the search if the node wasn't a leaf
                        // continuing the search
                        GetChildFaceNodesRecursive(
                            depth,
                            newChildIndex,
                            nodeRefferenceIndex, shiftBuffer, axis,
                            indexI, indexJ
                        );
                    }
                }
            }
            
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
                childOffsetX = static_cast <int> ((posX - baseX) / cellSizeX);
                childOffsetY = static_cast <int> ((posY - baseY) / cellSizeY);
                childOffsetZ = static_cast <int> ((posZ - baseZ) / cellSizeZ);

                childIndex = offsetOrderChildIndex(childOffsetX, childOffsetY, childOffsetZ);  // the child's index
                
                lastNodeIndex = nodeIndex;
                nodeIndex = childPointReferences(nodeIndex, childIndex);  // the new node index
                if (nodeIndex < 0) {  // making sure the node won't go off of an empty cell
                    nodeIndex = lastNodeIndex;
                    break;
                }

                // adjusting the base position for the corner of the cell
                baseX += cellSizeX * static_cast <double> (childOffsetX);
                baseY += cellSizeY * static_cast <double> (childOffsetY);
                baseZ += cellSizeZ * static_cast <double> (childOffsetZ);

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
            // getting the leaf node to begin the search (and get the path to back-trace efficently)
            int leafNodeIndex = GetLeafIndex(samplePositionX, samplePositionY, samplePositionZ);

            return NearestNeighborSearch(leafNodeIndex, samplePositionX, samplePositionY, samplePositionZ);  // returning the distance
        }

        // preforms a nearest neighbor search based on a given point
        public: double NearestNeighborSearch (int leafNodeIndex, double samplePositionX, double samplePositionY, double samplePositionZ)
        {
            double distance = std::numeric_limits<double>::max();  // the minimum distance found
            
            // checking the initial point
            BranchSearch(depthIndexBufferSearch(lastLeafDepth), -1, distance, samplePositionX, samplePositionY, samplePositionZ);
            
            // itterating over the neighboring cells
            int currentDepth = lastLeafDepth;
            for (int i = 0; i < lastLeafDepth; i++)
            {
                if (distance < std::numeric_limits<double>::max()) break;  // if any point is found the distance is returned

                // starting a new search through the whole branch
                currentDepth--;
                BranchSearch(depthIndexBufferSearch(currentDepth), depthIndexBufferSearch(currentDepth + 1), distance, samplePositionX, samplePositionY, samplePositionZ);
            }

            // returning the distance
            return sqrt(distance);
        }

        // searching down the entirety of a branch      depth-first approach (goes fully down the first node, then expands to the surrounding nodes and goes to the max depth, and so on)
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

