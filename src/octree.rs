#![allow(non_snake_case)]

// some of the methods aren't being used,
// but are implimented as they may be usefull
// in future projects
#![allow(dead_code)]

use crate::priorityQueue;

const MAX_BINARY_SEATCH_ITTERATIONS: isize = 250;


pub fn BinarySearch (points: &Vec <usize>, searchValue: &usize) -> Option <usize> {
    let mut currentIndex: usize = 0;
    let mut dividedSize = points.len();

    let mut halfWidth: usize;
    // edit this?
    for _i in 0..MAX_BINARY_SEATCH_ITTERATIONS {
        halfWidth = dividedSize / 2;
        dividedSize -= halfWidth;

        if let Some(middleValue) = points.get(currentIndex + halfWidth) {
            if middleValue == searchValue {
                return Some(currentIndex + halfWidth);
            }
            if middleValue < searchValue {
                currentIndex += halfWidth;
            }
        }
    }
    None
}


const CHILD_OFFSET_ORDER: [(usize, usize, usize); 8] = [
    (1, 0, 0),
    (1, 1, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 0, 1),
    (1, 1, 1),
    (0, 1, 1),
    (0, 0, 0)
];

lazy_static::lazy_static! {
    static ref CHILD_OFFSET_ORDER_INDEX:
        std::collections::HashMap <[usize; 3], usize> = {
            std::collections::HashMap::from([
                ([1, 0, 0], 0),
                ([1, 1, 0], 1),
                ([0, 1, 0], 2),
                ([0, 0, 1], 3),
                ([1, 0, 1], 4),
                ([1, 1, 1], 5),
                ([0, 1, 1], 6),
                ([0, 0, 0], 7),
            ])
        };
}


pub struct Octree {
    childPointReferences: Vec <[Option <usize>; 8]>,
    positionIndexesPlusOne: Vec <Vec <usize>>,
    leafNodePositions: Vec <Option <(f64, f64, f64)>>,
    leafNodeDepths: Vec <Option <usize>>,
    numberOfLeafNodes: usize,
    offsetX: f64,
    offsetY: f64,
    offsetZ: f64,
    sizeX: f64,
    sizeY: f64,
    sizeZ: f64,
    rootNodeReferenceIndex: usize,
    maxDepth: usize,
    depthSizeScalars: Vec <f64>,
    depthIndexBufferSearch: Vec <usize>,
    lastLeafDepth: usize,
    cellNeighborRefferences: Vec <Vec <usize>>,
    leafNodes: Vec <usize>,

    // corner references...    using morton codes floor((pos - minTreeBoundingPos) / smallestCellSize)
    cornerPoints: std::collections::HashMap <u128, (f64, f64, f64)>,
    
    // the indexes in the corner points hash-map
    // a hash-map is needed to allow for checking
    // for duplicate points
    nodeCornerReferences: Vec <Vec <u128>>,
}


impl Octree {
    pub fn new ((xOffset, yOffset, zOffset): (f64, f64, f64),
                (xSize, ySize, zSize): (f64, f64, f64),
                maximumDepth: usize) -> Self {
        Octree {
            childPointReferences: vec!(),
            positionIndexesPlusOne: vec!(),
            leafNodePositions: vec!(),
            leafNodeDepths: vec!(),
            numberOfLeafNodes: 0,
            offsetX: xOffset,
            offsetY: yOffset,
            offsetZ: zOffset,
            sizeX: xSize,
            sizeY: ySize,
            sizeZ: zSize,
            rootNodeReferenceIndex: 0,
            maxDepth: maximumDepth,
            depthSizeScalars: {
                let mut depthScalars: Vec <f64> = vec!();
                for depth in 0..=maximumDepth {
                    depthScalars.push(
                        1.0 / 2.0f64.powf(depth as f64)
                    );
                }
                depthScalars
            },
            depthIndexBufferSearch: vec!(),
            lastLeafDepth: 0,
            cellNeighborRefferences: vec!(),
            leafNodes: vec!(),
            cornerPoints: std::collections::HashMap::new(),
            nodeCornerReferences: vec!(),  // one for every cell (representing it's 8 corners)
        }
    }


    pub fn GetCornerMortonCode (&self, posX: f64, posY: f64, posZ: f64) -> u128 {
        // getting the intager position for the morton code
        // floor((pos - minTreeBoundingPos) / smallestCellSize)
        let x = ((posX - self.offsetX) / self.depthSizeScalars[self.maxDepth]) as u32;
        let y = ((posY - self.offsetY) / self.depthSizeScalars[self.maxDepth]) as u32;
        let z = ((posZ - self.offsetZ) / self.depthSizeScalars[self.maxDepth]) as u32;

        // getting the morton code
        let mut mortonCode: u128 = 0;
        for i in 0..32 {  // 64-bit integer
            mortonCode |= (((x >> i) & 1) as u128) << (3 * i);
            mortonCode |= (((y >> i) & 1) as u128) << (3 * i + 1);
            mortonCode |= (((z >> i) & 1) as u128) << (3 * i + 2);
        } mortonCode
    }

    pub fn GenerateCornerPointsChache (&mut self) -> usize{
        let mut totalCornerPoints = 0usize;

        let mut mortonCode: u128;
        let mut depthScalar: f64;
        let mut cellSize: (f64, f64, f64);
        let mut position: (f64, f64, f64);
        let mut cornerOffset: (usize, usize, usize);
        // going through every leaf node
        // leaf index is being used; use a binary search to go
        // from a node index back to leaf index (not needed inside here)
        for (leafIndex, nodeIndex) in self.leafNodes.iter().enumerate() {
            position = self.GetLeafPosition(*nodeIndex).
                    expect("Failed to get leaf node position...");
            depthScalar = self.depthSizeScalars[
                self.leafNodeDepths[*nodeIndex].
                    expect("Failed to get node depth...")
            ];
            cellSize = (
                self.sizeX * depthScalar,
                self.sizeY * depthScalar,
                self.sizeZ * depthScalar,
            );

            self.nodeCornerReferences.push(vec!());
            
            // getting the 8 corners
            for i in 0..8 {
                cornerOffset = CHILD_OFFSET_ORDER[i];

                // getting the morton code
                mortonCode = self.GetCornerMortonCode(
                    position.0 + cellSize.0 * cornerOffset.0 as f64,
                    position.1 + cellSize.1 * cornerOffset.1 as f64,
                    position.2 + cellSize.2 * cornerOffset.2 as f64,
                );
                if self.cornerPoints.insert(mortonCode, position).is_none() {
                    totalCornerPoints += 1;
                }
                self.nodeCornerReferences[leafIndex].push(mortonCode);
            }
        } totalCornerPoints
    }


    pub fn SubDivide (&mut self, points: &Vec <(f64, f64, f64)>) {
        self.numberOfLeafNodes = 0;

        self.rootNodeReferenceIndex = self.SubDivideRecursive (
            points,
            (0, 0, 0),
            (self.offsetX, self.offsetY, self.offsetZ),
            1
        );

        // I'll leave this for now
        println!("Num cells: {}      root node: {}", self.childPointReferences.len(), self.rootNodeReferenceIndex);
    }

    fn SubDivideRecursive ( &mut self,
                            points: &Vec <(f64, f64, f64)>,
                            shift: (usize, usize, usize),
                            position: (f64, f64, f64),
                            depth: usize) -> usize {
        
        let scaledDepth = self.depthSizeScalars[depth - 1];
        let nodeSize: (f64, f64, f64) = (
            self.sizeX * scaledDepth,
            self.sizeY * scaledDepth,
            self.sizeZ * scaledDepth
        );
        let newPosition: (f64, f64, f64) = (
            position.0 + nodeSize.0 * shift.0 as f64,
            position.1 + nodeSize.1 * shift.1 as f64,
            position.2 + nodeSize.2 * shift.2 as f64
        );

        if depth == self.maxDepth {
            self.leafNodes.push(self.childPointReferences.len());
            self.leafNodePositions.push(Some((
                newPosition.0 + nodeSize.0 * 0.5,
                newPosition.1 + nodeSize.1 * 0.5,
                newPosition.2 + nodeSize.2 * 0.5
            )));
            self.leafNodeDepths.push(Some(depth));
            self.childPointReferences.push([ None; 8 ]);

            self.positionIndexesPlusOne.push(vec!());
            for (i, point) in points.iter().enumerate() {
                if  point.0 >= newPosition.0 && point.0 < newPosition.0 + nodeSize.0 &&
                    point.1 >= newPosition.1 && point.1 < newPosition.1 + nodeSize.1 &&
                    point.2 >= newPosition.2 && point.2 < newPosition.2 + nodeSize.2 {

                    if let Some(buffer) = self.positionIndexesPlusOne.last_mut() {
                        buffer.push(i + 1);
                    }
                }
            }

            // returning the node index
            self.cellNeighborRefferences.push(vec!());
            return self.childPointReferences.len() - 1;
        }

        let mut lastIndex: Option <usize> = None;
        let mut numberBoundingPoints: usize = 0;
        for (i, point) in points.iter().enumerate() {
            if  point.0 >= newPosition.0 && point.0 < newPosition.0 + nodeSize.0 &&
                point.1 >= newPosition.1 && point.1 < newPosition.1 + nodeSize.1 &&
                point.2 >= newPosition.2 && point.2 < newPosition.2 + nodeSize.2 {
            
                lastIndex = Some(i);
                numberBoundingPoints += 1;
                if numberBoundingPoints > 1 {  break;  }
            }
        }

        if numberBoundingPoints == 0 {
            self.leafNodes.push(self.childPointReferences.len());
            self.leafNodePositions.push(Some((
                newPosition.0 + nodeSize.0 * 0.5,
                newPosition.1 + nodeSize.1 * 0.5,
                newPosition.2 + nodeSize.2 * 0.5
            )));
            self.leafNodeDepths.push(Some(depth));
            self.childPointReferences.push([ None; 8 ]);

            self.positionIndexesPlusOne.push(vec!());
            self.cellNeighborRefferences.push(vec!());

            return self.childPointReferences.len() - 1;
        }

        if numberBoundingPoints == 1 {
            self.leafNodes.push(self.childPointReferences.len());
            self.leafNodePositions.push(Some((
                newPosition.0 + nodeSize.0 * 0.5,
                newPosition.1 + nodeSize.1 * 0.5,
                newPosition.2 + nodeSize.2 * 0.5
            )));
            self.leafNodeDepths.push(Some(depth));
            self.childPointReferences.push([ None; 8 ]);

            self.positionIndexesPlusOne.push(vec![
                lastIndex.expect("Failed to find point")
            ]);
            self.cellNeighborRefferences.push(vec!());

            return self.childPointReferences.len() - 1;
        }

        let newChildIndexes: [Option <usize>; 8] = [
            Some(self.SubDivideRecursive(points, CHILD_OFFSET_ORDER[0], newPosition, depth + 1)),
            Some(self.SubDivideRecursive(points, CHILD_OFFSET_ORDER[1], newPosition, depth + 1)),
            Some(self.SubDivideRecursive(points, CHILD_OFFSET_ORDER[2], newPosition, depth + 1)),
            Some(self.SubDivideRecursive(points, CHILD_OFFSET_ORDER[3], newPosition, depth + 1)),
            Some(self.SubDivideRecursive(points, CHILD_OFFSET_ORDER[4], newPosition, depth + 1)),
            Some(self.SubDivideRecursive(points, CHILD_OFFSET_ORDER[5], newPosition, depth + 1)),
            Some(self.SubDivideRecursive(points, CHILD_OFFSET_ORDER[6], newPosition, depth + 1)),
            Some(self.SubDivideRecursive(points, CHILD_OFFSET_ORDER[7], newPosition, depth + 1))
        ];

        self.childPointReferences.push(newChildIndexes);
        self.positionIndexesPlusOne.push(vec!());
        self.leafNodeDepths.push(None);
        self.leafNodePositions.push(None);
        self.cellNeighborRefferences.push(vec!());

        self.childPointReferences.len() - 1
    }

    
    pub fn GenerateNodeReferenceChache (&mut self) {
        self.NodeReferenceCacheRecursiveAscent(0, self.rootNodeReferenceIndex);
    }

    fn NodeReferenceCacheRecursiveAscent (&mut self, depth: usize, childIndex: usize) {
        for i in 0..8 {
            if let Some(child) = self.childPointReferences[childIndex][i] {
                self.NodeReferenceCacheRecursiveAscent(
                    depth + 1,
                    child
                );
            } else {
                self.GenerateLeafNeighbors(depth, childIndex);
                return;
            }
        }
    }

    fn GenerateLeafNeighbors (&mut self, depth: usize, childIndex: usize) {
        let scaledDepth = self.depthSizeScalars[depth];
        let nodeSize = (
            self.sizeX * scaledDepth,
            self.sizeY * scaledDepth,
            self.sizeZ * scaledDepth,
        );

        let mut ascentNodeIndex = self.rootNodeReferenceIndex;

        let mut base = (
            self.offsetX, self.offsetY, self.offsetZ
        );

        let mut childOffset: [usize; 3];
        let mut widthScalar: f64;
        let mut cellSize: (f64, f64, f64);
        let mut ascentChildIndex: usize;

        self.depthIndexBufferSearch.clear();
        self.depthIndexBufferSearch.push(ascentNodeIndex);

        let mut shiftBuffer: Vec <[usize; 3]> = vec!();
        shiftBuffer.push([0usize; 3]);

        let searchPosition =
            self.GetLeafPosition(childIndex).
                expect("Failed to get leaf node position");
        
        for i in 1..depth {
            widthScalar = self.depthSizeScalars[i];
            cellSize = (
                self.sizeX * widthScalar,
                self.sizeY * widthScalar,
                self.sizeZ * widthScalar,
            );

            childOffset = [
                ((searchPosition.0 - base.0) / cellSize.0) as usize,
                ((searchPosition.1 - base.1) / cellSize.1) as usize,
                ((searchPosition.2 - base.2) / cellSize.2) as usize,
            ];

            shiftBuffer.push(childOffset);

            ascentChildIndex = *CHILD_OFFSET_ORDER_INDEX.get(&childOffset)
                               .expect("Invalid child offset");
            ascentNodeIndex = self.childPointReferences[ascentNodeIndex][ascentChildIndex]
                              .expect("Failed to gather child point reference");
            
            base = (
                base.0 + cellSize.0 * childOffset[0] as f64,
                base.1 + cellSize.1 * childOffset[1] as f64,
                base.2 + cellSize.2 * childOffset[2] as f64,
            );

            self.depthIndexBufferSearch.push(ascentNodeIndex);
        }

        self.depthIndexBufferSearch.push(childIndex);

        childOffset = [
            ((searchPosition.0 - base.0) / nodeSize.0) as usize,
            ((searchPosition.1 - base.1) / nodeSize.1) as usize,
            ((searchPosition.2 - base.2) / nodeSize.2) as usize,
        ];
        shiftBuffer.push(childOffset);

        self.GetChildFaceNodesRecursive(
            depth,
            self.childPointReferences[self.depthIndexBufferSearch[depth - 1]][
                    *CHILD_OFFSET_ORDER_INDEX.get(
                        &[1 - childOffset[0], childOffset[1], childOffset[2]]
                    ).expect("Invalid offset")
                ].expect("Failed to get child"),
            1 - childOffset[0],
            0,
            1,
            2
        );
        self.GetChildFaceNodesRecursive(
            depth,
            self.childPointReferences[self.depthIndexBufferSearch[depth - 1]][
                    *CHILD_OFFSET_ORDER_INDEX.get(
                        &[childOffset[0], 1 - childOffset[1], childOffset[2]]
                    ).expect("Invalid offset")
                ].expect("Failed to get child"),
            1 - childOffset[1],
            1,
            0,
            2
        );
        self.GetChildFaceNodesRecursive(
            depth,
            self.childPointReferences[self.depthIndexBufferSearch[depth - 1]][
                    *CHILD_OFFSET_ORDER_INDEX.get(
                        &[childOffset[0], childOffset[1], 1 - childOffset[2]]
                    ).expect("Invalid offset")
                ].expect("Failed to get child"),
            1 - childOffset[2],
            2,
            0,
            1
        );

        self.FindChildFace(depth, &shiftBuffer, 0);
        self.FindChildFace(depth, &shiftBuffer, 1);
        self.FindChildFace(depth, &shiftBuffer, 2);
    }

    fn FindChildFace (&mut self, depth: usize, shiftBuffer: &Vec <[usize; 3]>, axis: usize) {
        let lookingForShift = 1 - shiftBuffer[depth][axis];

        for i in (0..depth).rev() {
            if shiftBuffer[i][axis] == lookingForShift {
                let offsetIndex: usize;
                let indexI: usize;
                let indexJ: usize;

                match axis {
                    0 => {
                        offsetIndex = *CHILD_OFFSET_ORDER_INDEX.get(&[
                            shiftBuffer[depth][axis], 0, 0
                        ]).expect("Invalid shift");
                        indexI = 1; indexJ = 2;
                    },
                    1 => {
                        offsetIndex = *CHILD_OFFSET_ORDER_INDEX.get(&[
                            0, shiftBuffer[depth][axis], 0
                        ]).expect("Invalid shift");
                        indexI = 0; indexJ = 2;
                    },
                    _ => {
                        offsetIndex = *CHILD_OFFSET_ORDER_INDEX.get(&[
                            0, 0, shiftBuffer[depth][axis]
                        ]).expect("Invalid shift");
                        indexI = 0; indexJ = 1;
                    },
                }

                let newChildIndex =
                    self.childPointReferences[self.depthIndexBufferSearch[i]]
                    [offsetIndex].expect("Failed to get child");
                
                self.GetChildFaceNodesRecursive (
                    depth,
                    newChildIndex,
                    shiftBuffer[depth][axis],
                    axis,
                    indexI,
                    indexJ,
                );
                
                return;
            }
        }
    }

    fn GetChildFaceNodesRecursive (
        &mut self, depth: usize, childIndex: usize,
        shift: usize, axis: usize,
        indexI: usize, indexJ: usize) {

        let mut offsetKey = [0usize, 0usize, 0usize];
        offsetKey[axis] = shift;

        let mut newChildIndex: Option <usize>;
        for i in 0..=1 {
            for j in 0..=1 {
                offsetKey[indexI] = i;
                offsetKey[indexJ] = j;

                newChildIndex = self.childPointReferences
                    [childIndex][*CHILD_OFFSET_ORDER_INDEX.get(
                        &offsetKey).expect("Invalid shift")];
                
                if newChildIndex.is_none() {
                    self.cellNeighborRefferences[self.depthIndexBufferSearch[depth]].
                            push(childIndex);
                    return;
                }

                self.GetChildFaceNodesRecursive(
                    depth,
                    newChildIndex.unwrap(),
                    shift,
                    axis,
                    indexI,
                    indexJ
                );
            }
        }
    }


    pub fn GetLeafIndex (&mut self, inputPos: (f64, f64, f64)) -> usize {
        let pos = (
            inputPos.0.max(self.offsetX).min(self.sizeX + self.offsetX),
            inputPos.1.max(self.offsetY).min(self.sizeY + self.offsetY),
            inputPos.2.max(self.offsetZ).min(self.sizeZ + self.offsetZ),
        );

        let mut nodeIndex = Some(self.rootNodeReferenceIndex);
        let mut lastNodeIndex: Option <usize>;

        let mut base = (
            self.offsetX, self.offsetY, self.offsetZ
        );

        let mut childOffset: (usize, usize, usize);
        let mut childIndex: usize;
        let mut widthScalar: f64;

        let mut nodeSize: (f64, f64, f64);

        self.depthIndexBufferSearch.clear();
        self.depthIndexBufferSearch.push(nodeIndex.unwrap());

        for i in 1..self.maxDepth {
            self.lastLeafDepth = i;

            widthScalar = self.depthSizeScalars[i];
            nodeSize = (
                self.sizeX * widthScalar,
                self.sizeY * widthScalar,
                self.sizeZ * widthScalar,
            );

            childOffset = (
                ((pos.0 - base.0) / nodeSize.0) as usize,
                ((pos.1 - base.1) / nodeSize.1) as usize,
                ((pos.2 - base.2) / nodeSize.2) as usize,
            );

            childIndex = *CHILD_OFFSET_ORDER_INDEX.get(&[childOffset.0, childOffset.1, childOffset.2]).
                expect("Invalid offset");
            
            lastNodeIndex = nodeIndex;
            nodeIndex = self.childPointReferences[nodeIndex.unwrap()][childIndex];
            if nodeIndex.is_none() {
                nodeIndex = lastNodeIndex;
                break;
            }

            base = (
                base.0 + nodeSize.0 * childOffset.0 as f64,
                base.1 + nodeSize.1 * childOffset.1 as f64,
                base.2 + nodeSize.2 * childOffset.2 as f64,
            );

            // do i actually need to check if this position has any points?
            // wouldn't the previous if statement cover that?

            self.depthIndexBufferSearch.push(nodeIndex.unwrap());
        }

        self.depthIndexBufferSearch.push(nodeIndex.unwrap());

        nodeIndex.unwrap()
    }


    fn GetSquaredDistance (&self, point1: (f64, f64, f64), point2: (f64, f64, f64)) -> f64 {
        let deltaPos = (
            point1.0 - point2.0,
            point1.1 - point2.1,
            point1.2 - point2.2,
        );

        deltaPos.0*deltaPos.0 + deltaPos.1*deltaPos.1 + deltaPos.2*deltaPos.2
    }

    pub fn NearestNeighborSearch (&self,
                                  points: &Vec <(f64, f64, f64)>,
                                  leafNodeIndex: usize,
                                  samplePosition: (f64, f64, f64)
            ) -> Option <f64> {
        
        let mut queue = priorityQueue::MinHeapBinaryTree::new();

        let queryPoint = (
            samplePosition.0, samplePosition.1, samplePosition.2
        );

        queue.Push(
            (0.0, leafNodeIndex)
        );

        let mut distance: f64;
        let mut nodeWidth: f64;
        let mut minNodeDst: f64;
        let mut depthScalar: f64;
        let mut currentNode: (f64, usize);
        let mut nodeWidthXYZ: (f64, f64, f64);

        let mut minDst = f64::MAX;

        let mut visited: std::collections::HashMap <usize, bool> =
            std::collections::HashMap::new();
        visited.insert(leafNodeIndex, true);

        while !queue.IsEmpty() {
            // the queue has to not be empty to get here
            currentNode = queue.Pop().unwrap();

            depthScalar = self.depthSizeScalars[
                self.leafNodeDepths[currentNode.1]
                .expect("Not a leaf node")];
            nodeWidthXYZ = (
                self.sizeX * depthScalar,
                self.sizeY * depthScalar,
                self.sizeZ * depthScalar
            );
            nodeWidth = nodeWidthXYZ.0*nodeWidthXYZ.0 +
                        nodeWidthXYZ.1*nodeWidthXYZ.1 +
                        nodeWidthXYZ.2*nodeWidthXYZ.2;
            
            minNodeDst = currentNode.0 + nodeWidth * 0.25 -
                         (currentNode.0 * nodeWidth).sqrt();
            
            if minNodeDst >= minDst {
                return Some(minDst.sqrt())
            }

            visited.insert(currentNode.1, true);

            for positionIndex in &self.positionIndexesPlusOne[currentNode.1] {
                let point = points[*positionIndex];
                distance = self.GetSquaredDistance(point, samplePosition);
                minDst = minDst.min(distance);
            }

            for neighborIndex in &self.cellNeighborRefferences[currentNode.1] {
                if visited.contains_key(neighborIndex) {
                    continue;
                }

                distance = self.GetSquaredDistance(
                    self.GetLeafPosition(*neighborIndex).
                                expect("Failed to get leaf position"),
                    queryPoint
                );

                queue.Push((
                    distance, *neighborIndex
                ));
            }
        }

        Some(minDst.sqrt())
        //None    ig it can actually end up searching every node.....
    }


    pub fn GetLeafNodes (&self) -> Vec <usize> {
        self.leafNodes.clone()
    }

    pub fn GetLeafPosition (&self, nodeIndex: usize) -> Option <(f64, f64, f64)> {
        if let Some(node) = self.leafNodePositions.get(nodeIndex) {
            if node.is_none() { return None; }
            return *node;
        } None
    }

}


