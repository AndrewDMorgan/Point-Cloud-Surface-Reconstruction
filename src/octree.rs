#![allow(non_snake_case)]

// some of the methods aren't being used,
// but are implemented as they may be useful
// in future projects or other inclusions
// of the octree data structure
#![allow(dead_code)]

use crate::priorityQueue;

const MAX_BINARY_SEARCH_ITERATIONS: isize = 10000;


pub fn BinarySearch <T: Eq + std::cmp::PartialOrd> (points: &Vec <T>, searchValue: &T) -> Option <usize> {
    let mut currentIndex: usize = 0;
    let mut dividedSize = points.len();

    let mut halfWidth: usize;
    // edit this?
    for _i in 0..MAX_BINARY_SEARCH_ITERATIONS {
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
    positionIndexes: Vec <Vec <usize>>,
    nodePositions: Vec <(f64, f64, f64)>,
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
    cellNeighborReferences: Vec <Vec <usize>>,
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
            positionIndexes: vec!(),
            nodePositions: vec!(),
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
            cellNeighborReferences: vec!(),
            leafNodes: vec!(),
            cornerPoints: std::collections::HashMap::new(),
            nodeCornerReferences: vec!(),  // one for every cell (representing its 8 corners)
        }
    }



    pub fn GetCornerPointReferences (&self, leafNode: usize) -> Vec <u128> {
        self.nodeCornerReferences[
            self.GetLeafNodeReverseIndex(leafNode)
        ].clone()
    }

    pub fn GetCornerPointKeys (&self) -> Vec <u128> {
        let mut keys: Vec <u128> = vec!();
        for corner in self.cornerPoints.keys() {
            keys.push(*corner);
        } keys
    }

    pub fn LookupCornerPoint (&self, id: u128) -> (f64, f64, f64) {
        *self.cornerPoints.get(&id).
            expect("Failed to gather corner point")
    }



    // TEMP remove once debugged
    pub fn GetAllNeighborCornersTemp(&self, leafNode: usize) -> Vec <(f64, f64, f64)> {
        let leafIndex = self.GetLeafNodeReverseIndex(leafNode);
        let mut points: Vec <(f64, f64, f64)> = vec!();

        for cornerCode in &self.nodeCornerReferences[leafIndex] {
            points.push(
                *self.cornerPoints.get(cornerCode).expect("Failed...")
            );
        }

        println!(": {}", self.cellNeighborReferences[leafNode].len());
        for neighborIndex in &self.cellNeighborReferences[leafNode] {
            for cornerCode in &self.nodeCornerReferences[
                self.GetLeafNodeReverseIndex(*neighborIndex)
            ] {
                points.push(
                    *self.cornerPoints.get(cornerCode).expect("Failed...")
                );
            }
        }

        points
    }





    fn GetCornerMortonCode (&self, posX: f64, posY: f64, posZ: f64) -> u128 {
        // getting the integer position for the morton code
        // floor((pos - minTreeBoundingPos) / smallestCellSize)
        // is the tiny offset needed to prevent float point errors?

        // 1, 0.5, 0.25, 0.125       1/2^x    1/1  1/2  1/4  1/8
        // any smaller unit should allow bigger units to fit in correctly so that should be fine
        // offset them by half a unit to ensure floating point errors don't add up....
        // that didn't help..........

        let x = ((posX - self.offsetX + self.sizeX*self.depthSizeScalars[self.maxDepth]*0.5) / self.depthSizeScalars[self.maxDepth]) as u32;
        let y = ((posY - self.offsetY + self.sizeX*self.depthSizeScalars[self.maxDepth]*0.5) / self.depthSizeScalars[self.maxDepth]) as u32;
        let z = ((posZ - self.offsetZ + self.sizeX*self.depthSizeScalars[self.maxDepth]*0.5) / self.depthSizeScalars[self.maxDepth]) as u32;

        self.GetMortonCode(x, y, z)
    }

    pub fn GenerateCornerPointsCache(&mut self) -> usize{
        let mut totalCornerPoints = 0usize;

        let mut mortonCode: u128;
        let mut depthScalar: f64;
        let mut cellSize: (f64, f64, f64);
        let mut position: (f64, f64, f64);
        let mut newPosition: (f64, f64, f64);
        let mut cornerOffset: (usize, usize, usize);
        // going through every leaf node
        // leaf index is being used; use a binary search to go
        // from a node index back to leaf index (not needed inside here)
        for (leafIndex, nodeIndex) in self.leafNodes.iter().enumerate() {
            position = self.GetNodePosition(*nodeIndex);
            depthScalar = self.depthSizeScalars[
                self.leafNodeDepths[*nodeIndex].
                    expect("Failed to get node depth...")
            ];
            cellSize = (
                self.sizeX * depthScalar * 0.5,
                self.sizeY * depthScalar * 0.5,
                self.sizeZ * depthScalar * 0.5,
            );

            self.nodeCornerReferences.push(vec!());
            
            // getting the 8 corners
            for i in 0..8 {
                cornerOffset = CHILD_OFFSET_ORDER[i];

                // getting the morton code
                newPosition = (
                    position.0 + cellSize.0 * (cornerOffset.0 as f64 * 2.0 - 1.0),
                    position.1 + cellSize.1 * (cornerOffset.1 as f64 * 2.0 - 1.0),
                    position.2 + cellSize.2 * (cornerOffset.2 as f64 * 2.0 - 1.0),
                );
                mortonCode = self.GetCornerMortonCode(newPosition.0, newPosition.1, newPosition.2);
                if self.cornerPoints.insert(mortonCode, newPosition).is_none() {
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
            0 // starting at 0 doesn't change anything it seems?
        );

        // I'll leave this for now
        println!("Num cells: {}      root node: {}", self.childPointReferences.len(), self.rootNodeReferenceIndex);
    }

    fn SubDivideRecursive ( &mut self,
                            points: &Vec <(f64, f64, f64)>,
                            shift: (usize, usize, usize),
                            position: (f64, f64, f64),
                            depth: usize) -> usize {
        
        let scaledDepth = self.depthSizeScalars[depth];
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
            self.nodePositions.push((
                newPosition.0 + nodeSize.0 * 0.5,
                newPosition.1 + nodeSize.1 * 0.5,
                newPosition.2 + nodeSize.2 * 0.5
            ));
            self.leafNodeDepths.push(Some(depth));
            self.childPointReferences.push([ None; 8 ]);

            self.positionIndexes.push(vec!());
            for (i, point) in points.iter().enumerate() {
                if  point.0 >= newPosition.0 && point.0 < newPosition.0 + nodeSize.0 &&
                    point.1 >= newPosition.1 && point.1 < newPosition.1 + nodeSize.1 &&
                    point.2 >= newPosition.2 && point.2 < newPosition.2 + nodeSize.2 {

                    if let Some(buffer) = self.positionIndexes.last_mut() {
                        buffer.push(i);
                    }
                }
            }

            // returning the node index
            self.cellNeighborReferences.push(vec!());
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
            self.nodePositions.push((
                newPosition.0 + nodeSize.0 * 0.5,
                newPosition.1 + nodeSize.1 * 0.5,
                newPosition.2 + nodeSize.2 * 0.5
            ));
            self.leafNodeDepths.push(Some(depth));
            self.childPointReferences.push([ None; 8 ]);

            self.positionIndexes.push(vec!());
            self.cellNeighborReferences.push(vec!());

            return self.childPointReferences.len() - 1;
        }

        if numberBoundingPoints == 1 {
            self.leafNodes.push(self.childPointReferences.len());
            self.nodePositions.push((
                newPosition.0 + nodeSize.0 * 0.5,
                newPosition.1 + nodeSize.1 * 0.5,
                newPosition.2 + nodeSize.2 * 0.5
            ));
            self.leafNodeDepths.push(Some(depth));
            self.childPointReferences.push([ None; 8 ]);

            self.positionIndexes.push(vec![
                lastIndex.expect("Failed to find point")
            ]);
            self.cellNeighborReferences.push(vec!());

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
        self.positionIndexes.push(vec!());
        self.leafNodeDepths.push(None);
        self.nodePositions.push((
            newPosition.0 + nodeSize.0 * 0.5,
            newPosition.1 + nodeSize.1 * 0.5,
            newPosition.2 + nodeSize.2 * 0.5
        ));
        self.cellNeighborReferences.push(vec!());

        self.childPointReferences.len() - 1
    }

    
    pub fn GenerateNodeReferenceCache(&mut self) {
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

        let mut childOffset: [usize; 3] = [0usize; 3];
        let mut widthScalar: f64;
        let mut cellSize: (f64, f64, f64);
        let mut ascentChildIndex: usize;

        self.depthIndexBufferSearch.clear();
        self.depthIndexBufferSearch.push(ascentNodeIndex);

        let mut shiftBuffer: Vec <[usize; 3]> = vec!();
        shiftBuffer.push([0usize; 3]);

        let searchPosition =
            self.GetNodePosition(childIndex);
        
        for i in 1..=depth {
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

            base = (
                base.0 + cellSize.0 * childOffset[0] as f64,
                base.1 + cellSize.1 * childOffset[1] as f64,
                base.2 + cellSize.2 * childOffset[2] as f64,
            );

            shiftBuffer.push(childOffset);

            ascentChildIndex = *CHILD_OFFSET_ORDER_INDEX.get(&childOffset)
                               .expect("Invalid child offset");
            ascentNodeIndex = self.childPointReferences[ascentNodeIndex][ascentChildIndex]
                              .expect("Failed to gather child point reference");

            self.depthIndexBufferSearch.push(ascentNodeIndex);
        }

        // this should be true? It should go too high or wide or anything right?
        self.GetChildFaceNodesRecursive(
            depth,
            self.childPointReferences[self.depthIndexBufferSearch[depth - 1]][
                    *CHILD_OFFSET_ORDER_INDEX.get(
                        &[1 - childOffset[0], childOffset[1], childOffset[2]]
                    ).expect("Invalid offset")
                ].expect("Failed to get child"),
            childOffset[0],
            0,
            1,
            2,
            depth,
            (self.offsetX, self.offsetY, self.offsetZ),
            (self.sizeX, self.sizeY, self.sizeZ)
        );

        // this should be true? It should go too high or wide or anything right?
        self.GetChildFaceNodesRecursive(
            depth,
            self.childPointReferences[self.depthIndexBufferSearch[depth - 1]][
                    *CHILD_OFFSET_ORDER_INDEX.get(
                        &[childOffset[0], 1 - childOffset[1], childOffset[2]]
                    ).expect("Invalid offset")
                ].expect("Failed to get child"),
            childOffset[1],
            1,
            0,
            2,
            depth,
            (self.offsetX, self.offsetY, self.offsetZ),
            (self.sizeX, self.sizeY, self.sizeZ)
        );

        // this should be true? It should go too high or wide or anything right?
        self.GetChildFaceNodesRecursive(
            depth,
            self.childPointReferences[self.depthIndexBufferSearch[depth - 1]][
                    *CHILD_OFFSET_ORDER_INDEX.get(
                        &[childOffset[0], childOffset[1], 1 - childOffset[2]]
                    ).expect("Invalid offset")
                ].expect("Failed to get child"),
            childOffset[2],
            2,
            0,
            1,
            depth,
            (self.offsetX, self.offsetY, self.offsetZ),
            (self.sizeX, self.sizeY, self.sizeZ)
        );

        self.FindChildFace(
            depth, &shiftBuffer, 0,
            (self.offsetX, base.1, base.2),
            (self.sizeX, nodeSize.1, nodeSize.2)
        );
        self.FindChildFace(
            depth, &shiftBuffer, 1,
            (base.0, self.offsetY, base.2),
            (nodeSize.0, self.sizeY, nodeSize.2)
        );
        self.FindChildFace(
            depth, &shiftBuffer, 2,
            (base.0, base.1, self.offsetZ),
            (nodeSize.0, nodeSize.1, self.sizeZ)
        );
    }

    fn FindChildFace (&mut self, depth: usize, shiftBuffer: &Vec <[usize; 3]>, axis: usize, pos: (f64, f64, f64), size: (f64, f64, f64)) {
        let lookingForShift = shiftBuffer[depth][axis];

        for d in (0..(depth-1)).rev() {
            if shiftBuffer[d][axis] == lookingForShift {
                let newChildIndex: usize;
                let offsetIndex: usize;
                let indexI: usize;
                let indexJ: usize;

                match axis {
                    0 => {
                        offsetIndex =
                            *CHILD_OFFSET_ORDER_INDEX.get(&[
                                lookingForShift, 1-shiftBuffer[d][1], 1-shiftBuffer[d][2]
                            ]).expect("Invalid shift");
                        indexI = 1; indexJ = 2;
                    },
                    1 => {
                        offsetIndex =
                            *CHILD_OFFSET_ORDER_INDEX.get(&[
                                1-shiftBuffer[d][0], lookingForShift, 1-shiftBuffer[d][2]
                            ]).expect("Invalid shift");
                        indexI = 0; indexJ = 2;
                    },
                    _ => {
                        offsetIndex =
                            *CHILD_OFFSET_ORDER_INDEX.get(&[
                                1-shiftBuffer[d][0], 1-shiftBuffer[d][1], lookingForShift
                            ]).expect("Invalid shift");
                        indexI = 0; indexJ = 1;
                    },
                }
                
                newChildIndex =
                    self.childPointReferences[self.depthIndexBufferSearch[d]]
                        [offsetIndex].expect("Failed to get child");
                
                self.GetChildFaceNodesRecursive (
                    depth,
                    newChildIndex,
                    1 - shiftBuffer[depth][axis],
                    axis,
                    indexI,
                    indexJ,
                    d,
                    pos,
                    size
                );

                return;
            }
        }
    }

    fn GetChildFaceNodesRecursive (
        &mut self, depth: usize, childIndex: usize,
        shift: usize, axis: usize,
        indexI: usize, indexJ: usize,
        currentDepth: usize,
        pos: (f64, f64, f64), size: (f64, f64, f64)
    ) {

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
                    let widthScalar = self.depthSizeScalars[currentDepth + 1];
                    let cellSize = (
                        self.sizeX * widthScalar,
                        self.sizeY * widthScalar,
                        self.sizeZ * widthScalar,
                    );
                    let cellPosition = self.GetNodePosition(childIndex);

                    if ((cellPosition.0 - cellSize.0*0.5 < pos.0 + size.0 && cellPosition.0 + cellSize.0*0.5 > pos.0 ) ||
                        (pos.0 < cellPosition.0 + cellSize.0*0.5 && pos.0 + size.0 > cellPosition.0 - cellSize.0*0.5)) &&
                       ((cellPosition.1 - cellSize.1*0.5 < pos.1 + size.1 && cellPosition.1 + cellSize.1*0.5 > pos.1 ) ||
                        (pos.1 < cellPosition.1 + cellSize.1*0.5 && pos.1 + size.1 > cellPosition.1 - cellSize.1*0.5)) &&
                       ((cellPosition.2 - cellSize.2*0.5 < pos.2 + size.2 && cellPosition.2 + cellSize.2*0.5 > pos.2) ||
                        (pos.2 < cellPosition.2 + cellSize.2*0.5 && pos.2 + size.2 > cellPosition.2 - cellSize.2*0.5)
                    ) {
                        self.cellNeighborReferences[self.depthIndexBufferSearch[depth]].
                                push(childIndex);
                    }
                    return;
                }

                self.GetChildFaceNodesRecursive(
                    depth,
                    newChildIndex.unwrap(),
                    shift,
                    axis,
                    indexI,
                    indexJ,
                    currentDepth + 1,
                    pos,
                    size
                );
            }
        }
    }


    pub fn GetLeafIndex (&mut self, inputPos: &(f64, f64, f64)) -> usize {
        let pos = (
            // magic number is for precision to prevent errors at the far extremes
            inputPos.0.max(self.offsetX).min(self.sizeX + self.offsetX - 0.00001),
            inputPos.1.max(self.offsetY).min(self.sizeY + self.offsetY - 0.00001),
            inputPos.2.max(self.offsetZ).min(self.sizeZ + self.offsetZ - 0.00001),
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

        for i in 1..=self.maxDepth {
            self.lastLeafDepth = i;

            widthScalar = self.depthSizeScalars[i];  // is this right?????? -1 or no -1?????
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

            //println!("({}) Offset: {}, {}, {}    Pos: {}, {}, {}    Base: {}, {}, {}    Size: {}, {}, {}", i, childOffset.0, childOffset.1, childOffset.2, pos.0, pos.1, pos.2, base.0, base.1, base.2, nodeSize.0, nodeSize.1, nodeSize.2);

            childIndex = *CHILD_OFFSET_ORDER_INDEX.get(&[childOffset.0, childOffset.1, childOffset.2]).
                expect("Invalid offset");
            
            base = (
                base.0 + nodeSize.0 * childOffset.0 as f64,
                base.1 + nodeSize.1 * childOffset.1 as f64,
                base.2 + nodeSize.2 * childOffset.2 as f64,
            );

            lastNodeIndex = nodeIndex;
            nodeIndex = self.childPointReferences[nodeIndex.unwrap()][childIndex];
            if nodeIndex.is_none() {
                nodeIndex = lastNodeIndex;
                break;
            }

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

        //let mut numSeen = 1;
        //let mut numPopped = 0;

        while !queue.IsEmpty() {
            //numPopped += 1;
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
                //println!("Seen: {}", numSeen);
                return Some(minDst.sqrt())
            }

            visited.insert(currentNode.1, true);

            for positionIndex in &self.positionIndexes[currentNode.1] {
                let point = points[*positionIndex];
                distance = self.GetSquaredDistance(point, samplePosition);
                minDst = minDst.min(distance);
            }

            for neighborIndex in &self.cellNeighborReferences[currentNode.1] {
                if visited.contains_key(neighborIndex) {
                    continue;
                }

                distance = self.GetSquaredDistance(
                    self.GetNodePosition(*neighborIndex),
                    queryPoint
                );

                // hopefully this isn't broken and works. Should speed up the program bc/ the math
                // I'd think is much quicker than the insertion algorithm
                /*depthScalar = self.depthSizeScalars[
                self.leafNodeDepths[*neighborIndex]
                    .expect("Not a leaf node")];
                nodeWidthXYZ = (
                    self.sizeX * depthScalar,
                    self.sizeY * depthScalar,
                    self.sizeZ * depthScalar
                );
                nodeWidth = nodeWidthXYZ.0*nodeWidthXYZ.0 +
                            nodeWidthXYZ.1*nodeWidthXYZ.1 +
                            nodeWidthXYZ.2*nodeWidthXYZ.2;
                
                minNodeDst = distance + nodeWidth * 0.25 -
                            (distance * nodeWidth).sqrt();
                
                if minNodeDst < minDst || true {  // remove || true  once I know if it works or not
                    queue.Push((
                        distance, *neighborIndex
                    ));
                }// */
                //numSeen += 1;
                queue.Push((
                    distance, *neighborIndex
                ));
            }
        }

        //println!("Ran out; Seen: {}      MD: {}      Popped: {}", numSeen, minDst.sqrt(), numPopped);

        Some(minDst.sqrt())
        //None    ig it can actually end up searching every node.....
    }


    pub fn GetLeafNodes (&self) -> Vec <usize> {
        self.leafNodes.clone()
    }

    pub fn GetCornerPositionIndex (&self, leafNodeIndex: usize, cornerIndex: usize) -> usize {
        // not sure if the keys are in order; if so use a binary search for better performance
        let mortonCode = self.nodeCornerReferences[leafNodeIndex][cornerIndex];
        for (i, key) in self.cornerPoints.keys().enumerate() {
            if *key == mortonCode {
                return i;
            }
        } 0  // error?  shouldn't ever reach here? I give up...  :(
    }

    pub fn GetNodePosition (&self, nodeIndex: usize) -> (f64, f64, f64) {
        self.nodePositions[nodeIndex]
        //if let Some(node) = self.nodePositions.get(nodeIndex) {
        //    if node.is_none() { return None; }
        //    return *node;
        //} None
    }

    pub fn GetLeafNodeReverseIndex (&self, leafNode: usize) -> usize {
        self.leafNodes.binary_search(&leafNode).expect("failed to get the point")
        /*BinarySearch(&self.leafNodes, &leafNode).  // not sure if this is working?
            expect("Failed to get point")*/
    }

    // does this need to be public?
    pub fn GetMortonCode (&self, xi: u32, yi: u32, zi: u32) -> u128 {
        // getting the morton code
        let mut x = xi;
        let mut y = yi;
        let mut z = zi;
        
        // Interleave x bits (starting from the least significant bits)
        x = (x | (x << 16)) & 0x030000FF;
        x = (x | (x << 8)) & 0x0300F00F;
        x = (x | (x << 4)) & 0x030C30C3;
        x = (x | (x << 2)) & 0x09249249;
        
        // Interleave y bits (shifted by 1 bit to align with interleaving of x)
        y = (y | (y << 16)) & 0x030000FF;
        y = (y | (y << 8)) & 0x0300F00F;
        y = (y | (y << 4)) & 0x030C30C3;
        y = (y | (y << 2)) & 0x09249249;
        
        // Interleave z bits (shifted by 2 bits to align with interleaving of x)
        z = (z | (z << 16)) & 0x030000FF;
        z = (z | (z << 8)) & 0x0300F00F;
        z = (z | (z << 4)) & 0x030C30C3;
        z = (z | (z << 2)) & 0x09249249;
        
        // Combine the interleaved bits of x, y, and z
        (x as u128) | ((y as u128) << 1) | ((z as u128) << 2)
        /*let mut mortonCode: u128 = 0;
        for i in 0..32 {  // 64-bit integer
            mortonCode |= (((xi >> i) & 1) as u128) << (3 * i);
            mortonCode |= (((y >> i) & 1) as u128) << (3 * i + 1);
            mortonCode |= (((z >> i) & 1) as u128) << (3 * i + 2);
        } mortonCode*/
    }

    // has to clone them sadly.... (not sure how else to do it)
    pub fn GetNodeNeighbors (&self, nodeIndex: usize) -> Vec <usize> {
        self.cellNeighborReferences[nodeIndex].clone()
    }

}


