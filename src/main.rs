/*
Fix sign calculation :(

Improve the octree subdivision system to be smarter about searching for intersections with boxes
Start using the node corners for the array in distance field calculations
Impliment dual contouring or something like it for water tight iso contour extraction

Add distance check before pushing items to the priority queue
    Many of the items being pushed could be cancled out early
    It costs one square root operation and saves many, many checks

Finish fixing the traversal caching
    The bug is narrowed down to the culling system in the final stage
    All neighbors are being correctly generated it seems (test size 1 node so errors are possible)
    However the wrong neighbors are being culled, and some aren't being culled

*/

// snake case is just bad
#![allow(non_snake_case)]

// some of the constants use this to make them readable
#![allow(clippy::eq_op)]

mod priorityQueue;
mod octree;

const LOADING_DISTANCE_FIELD_SAVE: bool = false;
const GENERATE_NEW_VERSION_SDF: bool = false;//true;
const USING_OLD_SDF: bool = true;

// commented out any items that have no implimentation or
// are never checked (todo's)
//const FIELD_FILE: &str = "";
//const SAVE_FILE: &str = "";

//const OBJ_SAVE_FILE: &str = "";

const LOAD_PCD: bool = false;
const WRITE_OCTREE_PCD: bool = true;
const OCTREE_PCD_FILE: &str = "assets/octree_neighbors.ply";

const WRITE_OUTPUT_PCD: bool = false;
const OCTREE_OUTPUT_PCD_FILE: &str = "assets/octree_output_testhallowsphere.ply";
//const PCD_FILE: &str = "";
//const POINTS_LOADED_PERCENT: f64 = 100.0 / (100.0);

const GENERATE_HALLOW_SPHERE: bool = true;
const GENERATE_SOLID_SPHERE: bool = false;
const GENERATE_CUBE: bool = false;

const MAX_OCTREE_DEPTH: usize = 10;

const SIMULATING_SURFACE_TENSION: bool = false;
/*const TENSION_ITTERATIONS: isize = 2500;
const VELOCITY: f64 = 0.75;
const DT: f64 = 0.0055 * 0.1;*/


const MAX_SAMPLING_SPACE: usize = GRID_SIZE.0 * GRID_SIZE.1 * GRID_SIZE.2;
const GRID_SIZE: (usize, usize, usize) = (50, 50, 50);


const AUTO_SET: bool = true;



//=========================================================================================================
//                                          Sign Calculations
//=========================================================================================================


fn CalculateSigns (setupParameters: &SetupParameters,
              octree: &mut octree::Octree,
              distanceField: &Vec <(f64, u128)>,
              signedField: &mut Vec <i8>,
              numberOfNodes: usize,
    ) {
    
    // umm...... search algerithm still has a couple random artifacts (3 in this case)    still has some? but it's better
    // the chaching system, or searching system may be borken...... :(    probably still true......
    // probably becoming random from overflow errors; very sad;       nope, it's in here. No overflow and the total dst is the same every time
    let mut cummulative = 0.0;
    for dst in distanceField {
        cummulative += dst.0;
    }

    println!("Total dist: {}", cummulative);

    // getting the base node at the corner
    let leafIndex = octree.GetLeafIndex(&(
        setupParameters.SampleSpaceBounds[0] - setupParameters.SampleSpaceOffset[0],
        setupParameters.SampleSpaceBounds[1] - setupParameters.SampleSpaceOffset[1],
        setupParameters.SampleSpaceBounds[2] - setupParameters.SampleSpaceOffset[2],
    ));
    let startingIndex = octree.GetCornerPositionIndex(leafIndex, 0);

    // the starting index is random but that makes sense considering the order of a hashmap is non-deterministic
    println!("Starting Leaf: {}", startingIndex);

    // creating buffers for solid points, and hallow points
    let mut solidEdgesBuffer: Vec <(usize, usize)> = vec!();
    let mut voidPointsBuffer = vec![
        (startingIndex, leafIndex)
    ];

    let mut sign = 1;

    let mut numSolid = 0;
    let mut numHallow = 0;
    let mut numWall = 0;
    let mut numFound = 0;

    let mut newPointsBuffer: Vec <(usize, usize)> = vec!();

    let mut leaf: usize;
    let mut cornerNode: usize;

    let mut taken: std::collections::HashMap <usize, bool> = std::collections::HashMap::new();
    taken.insert(startingIndex, true);

    // starting the itteration and going until every position has been filled
    let mut i = 0;
    while !voidPointsBuffer.is_empty() {
        println!("Itteration {}", i);
        i += 1;
        // emptying void points buffer
        while !newPointsBuffer.is_empty() || voidPointsBuffer.len() > 0 {
            while let Some(newVoidPoint) = newPointsBuffer.pop() {
                voidPointsBuffer.push(newVoidPoint);
            }
            
            println!("Len void: {}    Added: {}", voidPointsBuffer.len(), numFound);
            while let Some((nodeIndex, leafIndex)) = voidPointsBuffer.pop() {
                // updating the sign
                if signedField[nodeIndex] != 0 {  continue;  }
                
                //println!("({}): signedField: {}", nodeIndex, signedField[nodeIndex]);
                signedField[
                    nodeIndex
                ] = sign;

                numFound += 1;

                // getting neighboring cells
                for neighborIndex in octree.GetNodeNeighbors(leafIndex) {
                    leaf = octree.GetLeafNodeReverseIndex(neighborIndex);
                    for i in 0..8usize {
                        cornerNode = octree.GetCornerPositionIndex(leaf, i);
                        if taken.contains_key(&cornerNode) {  continue;  }
                        taken.insert(cornerNode, true);
                        
                        if distanceField[cornerNode].0 < setupParameters.IsoContourLevel {
                            solidEdgesBuffer.push((cornerNode, neighborIndex));
                            numSolid += 1;
                        } else {
                            newPointsBuffer.push((cornerNode, neighborIndex));
                            if sign == 1 {numHallow += 1;}
                            else {numWall += 1;}
                        }
                    }
                }
            }
        }
        
        println!("Past!!!");
        // emptying the solid edges buffer
        while !newPointsBuffer.is_empty() || solidEdgesBuffer.len() > 0 {
            while let Some(newEdgePoint) = newPointsBuffer.pop() {
                solidEdgesBuffer.push(newEdgePoint);
            }
            println!("Len void: {}    Added: {}", solidEdgesBuffer.len(), numFound);

            while let Some((nodeIndex, leafIndex)) = solidEdgesBuffer.pop() {
                if signedField[nodeIndex] != 0 {  continue;  }
                
                signedField[
                    nodeIndex
                ] = -1;

                numFound += 1;

                for neighborIndex in octree.GetNodeNeighbors(leafIndex) {
                    leaf = octree.GetLeafNodeReverseIndex(neighborIndex);
                    for i in 0..8usize {
                        cornerNode = octree.GetCornerPositionIndex(leaf, i);
                        if taken.contains_key(&cornerNode) {  continue;  }
                        taken.insert(cornerNode, true);
                        
                        if distanceField[cornerNode].0 < setupParameters.IsoContourLevel {
                            newPointsBuffer.push((cornerNode, neighborIndex));
                            numSolid += 1;
                        } else {
                            voidPointsBuffer.push((cornerNode, neighborIndex));
                            if sign == 1 {numHallow += 1;}
                            else {numWall += 1;}
                        }
                    }
                }
            }
        }

        sign *= -1;
    }

    println!("({}) Solid: {}    Voids;  Hallow: {}  Solid: {}", numFound, numSolid, numHallow, numWall);
}



/*
// the offsets for neighbors
const FLOOD_FILL_NEIGHBOR_OFFSETS: [(isize, isize, isize); 6] = [
    (1 , 0 , 0 ),
    (-1, 0 , 0 ),
    (0 , 1 , 0 ),
    (0 , -1, 0 ),
    (0 , 0 , 1 ),
    (0 , 0 , -1),
];

fn FloodFill (setupParameters: &SetupParameters,
              (startX, startY, startZ): (usize, usize, usize),
              distanceField: &[f64; MAX_SAMPLING_SPACE],
              signedGrid: &mut [i8; MAX_SAMPLING_SPACE],
              sign: i8) {

    // the stored neighboring points
    let mut neighbors: Vec <(usize, usize, usize)> = vec!();
    let mut newNeighbors: Vec <(usize, usize, usize)> = vec!();

    //let mut numNeighborCopies = [0usize; MAX_SAMPLING_SPACE];

    neighbors.push(
        (startX, startY, startZ)
    );
    let mut numNeighbors = 1usize;

    while numNeighbors > 0 {
        let mut numNewNeighbors = 0usize;

        for neighbor in &neighbors {
            let (indexX, indexY, indexZ) = *neighbor;
            let arrayIndex = indexX + (indexY + indexZ * GRID_SIZE.1) * GRID_SIZE.0;
            if distanceField[arrayIndex] >= setupParameters.IsoContourLevel{
                signedGrid[arrayIndex] = sign;

                for (offsetX, offsetY, offsetZ) in FLOOD_FILL_NEIGHBOR_OFFSETS {
                    let neighborXPos = (indexX as isize + offsetX).unsigned_abs();
                    let neighborYPos = (indexX as isize + offsetY).unsigned_abs();
                    let neighborZPos = (indexX as isize + offsetZ).unsigned_abs();

                    if neighborXPos < GRID_SIZE.0 &&
                        neighborYPos < GRID_SIZE.1 &&
                        neighborZPos < GRID_SIZE.2 {

                        let indexKey = neighborXPos + (neighborYPos + neighborZPos * GRID_SIZE.1) * GRID_SIZE.0;
                        if signedGrid[indexKey] == 0 {
                            numNewNeighbors += 1;
                            newNeighbors.push(
                                (neighborXPos, neighborYPos, neighborZPos)
                            );
                        }
                    }
                }
            }
        }
        
        if numNewNeighbors == 0 { return; }

        numNeighbors = numNewNeighbors;
        neighbors.clear();

        while let Some(neighbor) = newNeighbors.pop() {
            neighbors.push(neighbor);
        }
    }
}


fn CalculateSigns (setupParameters: &SetupParameters,
                   distanceField: &[f64; MAX_SAMPLING_SPACE],
                   signedGrid: &mut [i8; MAX_SAMPLING_SPACE]) {

    FloodFill(setupParameters, (0, 0, 0), distanceField, signedGrid, 1);

    for x in 0..GRID_SIZE.0 {
        for y in 0..GRID_SIZE.1 {
            let mut sign = 1i8;

            for z in 0..GRID_SIZE.2 {
                let index = x + (y + z * GRID_SIZE.1) * GRID_SIZE.0;
                if signedGrid[index] != 0 {
                    sign = signedGrid[index];
                } else if distanceField[index] >= setupParameters.IsoContourLevel {
                    sign *= -1;
                    FloodFill(setupParameters, (x, y, z), distanceField, signedGrid, sign);
                }
            }
        }
    }
}*/




//=========================================================================================================
//                                       Data Loading & Saving
//=========================================================================================================


// data loading functions
// add them here (don't feel like figuring it out right now)




//=========================================================================================================
//                                          Shape Generation
//=========================================================================================================


fn FibonacciSphere (samples: isize, pointCloud: &mut Vec <(f64, f64, f64)>,
                    scalingFactor: f64, offset: (f64, f64, f64)) {
    
    let phi = std::f64::consts::PI * (std::f64::consts::FRAC_1_SQRT_2 - 1.0);  // the magic number is sqrt 0.5

    for i in 0..samples {
        let y = 1.0 - (i as f64 / (samples - 1) as f64) * 2.0;
        let radius = (1.0 - y*y).sqrt();

        let theta = phi * i as f64;

        let x = theta.cos() * radius;
        let z = theta.sin() * radius;

        pointCloud.push((
            x * scalingFactor + offset.0,
            y * scalingFactor + offset.1,
            z * scalingFactor + offset.2,
        ));
    }
}



//=========================================================================================================
//                                          Main Script
//=========================================================================================================


pub struct SetupParameters {
    SampleSpaceBounds: [f64; 3],
    SampleSpaceOffset: [f64; 3],
    IsoContourLevel: f64,
    InverseDeltaX: f64,
    InverseDeltaY: f64,
    InverseDeltaZ: f64,
}


fn main() {

    let mut initoctree: Option <octree::Octree> = None;

    let mut setupParameters = SetupParameters {
        SampleSpaceBounds: [0.0, 0.0, 0.0],
        SampleSpaceOffset: [0.0, 0.0, 0.0],
        IsoContourLevel: 4.75,// * 0.25,
        InverseDeltaX: 1.0 / (1.0),
        InverseDeltaY: 1.0 / (1.0),
        InverseDeltaZ: 1.0 / (1.0),
    };

    //let mut distanceField = [0.0f64; MAX_SAMPLING_SPACE];
    let mut octreeDistanceGrid: Vec <(f64, u128)> = vec!();

    if LOADING_DISTANCE_FIELD_SAVE {
        //...
    } else {
        let mut pointCloud: Vec <(f64, f64, f64)> = vec!();

        if LOAD_PCD {
            //...
        } else {
            if GENERATE_HALLOW_SPHERE {
                FibonacciSphere(575*1, &mut pointCloud, 34.0, (50.0, 50.0, 50.0));
                FibonacciSphere(250*1, &mut pointCloud, 15.0, (50.0, 50.0, 50.0));
            }
            if GENERATE_SOLID_SPHERE {
                FibonacciSphere(200, &mut pointCloud, 14.0, (50.0, 50.0, 50.0));
            }
            if GENERATE_CUBE {
                for x in 0..8 {
                    for y in 0..8 {
                        for z in 0..8 {
                            if (x == 7 || y == 7 || z == 7) &&
                               !pointCloud.contains(&(x as f64, y as f64, z as f64)) {
                                pointCloud.push(
                                    (x as f64, y as f64, z as f64)
                                );
                            }
                        }
                    }
                }
            }

            if AUTO_SET {
                let (mut lowestX, mut lowestY, mut lowestZ): (f64, f64, f64) =
                    (f64::MAX, f64::MAX, f64::MAX);
                let (mut largestX, mut largestY, mut largestZ): (f64, f64, f64) =
                    (f64::MIN, f64::MIN, f64::MIN);

                for (x, y, z) in &pointCloud {
                    lowestX = lowestX.min(*x);
                    lowestY = lowestY.min(*y);
                    lowestZ = lowestZ.min(*z);

                    largestX = largestX.max(*x);
                    largestY = largestY.max(*y);
                    largestZ = largestZ.max(*z);
                }

                let mut cloudSizeX = largestX - lowestX;
                let mut cloudSizeY = largestY - lowestY;
                let mut cloudSizeZ = largestZ - lowestZ;

                setupParameters.InverseDeltaX = cloudSizeX / GRID_SIZE.0 as f64;
                setupParameters.InverseDeltaY = cloudSizeY / GRID_SIZE.1 as f64;
                setupParameters.InverseDeltaZ = cloudSizeZ / GRID_SIZE.2 as f64;

                let newIsoContourLevel = setupParameters.IsoContourLevel * (
                    setupParameters.InverseDeltaX + setupParameters.InverseDeltaY + setupParameters.InverseDeltaZ
                ) / 3.0;

                cloudSizeX += newIsoContourLevel * 4.0;
                cloudSizeY += newIsoContourLevel * 4.0;
                cloudSizeZ += newIsoContourLevel * 4.0;

                setupParameters.SampleSpaceOffset = [
                    lowestX - newIsoContourLevel * 2.0,
                    lowestY - newIsoContourLevel * 2.0,
                    lowestZ - newIsoContourLevel * 2.0,
                ];

                setupParameters.InverseDeltaX = cloudSizeX / GRID_SIZE.0 as f64;
                setupParameters.InverseDeltaY = cloudSizeY / GRID_SIZE.1 as f64;
                setupParameters.InverseDeltaZ = cloudSizeZ / GRID_SIZE.2 as f64;

                let dAvg = (setupParameters.InverseDeltaX + setupParameters.InverseDeltaY + setupParameters.InverseDeltaZ) / 3.0;
                setupParameters.IsoContourLevel *= dAvg;

                setupParameters.InverseDeltaX = 1.0 / setupParameters.InverseDeltaX;
                setupParameters.InverseDeltaY = 1.0 / setupParameters.InverseDeltaY;
                setupParameters.InverseDeltaZ = 1.0 / setupParameters.InverseDeltaZ;

                setupParameters.SampleSpaceBounds = [
                    cloudSizeX, cloudSizeY, cloudSizeZ
                ];

                println!("Cloud size: {}, {}, {}", GRID_SIZE.0, GRID_SIZE.1, GRID_SIZE.2);
                println!("Cloud offset: {}, {}, {}", setupParameters.SampleSpaceOffset[0], setupParameters.SampleSpaceOffset[1], setupParameters.SampleSpaceOffset[2]);
                println!("Deltas: {}, {}, {}", 1.0/setupParameters.InverseDeltaX, 1.0/setupParameters.InverseDeltaY, 1.0/setupParameters.InverseDeltaZ);
                println!("Lowest: {}, {}, {}", lowestX, lowestY, lowestZ);
                println!("largest: {}, {}, {}", largestX, largestY, largestZ);
                println!("Iso contour level: {}", setupParameters.IsoContourLevel);
            }

            let mut octree = octree::Octree::new(
                (setupParameters.SampleSpaceOffset[0], setupParameters.SampleSpaceOffset[1], setupParameters.SampleSpaceOffset[2]),
                (setupParameters.SampleSpaceBounds[0], setupParameters.SampleSpaceBounds[1], setupParameters.SampleSpaceBounds[2]),
                MAX_OCTREE_DEPTH
            );
            octree.SubDivide(&pointCloud);
            octree.GenerateNodeReferenceCache();
            let numCorners = octree.GenerateCornerPointsCache();
            println!("Number of corners: {}", numCorners);

            let ip = (35.23, 43.35, 54.33);
            let leaf = octree.GetLeafIndex(&ip);
            
            if WRITE_OCTREE_PCD {
                // outputing the points to a point cloud file
                let samples = 250; // 250
                let sphereSize = 0.5;  // 0.05
                let mut outputString = format!("ply\nformat ascii 1.0\nelement vertex {}\nproperty float x\nproperty float y\nproperty float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nend_header\n", numCorners*(samples+1)).to_string();
                for corner in octree.GetCornerPointKeys() {
                    let corner = octree.LookupCornerPoint(corner);
                    // generating a sphere around the point so it renders bigger
                    for corner in octree.GetAllNeighborCornersTemp(leaf) {
                    let mut sphere: Vec<(f64, f64, f64)> = vec![
                        corner
                    ];
                    FibonacciSphere(samples as isize, &mut sphere, sphereSize, corner);
                    for point in sphere {
                        outputString += &format!("{} {} {} 100 25 25\n", point.0, point.1, point.2);
                    }} break;
                }
                std::fs::write(OCTREE_PCD_FILE, outputString).expect("Unable to write file");
            }

            for cornerId in octree.GetCornerPointKeys() {
                let samplePosition = octree.LookupCornerPoint(cornerId);
                let leafNodeIndex = octree.GetLeafIndex(&(
                    samplePosition.0 - 0.001,
                    samplePosition.1 - 0.001,
                    samplePosition.2 - 0.001
                ));
                let dst = octree.NearestNeighborSearch(&pointCloud, leafNodeIndex, samplePosition).
                expect("Failed to get distance");
                octreeDistanceGrid.push(
                    (octree.NearestNeighborSearch(&pointCloud, leafNodeIndex, samplePosition).
                        expect("Failed to get distance"),
                    cornerId)
                );
                //println!("dst: {}", dst);
            }

            let mut octreeSignedGrid: Vec <i8> = vec![0i8; numCorners];
            CalculateSigns(&setupParameters, &mut octree, &octreeDistanceGrid, &mut octreeSignedGrid, numCorners);
            for (index, point) in octreeDistanceGrid.iter_mut().enumerate() {
                point.0 *= octreeSignedGrid[index] as f64;
            }

            /*let corners = octree.GetAllNeighborConersTemp(leafNodeIndex);
            for corner in corners {
                println!("[{}, {}, {}],", corner.0, corner.1, corner.2);
            }*/
            
            /*for _i in 0..1 {
                let dst = octree.NearestNeighborSearch(
                    &pointCloud, leafNodeIndex, inputPos
                );
                println!("dst: {}", dst.expect("Failed to get distance"));
            }*/

            //let mut octreeNeighborDistance: Vec <f64> = vec!();     would need editing....
            
            /*for node in octree.GetLeafNodes() { .   will need editing and be based on corners instead
                let position = octree.GetLeafPosition(node).
                            expect("Failed to get position");
                let leafNodeIndex = octree.GetLeafIndex(position);
                octreeNeighborDistance.push(
                    octree.NearestNeighborSearch(&pointCloud,
                                                leafNodeIndex,
                                                position).
                                                expect("Failed to find distance")
                );
            }*/

            /*for x_ in 0..GRID_SIZE.0 {
                for y_ in 0..GRID_SIZE.1 {
                    for z_ in 0..GRID_SIZE.2 {
                        let xyz = (
                            x_ as f64 / setupParameters.InverseDeltaX + setupParameters.SampleSpaceOffset[0],
                            y_ as f64 / setupParameters.InverseDeltaY + setupParameters.SampleSpaceOffset[1],
                            z_ as f64 / setupParameters.InverseDeltaZ + setupParameters.SampleSpaceOffset[2],
                        );
                        
                        let leafNodeIndex = octree.GetLeafIndex(&xyz);
                        let distance = octree.NearestNeighborSearch(&pointCloud, leafNodeIndex, xyz);
                        distanceField[x_ + (y_ + z_ * GRID_SIZE.1) * GRID_SIZE.0] =
                            distance.expect("Failed to get distance");
                    }
                }
            }*/

            // temporary
            initoctree = Some(octree);

        }

        /*if USING_OLD_SDF {
            for x_ in 0..GRID_SIZE.0 {
                for y_ in 0..GRID_SIZE.1 {
                    for z_ in 0..GRID_SIZE.2 {
                        let index = x_ + (y_ + z_ * GRID_SIZE.1) * GRID_SIZE.0;
                        distanceField[index] = distanceField[index].abs();
                    }
                }
            }
        }

        if GENERATE_NEW_VERSION_SDF {
            
            let mut signedGrid = [0i8; MAX_SAMPLING_SPACE];

            // sign calculations are entirely borken..... (this error is propagating)
            CalculateSigns(&setupParameters, &distanceField, &mut signedGrid);
            
            let mut surfacePoints: Vec <(f64, f64, f64)> = vec!();

            for x in 1..GRID_SIZE.0-1 {
                for y in 1..GRID_SIZE.1-1 {
                    for z in 1..GRID_SIZE.2-1 {
                        if signedGrid[(x - 1) + (y + z * GRID_SIZE.1) * GRID_SIZE.0] > 0 ||
                           signedGrid[(x + 1) + (y + z * GRID_SIZE.1) * GRID_SIZE.0] > 0 ||
                           signedGrid[x + ((y - 1) + z * GRID_SIZE.1) * GRID_SIZE.0] > 0 ||
                           signedGrid[x + ((y + 1) + z * GRID_SIZE.1) * GRID_SIZE.0] > 0 ||
                           signedGrid[x + (y + (z - 1) * GRID_SIZE.1) * GRID_SIZE.0] > 0 ||
                           signedGrid[x + (y + (z + 1) * GRID_SIZE.1) * GRID_SIZE.0] > 0
                           {
                            let index = x + (y + z * GRID_SIZE.1) * GRID_SIZE.0;
                            if signedGrid[index] < 0 || distanceField[index] < setupParameters.IsoContourLevel {
                                surfacePoints.push(
                                    (x as f64, y as f64, z as f64)
                                );
                            }
                        }
                    }
                }
            }

            println!("Ready for new one: {}", surfacePoints.len());
            let mut octree = octree::Octree::new(
                (setupParameters.SampleSpaceOffset[0], setupParameters.SampleSpaceOffset[1], setupParameters.SampleSpaceOffset[2]),
                (setupParameters.SampleSpaceBounds[0], setupParameters.SampleSpaceBounds[1], setupParameters.SampleSpaceBounds[2]),
                MAX_OCTREE_DEPTH
            );
            octree.SubDivide(&surfacePoints);
            octree.GenerateNodeReferenceChache();
            let numCorners = octree.GenerateCornerPointsChache();
            println!("Number of corners: {}", numCorners);

            let oldIsoContourLevel = setupParameters.IsoContourLevel;
            if AUTO_SET && !LOADING_DISTANCE_FIELD_SAVE {
                setupParameters.IsoContourLevel *= 3.0 /
                    (1.0/setupParameters.InverseDeltaX +
                     1.0/setupParameters.InverseDeltaY +
                     1.0/setupParameters.InverseDeltaZ);
            }

            for x in 0..GRID_SIZE.0 {
                for y in 0..GRID_SIZE.1 {
                    for z in 0..GRID_SIZE.2 {
                        let position = (x as f64, y as f64, z as f64);
                        let leafNodeIndex = octree.GetLeafIndex(&position);
                        let distance = octree.NearestNeighborSearch(&surfacePoints,
                                leafNodeIndex, position).
                                expect("Failed to get distance");
                        
                        let index = x + (y + z * GRID_SIZE.1) * GRID_SIZE.0;
                        let mut sign = signedGrid[index];
                        
                        if distanceField[index] < oldIsoContourLevel && sign > -1 {
                            sign = -1;
                        }
                        if sign == 0 {
                            sign = 1;
                        }

                        distanceField[index] = distance *
                                    (sign * 0 + 1) as f64 +
                                    setupParameters.IsoContourLevel;
                    }
                }
            }

            if SIMULATING_SURFACE_TENSION {
                // do something
            }

            if !GENERATE_NEW_VERSION_SDF && USING_OLD_SDF {
                // do something
            }

            // do marching cubes or dual contouring... (todo)

        }*/
    }

    // testing rendering
    /*let slicePositionZ = 25usize;
    for x in 0..GRID_SIZE.0 {
        for y in 0..GRID_SIZE.1 {
            //if distanceField[x + (y + slicePositionZ * GRID_SIZE.1) * GRID_SIZE.0] < setupParameters.IsoContourLevel * 0.5 {
            //    print!("##");
            //} else if distanceField[x + (y + slicePositionZ * GRID_SIZE.1) * GRID_SIZE.0] < setupParameters.IsoContourLevel {
            //    print!("//");
            //} else {
            //    print!("  ");
            //}
            let dst = distanceField[x + (y + slicePositionZ * GRID_SIZE.1) * GRID_SIZE.0];
            print!("{}{} ", std::cmp::min((dst/setupParameters.IsoContourLevel*2.0) as usize, 9), std::cmp::min((dst/setupParameters.IsoContourLevel*2.0) as usize, 9));
        }
        println!();
    }*/


    // finish this... should provide some visualization of what's going on
    // saving the result (as a point cloud for now...)
    if WRITE_OUTPUT_PCD {
        let mut numSolid = 0;
        // outputing the points to a point cloud file
        let samples = 250; // 250    5
        let sphereSize = 0.5;  // 0.05     0.1
        let mut outputString = "".to_string();
        let cornersTree = initoctree.unwrap();
        for (distance, id) in octreeDistanceGrid {
            
            if distance < setupParameters.IsoContourLevel {
                let cornerPosition = cornersTree.LookupCornerPoint(id);

                // generating a sphere around the point so it renders bigger
                let mut sphere: Vec<(f64, f64, f64)> = vec![
                    cornerPosition
                ];
                FibonacciSphere(samples as isize, &mut sphere, sphereSize, cornerPosition);
                for point in sphere {
                    outputString += &format!("{} {} {} 100 25 25\n", point.0, point.1, point.2);
                }
                numSolid += 1
            }
        }

        let header = format!("ply\nformat ascii 1.0\nelement vertex {}\nproperty float x\nproperty float y\nproperty float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nend_header\n", numSolid*(samples+1));
        let finalString = header + outputString.as_str();

        println!("Num points: {}", numSolid * (samples + 1));

        std::fs::write(OCTREE_OUTPUT_PCD_FILE, finalString).expect("Unable to write file");
    }

}


