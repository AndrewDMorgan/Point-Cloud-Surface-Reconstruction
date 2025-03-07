/*
Improve the octree subdivision system to be smarter about searching for intersections with boxes
Start using the node corners for the array in distance field calculations
Impliment dual contouring or something like it for water tight iso contour extraction
*/

// snake case is just bad
#![allow(non_snake_case)]

// some of the constants use this to make them readable
#![allow(clippy::eq_op)]

mod priorityQueue;
mod octree;

const LOADING_DISTANCE_FIELD_SAVE: bool = false;
const GENERATE_NEW_VERSION_SDF: bool = true;
const USING_OLD_SDF: bool = true;

// commented out any items that have no implimentation or
// are never checked (todo's)
//const FIELD_FILE: &str = "";
//const SAVE_FILE: &str = "";

//const OBJ_SAVE_FILE: &str = "";

const LOAD_PCD: bool = false;
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
}




//=========================================================================================================
//                                       Data Loading & Saving
//=========================================================================================================


// data loading functions
// add them here (don't feel like figuring it out right now)




//=========================================================================================================
//                                          Shape Generation
//=========================================================================================================


fn FibonacciSphere (samples: isize, pointCloud: &mut Vec <(f64, f64, f64)>,
                    scalingFactor: f64, offset: f64) {
    
    let phi = std::f64::consts::PI * (std::f64::consts::FRAC_1_SQRT_2 - 1.0);  // the magic number is sqrt 0.5

    for i in 0..samples {
        let y = 1.0 - (i as f64 / (samples - 1) as f64) * 2.0;
        let radius = (1.0 - y*y).sqrt();

        let theta = phi * i as f64;

        let x = theta.cos() * radius;
        let z = theta.sin() * radius;

        pointCloud.push((
            x * scalingFactor + offset,
            y * scalingFactor + offset,
            z * scalingFactor + offset,
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

    let mut setupParameters = SetupParameters {
        SampleSpaceBounds: [101.0, 101.0, 101.0],
        SampleSpaceOffset: [0.0, 0.0, 0.0],
        IsoContourLevel: 6.0,
        InverseDeltaX: 1.0 / (1.0),
        InverseDeltaY: 1.0 / (1.0),
        InverseDeltaZ: 1.0 / (1.0),
    };

    let mut distanceField = [0.0f64; MAX_SAMPLING_SPACE];

    if LOADING_DISTANCE_FIELD_SAVE {
        //...
    } else {
        let mut pointCloud: Vec <(f64, f64, f64)> = vec!();

        if LOAD_PCD {
            //...
        } else {
            if GENERATE_HALLOW_SPHERE {
                FibonacciSphere(575, &mut pointCloud, 34.0, 50.0);
                FibonacciSphere(250, &mut pointCloud, 15.0, 50.0);
            } if GENERATE_SOLID_SPHERE {
                FibonacciSphere(200, &mut pointCloud, 14.0, 50.0);
            } if GENERATE_CUBE {
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

                setupParameters.InverseDeltaX = cloudSizeX / setupParameters.SampleSpaceBounds[0];
                setupParameters.InverseDeltaY = cloudSizeY / setupParameters.SampleSpaceBounds[1];
                setupParameters.InverseDeltaZ = cloudSizeZ / setupParameters.SampleSpaceBounds[2];

                setupParameters.IsoContourLevel *= (
                    setupParameters.InverseDeltaX + setupParameters.InverseDeltaY + setupParameters.InverseDeltaZ
                ) / 3.0;

                cloudSizeX += setupParameters.IsoContourLevel * 4.0;
                cloudSizeY += setupParameters.IsoContourLevel * 4.0;
                cloudSizeZ += setupParameters.IsoContourLevel * 4.0;

                setupParameters.SampleSpaceOffset = [
                    lowestX - setupParameters.IsoContourLevel * 2.0,
                    lowestY - setupParameters.IsoContourLevel * 2.0,
                    lowestZ - setupParameters.IsoContourLevel * 2.0,
                ];

                let dAvg = (setupParameters.InverseDeltaX + setupParameters.InverseDeltaY + setupParameters.InverseDeltaZ) / 3.0;
                setupParameters.IsoContourLevel *= dAvg;

                setupParameters.InverseDeltaX = 1.0 / (cloudSizeX / setupParameters.SampleSpaceBounds[0]);
                setupParameters.InverseDeltaY = 1.0 / (cloudSizeY / setupParameters.SampleSpaceBounds[1]);
                setupParameters.InverseDeltaZ = 1.0 / (cloudSizeZ / setupParameters.SampleSpaceBounds[2]);

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
            octree.GenerateNodeReferenceChache();
            let numCorners = octree.GenerateCornerPointsChache();
            println!("Number of corners: {}", numCorners);

            // 0, 0, 0 should be 52.6968; this error is in the algerithms--the c++ code does the same
            let inputPos = (7.970978545547169, 31.27784731879688, 7.970978545547169);
            let leafNodeIndex = octree.GetLeafIndex(inputPos);
            for _i in 0..1 {
                let dst = octree.NearestNeighborSearch(
                    &pointCloud, leafNodeIndex, inputPos
                );
                println!("dst: {}", dst.expect("Failed to get distance"));
            }

            //let mut octreeNeighborDistance: Vec <f64> = vec!();
            
            /*for node in octree.GetLeafNodes() {
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

            for x_ in 0..GRID_SIZE.0 {
                for y_ in 0..GRID_SIZE.1 {
                    for z_ in 0..GRID_SIZE.2 {
                        let xyz = (
                            x_ as f64 / setupParameters.InverseDeltaX + setupParameters.SampleSpaceOffset[0],
                            y_ as f64 / setupParameters.InverseDeltaX + setupParameters.SampleSpaceOffset[0],
                            z_ as f64 / setupParameters.InverseDeltaX + setupParameters.SampleSpaceOffset[0],
                        );

                        let leafNodeIndex = octree.GetLeafIndex(xyz);
                        let distance = octree.NearestNeighborSearch(&pointCloud, leafNodeIndex, xyz);
                        distanceField[x_ + (y_ + z_ * GRID_SIZE.1) * GRID_SIZE.0] =
                            distance.expect("Failed to get distance");
                    }
                }
            }
        }

        if USING_OLD_SDF {
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
                        let leafNodeIndex = octree.GetLeafIndex(position);
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
                                    sign as f64 +
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

        }
    }
}


