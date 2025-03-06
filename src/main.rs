#![allow(non_snake_case)]

/*
Octree construction should be working I think
The number of cells checks out at a minimum

Some of the other sections are also working so far
*/

mod octree;

const LOADING_DISTANCE_FIELD_SAVE: bool = false;
const GENERATE_NEW_SDF: bool = true;
const USING_OLD_SDF: bool = true;

const FIELD_FILE: &str = "";
const SAVE_FILE: &str = "";

const OBJ_SAVE_FILE: &str = "";

const LOAD_PCD: bool = false;
const PCD_FILE: &str = "";
const POINTS_LOADED_PERCENT: f64 = 100.0 / (100.0);

const GENERATE_HALLOW_SPHERE: bool = true;
const GENERATE_SOLID_SPHERE: bool = false;
const GENERATE_CUBE: bool = false;

const MAX_OCTREE_DEPTH: usize = 10;

const SHELL_CHUNK_SIZE: f64 = 10.0;

const SIMULATING_SURFACE_TENSION: bool = false;
const TENSION_ITTERATIONS: isize = 2500;
const VELOCITY: f64 = 0.75;

const GRID_SIZE: (usize, usize, usize) = (100, 100, 100);

const MAX_SAMPLING_SPACE: usize = GRID_SIZE.0 * GRID_SIZE.1 * GRID_SIZE.2;


const DT: f64 = 0.0055 * 0.1;

const AUTO_SET: bool = true;

const PI: f64 = 3.14159;



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

    let mut numNeighborCopies = [0usize; MAX_SAMPLING_SPACE];

    neighbors.push(
        (startX, startY, startZ)
    );
    let mut numNeighbors = 1usize;

    while numNeighbors > 0 {
        let mut numNewNeighbors = 0usize;

        for neighborIndex in 0..numNeighbors {
            let (indexX, indexY, indexZ) = neighbors[neighborIndex];
            let arrayIndex = indexX + (indexY + indexZ * GRID_SIZE.1) * GRID_SIZE.0;
            if signedGrid[arrayIndex] != 0 {
                if distanceField[arrayIndex] >= setupParameters.IsoContourLevel{
                    signedGrid[arrayIndex] = sign;

                    for i in 0..6 {
                        let (offsetX, offsetY, offsetZ) = FLOOD_FILL_NEIGHBOR_OFFSETS[i];
                        let neighborXPos = (indexX as isize + offsetX).unsigned_abs();
                        let neighborYPos = (indexX as isize + offsetY).unsigned_abs();
                        let neighborZPos = (indexX as isize + offsetZ).unsigned_abs();

                        if neighborXPos < GRID_SIZE.0 &&
                           neighborYPos < GRID_SIZE.1 &&
                           neighborYPos < GRID_SIZE.2 {

                            let indexKey = neighborXPos + (neighborYPos + neighborZPos * GRID_SIZE.1) * GRID_SIZE.0;
                            if numNeighborCopies[indexKey] == 0 {
                                numNewNeighbors += 1;
                                numNeighborCopies[indexKey] += 1;
                                newNeighbors.push(
                                    (neighborXPos, neighborYPos, neighborZPos)
                                );
                            }
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
                if signedGrid[index] >= 1 {
                    sign = signedGrid[index];
                } else{
                    if distanceField[index] >= setupParameters.IsoContourLevel {
                        sign *= -1;

                        FloodFill(setupParameters, (x, y, z), distanceField, signedGrid, sign);
                    }
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
    
    let phi = PI * (0.70710678118 - 1.0);  // the magic number is sqrt 0.5

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
    SampleSpaceSize: [f64; 3],
    SampleSpaceOffset: [f64; 3],
    IsoContourLevel: f64,
    InverseDeltaX: f64,
    InverseDeltaY: f64,
    InverseDeltaZ: f64,
}


fn main() {
    let mut setupParameters = SetupParameters {
        SampleSpaceSize: [101.0, 101.0, 101.0],
        SampleSpaceOffset: [0.0, 0.0, 0.0],
        IsoContourLevel: 6.0,
        InverseDeltaX: 1.0 / (1.0),
        InverseDeltaY: 1.0 / (1.0),
        InverseDeltaZ: 1.0 / (1.0),
    };

    let mut distanceField = [0.0; MAX_SAMPLING_SPACE];

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
                    (std::f64::MAX, std::f64::MAX, std::f64::MAX);
                let (mut largestX, mut largestY, mut largestZ): (f64, f64, f64) =
                    (std::f64::MIN, std::f64::MIN, std::f64::MIN);

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

                setupParameters.InverseDeltaX = cloudSizeX / setupParameters.SampleSpaceSize[0];
                setupParameters.InverseDeltaY = cloudSizeY / setupParameters.SampleSpaceSize[1];
                setupParameters.InverseDeltaZ = cloudSizeZ / setupParameters.SampleSpaceSize[2];

                setupParameters.IsoContourLevel = (
                    setupParameters.InverseDeltaX + setupParameters.InverseDeltaY + setupParameters.InverseDeltaZ
                ) / 3.0 * setupParameters.IsoContourLevel;

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

                setupParameters.InverseDeltaX = 1.0 / (cloudSizeX / setupParameters.SampleSpaceSize[0]);
                setupParameters.InverseDeltaY = 1.0 / (cloudSizeY / setupParameters.SampleSpaceSize[1]);
                setupParameters.InverseDeltaZ = 1.0 / (cloudSizeZ / setupParameters.SampleSpaceSize[2]);

                setupParameters.SampleSpaceSize = [
                    cloudSizeX, cloudSizeY, cloudSizeZ
                ];

                println!("Cloud size: {}, {}, {}", setupParameters.SampleSpaceSize[0], setupParameters.SampleSpaceSize[1], setupParameters.SampleSpaceSize[2]);
                println!("Cloud offset: {}, {}, {}", setupParameters.SampleSpaceOffset[0], setupParameters.SampleSpaceOffset[1], setupParameters.SampleSpaceOffset[2]);
                println!("Deltas: {}, {}, {}", 1.0/setupParameters.InverseDeltaX, 1.0/setupParameters.InverseDeltaY, 1.0/setupParameters.InverseDeltaZ);
                println!("Lowest: {}, {}, {}", lowestX, lowestY, lowestZ);
                println!("largest: {}, {}, {}", largestX, largestY, largestZ);
                println!("Iso contour level: {}", setupParameters.IsoContourLevel);
            }

        let mut octree = octree::Octree::new(
            (setupParameters.SampleSpaceOffset[0], setupParameters.SampleSpaceOffset[1], setupParameters.SampleSpaceOffset[2]),
            (setupParameters.SampleSpaceSize[0], setupParameters.SampleSpaceSize[1], setupParameters.SampleSpaceSize[2]),
            MAX_OCTREE_DEPTH
        );
        octree.SubDivide(&pointCloud);
        octree.GenerateNodeReferenceChache();

        //
        }
    }

}


