#![allow(non_snake_case)]

// some of the methods aren't being used,
// but are implimented as they may be usefull
// in future projects
#![allow(dead_code)]


enum FutureChildType {
    NextLayer,
    Incomplete,
    None,
}


struct Node {
    leftChild: Option <usize>,
    rightChild: Option <usize>,
    parent: Option <usize>,
    incompleteNodeLeft: Option <usize>,
    incompleteNodeRight: Option <usize>,
    typeOfFutureChildren: FutureChildType,
    depth: usize,
    // just swap this to reduce memory movement
    value: (f64, usize),
    changedL: bool,  // the left/right and if they moved forward
}


pub struct MinHeapBinaryTree {
    // stores the references to each child
    childReferences: Vec <Node>,

    // root node is always the first index (0)

    // this contains all incomplete nodes
    // these are the first place for a new value to be inserted
    incompleteNodes: Vec <(Option <usize>, usize)>,
    nextLayer: Vec <(Option <usize>, usize)>,  // the next layer of nodes (so the depth can expand)
    numberOfNodes: usize,

    // represents the deepest node
    // when a leaf node is inserted, if it's
    // shallower than this, push it to the start of the vector
    // and also push incomplete nodes
    // this should help balance the tree?
    deepest: usize,
}


impl MinHeapBinaryTree {
    pub fn new () -> Self {
        MinHeapBinaryTree {
            childReferences: vec!(),  // no references
            incompleteNodes: vec![(None, 0)],  // the root node
            nextLayer: vec!(),
            numberOfNodes: 0,
            deepest: 0,  // the root node
        }
    }

    pub fn Push (&mut self, value: (f64, usize)) {
        // gather the ne index for the node
        if let Some((parent, depth)) = self.incompleteNodes.pop() {
            if parent.is_none() && self.incompleteNodes.len() > 0 {
                // it wasn't the root node but rather an invalidated one; continuing the search for a valid incomplete node
                self.Push(value);
                return;
            }

            // getting the node
            let mut newNode = Node {
                leftChild: None,
                rightChild: None,
                parent,
                incompleteNodeLeft: None,
                incompleteNodeRight: None,
                typeOfFutureChildren: FutureChildType::Incomplete,
                depth,
                value,
                changedL: false,
            };

            // pushing the new future children
            if depth >= self.deepest {
                newNode.typeOfFutureChildren = FutureChildType::NextLayer;
                newNode.incompleteNodeLeft = Some(self.nextLayer.len());
                newNode.incompleteNodeRight = Some(self.nextLayer.len() + 1);
                self.nextLayer.push((Some(self.numberOfNodes), depth + 1));
                self.nextLayer.push((Some(self.numberOfNodes), depth + 1));
            } else {
                newNode.incompleteNodeLeft = Some(self.incompleteNodes.len());
                newNode.incompleteNodeRight = Some(self.incompleteNodes.len() + 1);
                self.incompleteNodes.push((Some(self.numberOfNodes), depth + 1));
                self.incompleteNodes.push((Some(self.numberOfNodes), depth + 1));
            }
            self.childReferences.push(newNode);
            self.numberOfNodes += 1;  // one node was added

            if parent.is_none() { return; }
            
            // updating the children and incomplete nodes
            let mut validParent = parent.unwrap();
            self.childReferences[validParent].typeOfFutureChildren = FutureChildType::None;
            if self.childReferences[validParent].leftChild.is_none() {
                self.childReferences[validParent].leftChild = Some(self.numberOfNodes - 1);
                self.childReferences[validParent].incompleteNodeLeft = None;
                self.childReferences[validParent].changedL = false;  // resetting
            } else {
                self.childReferences[validParent].rightChild = Some(self.numberOfNodes - 1);
                self.childReferences[validParent].incompleteNodeRight = None;
            }

            // swapping until the value satisfies the min-heap rules
            let mut parent: Option <usize>;

            let mut currentValue: (f64, usize);
            let mut currentIndex = self.numberOfNodes - 1;
            loop {
                currentValue = self.childReferences[currentIndex].value;
                parent = self.childReferences[currentIndex].parent;
                if parent.is_none() { return; }
                validParent = parent.unwrap();

                if value.0 < self.childReferences[validParent].value.0 && validParent != currentIndex {
                    // swapping the values and continuing
                    self.childReferences[currentIndex].value = self.childReferences[validParent].value;
                    self.childReferences[validParent].value = currentValue;
                    currentIndex = validParent;
                } else {
                    return;  // the min-heap rules are satisfied
                }
            }
        }

        // it could be that it was an invalid node
        if self.incompleteNodes.is_empty() {
            // moving all nodes from the next depth to the incompelte buffer
            while let Some(node) = self.nextLayer.pop() {
                if let Some(index) = node.0 {  // index == parent?
                    if let Some(_parent) = self.childReferences.get(index) {
                        if self.childReferences[index].incompleteNodeLeft == Some(self.nextLayer.len()) &&
                            !self.childReferences[index].changedL {
                            
                            self.childReferences[index].incompleteNodeLeft = Some(self.incompleteNodes.len());
                            self.childReferences[index].changedL = true;
                        } else {
                            self.childReferences[index].incompleteNodeRight = Some(self.incompleteNodes.len());
                        }
                        self.childReferences[index].typeOfFutureChildren = FutureChildType::Incomplete;
                        self.incompleteNodes.push((Some(index), self.incompleteNodes.len()));
                    }
                }
            }
            self.deepest += 1;  // the next layer has been reached
        }
        
        // calling insert to actually insert a value
        self.Push(value);
    }

    pub fn Pop (&mut self) -> Option <(f64, usize)> {
        // find a leaf node
        // take that leaf node, pop it off,
        // place that leaf node at the root replacing the current root
        // swap the root node down until it satisfies the rules

        // getting the leaf node, and then finding the next leaf node
        // does this work, or is this completely flawed??   (from basic testing it seems true? none of them ever had a child)
        let initialValue = Some(self.childReferences[0].value);
        if let Some(node) = self.childReferences.pop() {
            self.numberOfNodes -= 1;
            // updating the node's parent to not reference it and removing any
            // potential new positions from the vector of them
            if let Some(parent) = node.parent {  // 0 is the root, which would fail this check so it's safe to use as an exception case
                if self.childReferences[parent].leftChild == Some(self.childReferences.len()) {
                    self.childReferences[parent].leftChild = None;
                    self.childReferences[parent].incompleteNodeLeft = Some(self.incompleteNodes.len());
                } else {
                    self.childReferences[parent].rightChild = None;
                    self.childReferences[parent].incompleteNodeRight = Some(self.incompleteNodes.len());
                }

                // is this correct? it should be push the parent right? not the root?     seems to work.....
                self.incompleteNodes.push((Some(parent), node.depth));
            } else {
                self.childReferences.clear();
                self.incompleteNodes = vec![(None, 0)];
                self.numberOfNodes = 0;
                return initialValue;
            }
            
            // updating the deleted nodes incomplete references
            if matches!(node.typeOfFutureChildren, FutureChildType::NextLayer) {
                if let Some(index) = node.incompleteNodeLeft {
                    self.nextLayer[index] = (None, 0);
                } if let Some(index) = node.incompleteNodeRight {
                    self.nextLayer[index] = (None, 0);
                }
            } else {
                if let Some(index) = node.incompleteNodeLeft {
                    self.incompleteNodes[index] = (None, 0);
                }
                if let Some(index) = node.incompleteNodeRight {
                    self.incompleteNodes[index] = (None, 0);
                }
            }

            self.childReferences[0].value = node.value;

            // shifting the value down until the rules are satisfied
            let mut currentValue: (f64, usize);
            let mut currentNode = 0;
            loop {
                if let Some(child) = self.childReferences[currentNode].leftChild {
                    if self.childReferences[currentNode].value.0 > self.childReferences[child].value.0 {
                        currentValue = self.childReferences[currentNode].value;
                        self.childReferences[currentNode].value = self.childReferences[child].value;
                        self.childReferences[child].value = currentValue;
                        currentNode = child;
                    } else if let Some(child) = self.childReferences[currentNode].rightChild {
                        if self.childReferences[currentNode].value.0 > self.childReferences[child].value.0 {
                            currentValue = self.childReferences[currentNode].value;
                            self.childReferences[currentNode].value = self.childReferences[child].value;
                            self.childReferences[child].value = currentValue;
                            currentNode = child;
                        } else {
                            return initialValue;
                        }
                    } else {
                        return initialValue;
                    }
                } else {
                    // this is a leaf node
                    return initialValue;
                }
            }
        }

        None  // no valid nodes to be retrieved
    }


    pub fn IsEmpty (&self) -> bool {
        self.numberOfNodes == 0
    }

    pub fn Print (&self) {
        for node in &self.childReferences {
            println!("    Node: {}    value: {}    depth: {}    parent: {}    left: {}    right: {}", self.numberOfNodes - 1, node.value.0, node.depth, node.parent.unwrap_or(0), node.leftChild.unwrap_or(0), node.rightChild.unwrap_or(0));
        }
    }

    pub fn GetMin (&self) -> (f64, usize) {
        self.childReferences[0].value
    }

}


