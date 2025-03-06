#![allow(non_snake_case)]


struct Node {
    leftChild: Option <usize>,
    rightChild: Option <usize>,
    parent: Option <usize>,
    depth: usize,
    // just swap this to reduce memory movement
    value: f64,
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

    pub fn GetMin (&self) -> f64 {
        return self.childReferences[0].value;
    }

    pub fn Insert (&mut self, value: f64) {
        // gather the ne index for the node
        if let Some((parent, depth)) = self.incompleteNodes.pop() {
            // getting the node
            let newNode = Node {
                leftChild: None,
                rightChild: None,
                parent: parent,
                depth: depth,
                value: value
            };

            // pushing the new future children
            if depth >= self.deepest {
                self.nextLayer.push((Some(self.numberOfNodes), depth + 1));
                self.nextLayer.push((Some(self.numberOfNodes), depth + 1));
            } else {
                self.incompleteNodes.push((Some(self.numberOfNodes), depth + 1));
                self.incompleteNodes.push((Some(self.numberOfNodes), depth + 1));
            }
            self.childReferences.push(newNode);
            self.numberOfNodes += 1;  // two nodes were added

            if !parent.is_none() {
                let mut validParent = parent.unwrap();
                if self.childReferences[validParent].leftChild.is_none() {
                    self.childReferences[validParent].rightChild = Some(self.numberOfNodes);
                } else {
                    self.childReferences[validParent].leftChild = Some(self.numberOfNodes);
                }

                // swapping until the value satisfies the min-heap rules
                let mut parent: Option <usize>;

                let mut currentValue: f64;
                let mut currentIndex = self.numberOfNodes - 1;
                loop {
                    currentValue = self.childReferences[currentIndex].value;
                    parent = self.childReferences[currentIndex].parent;
                    if parent.is_none() { return; }
                    validParent = parent.unwrap();

                    if value < self.childReferences[validParent].value && validParent != currentIndex {
                        // swapping the values and continuing
                        self.childReferences[currentIndex].value = self.childReferences[validParent].value;
                        self.childReferences[validParent].value = currentValue;
                        currentIndex = validParent;
                    } else {

                        /*for node in &self.childReferences {
                            println!("    Node: {}    value: {}    depth: {}    parent: {}    left: {}    right: {}", self.numberOfNodes - 1, node.value, node.depth, node.parent.unwrap_or(0), node.leftChild.unwrap_or(0), node.rightChild.unwrap_or(0));
                        }*/

                        return;  // the min-heap rules are satisfied
                    }
                }
            }
        }

        // moving all nodes from the next depth to the incompelte buffer
        while let Some(nodeIndex) = self.nextLayer.pop() {
            self.incompleteNodes.push(nodeIndex);
        }
        self.deepest += 1;  // the next layer has been reached
        // calling insert to actually insert a value
        self.Insert(value);
    }

    pub fn Pop (&mut self) -> f64 {
        0.0
    }
}

