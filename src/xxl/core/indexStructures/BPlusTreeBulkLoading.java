/* XXL: The eXtensible and fleXible Library for data processing

Copyright (C) 2000-2011 Prof. Dr. Bernhard Seeger
                        Head of the Database Research Group
                        Department of Mathematics and Computer Science
                        University of Marburg
                        Germany

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 3 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library;  If not, see <http://www.gnu.org/licenses/>. 

    http://code.google.com/p/xxl/

*/

package xxl.core.indexStructures;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import xxl.core.collections.MapEntry;
import xxl.core.collections.containers.Container;
import xxl.core.functions.Constant;
import xxl.core.functions.Function;
import xxl.core.indexStructures.BPlusTree.IndexEntry;
import xxl.core.indexStructures.BPlusTree.Node;
import xxl.core.predicates.Predicate;

/**
 * This class provides functionality to bulk-load a BPlusTree.
 * The tree is created bottom-up.  
 */
public class BPlusTreeBulkLoading {
    /**
     * Bulkloading
     *  
     */
    public static final boolean ASCENDING_BULKLOAD = false;
    public static final boolean DESCENDING_BULKLOAD = true;
    
    protected BPlusTree btree; 
    protected Function<BPlusTree.Node,Container> determineTreeContainer;
    protected Predicate treeOverflows;
    protected ArrayList<MapEntry<Object,BPlusTree.Node>> treePath = null; //main buffer
    protected ArrayList<MapEntry<Object,BPlusTree.Node>> bufferPath = null; //buffer for sibling nodes
    protected boolean descending;

	/**
	 * Added by LX
	 */
	protected boolean soundRemedy = false;
	protected int remedySeed;
	protected int indexSize;
	protected ArrayList<Double> cdf = new ArrayList<>();
	protected int maxCap;
	protected int minCap;
	protected double initCap = 1.0;

	/**
	 * Bulk loads the given <tt>tree</tt>
	* with the given <tt>objects</tt> in ascending order. 
    * @param tree 
    * @param objects
    */
    public BPlusTreeBulkLoading(BPlusTree tree, Iterator objects){
    	 this(tree, objects, tree.getContainer);
    }

	/**
	 * Added by LX
	 */
	public BPlusTreeBulkLoading(BPlusTree tree, Iterator objects, int s, int i)
	{
		this(tree, objects, tree.getContainer);
	}
	/**
	 * Bulk loads the given <tt>tree</tt>
	 * with the given <tt>objects</tt> in ascending order.. 
	 * 
	 * @see BPlusTree#bulkLoad(BPlusTree, Iterator , Function&lt;BPlusTree.Node,Container&gt; , Predicate )
     * @param tree
     * @param objects
     * @param determineContainer
     */
    public BPlusTreeBulkLoading(BPlusTree tree, Iterator objects, Function<BPlusTree.Node,Container> determineContainer){
    	this(tree,  objects,  determineContainer,tree.overflows , BPlusTreeBulkLoading.ASCENDING_BULKLOAD);
    }
    
    /**
     * Bulk loads the given <tt>tree</tt>
 	* with the given <tt>objects</tt>. 
     * @param tree 
     * @param objects
     * @param order in which the data objects are sorted
     */
     public BPlusTreeBulkLoading(BPlusTree tree, Iterator objects, boolean order){
     	 this( tree, objects, tree.getContainer, order);
     }
     /**
      * Bulk loads the given <tt>tree</tt>
 	 * with the given <tt>objects</tt>. 
 	 * 
 	 * @see BPlusTree#bulkLoad(BPlusTree, Iterator , Function&lt;BPlusTree.Node,Container&gt; , Predicate )
      * @param tree
      * @param objects
      * @param determineContainer
      * @param order in which the data objects are sorted
      */
     public BPlusTreeBulkLoading(BPlusTree tree, Iterator objects, Function<BPlusTree.Node,Container> determineContainer, boolean order){
     	this(tree,  objects,  determineContainer, tree.overflows , order);
     }
     /**
      *  Bulk loads the given <tt>tree</tt>
      * with the given <tt>objects</tt> in ascending order 
      * 
      * @param tree
      * @param objects
      * @param determineContainer
      * @param overflows
      */
     public BPlusTreeBulkLoading(BPlusTree tree, Iterator objects, Function<BPlusTree.Node,Container> determineContainer, Predicate overflows){
    	 this(tree, objects,  determineContainer, BPlusTreeBulkLoading.ASCENDING_BULKLOAD  );
     }

     private void initializaCDF() {
     	// This assumes the min capacity is 0.5 and the max capacity is 1
		 maxCap = btree.singleLeafNodeCapacity;
		 minCap = (int) Math.ceil((maxCap+1) / 2);

		 cdf.add(0.0);
		 double prop, cumProp;
		 for (int i = minCap; i <= maxCap; i++) {
		 	prop = 1.0 * maxCap / ((i)*(i+1));
		 	cumProp = cdf.get(cdf.size() - 1);
		 	cdf.add(prop + cumProp);
		 }
	 }

	 private int getOneLoadSize(Random random) {
		 double p = random.nextDouble();
		 int j;
		 for (j = 1; j < cdf.size(); j++) {
			 if (p < cdf.get(j))
				 break;;
		 }
		 if (j >= cdf.size())
			 j = cdf.size() - 1;
		 return j-1+minCap;
	 }
    /**
     * Bulk loads the given <tt>tree</tt>
	 * with the given <tt>objects</tt>. 
     * Ensures that the right or left flank (depending on order)  does not contain less than D items
     * @see BPlusTree#adjustRightPath(Separator, int) 
     * @see BPlusTree#adjustLeftPath(Separator, int) 
     * @param tree
     * @param objects
     * @param determineContainer
     * @param overflows
     * @param order in which the data objects are sorted
     */
    public BPlusTreeBulkLoading(BPlusTree tree, Iterator objects, Function<BPlusTree.Node,Container> determineContainer, Predicate overflows, boolean order){
    	btree = tree;
    	this.determineTreeContainer =  determineContainer;
    	this.treeOverflows = overflows;
    	this.treePath = new ArrayList<MapEntry<Object,BPlusTree.Node>>();
    	this.bufferPath = new ArrayList<MapEntry<Object,BPlusTree.Node>>();
    	boolean isDuplicateEnabled = btree.duplicate;
    	Object first = null, last = null;
		
    	// new code 
    	this.descending = order;
    	
    	// insert all objects
		while (objects.hasNext()) {
			Object obj = objects.next();
			if (first == null) first = obj;
			last = obj;
			/**
			 * Added by LX
			 */
			btree.numOfDataInLeafNodes += 1;

			insertBulk(obj, 0);
		}

		// at this point we have a getpath which contains the right flank.
		// write unsaved nodes
		for (int level = 0; level < treePath.size();) {
			BPlusTree.IndexEntry indexEntry = saveBulk(treePath.get(level).getKey(), 
					treePath.get(level).getValue(),  isDuplicateEnabled);

			level++;
			
			if (level < treePath.size())
				insertBulk(indexEntry, level);
			else {
				tree.rootEntry = indexEntry;
				tree.rootDescriptor = tree.createKeyRange(tree.key(first), tree.key(last));
			}
		}
		adjustFlankPath(null, 0);
		
    }

	/**
	 * Added by LX
	 * Practical remedy: linear or random
	 */
	public BPlusTreeBulkLoading(BPlusTree tree, Iterator objects, Function<BPlusTree.Node,Container> determineContainer,
								Predicate overflows, boolean order, double initCap,
								int seed, boolean isLinear) throws FileNotFoundException {
		btree = tree;
		this.determineTreeContainer =  determineContainer;
		this.treeOverflows = overflows;
		this.treePath = new ArrayList<MapEntry<Object,BPlusTree.Node>>();
		this.bufferPath = new ArrayList<MapEntry<Object,BPlusTree.Node>>();
		boolean isDuplicateEnabled = btree.duplicate;
		Object first = null, last = null;

		// new code
		this.descending = order;

		/** Practical remedy */
		Random rand = new Random(seed);

		int counter;

		double maxUtil = 1.0;
		double minUtil = maxUtil - (maxUtil - initCap) * 2;
		double oddUtil = maxUtil, evenUtil = minUtil;
		int curNode = 1;
		String s;
		if (isLinear)
			s = "linear";
		else
			s = "random";

		while (objects.hasNext()) {
			int neededObj;
			if (isLinear) {
				if (curNode % 2 == 1) {
					neededObj = (int) Math.floor(btree.singleLeafNodeCapacity * oddUtil);
					oddUtil -= 0.01;
					oddUtil = oddUtil < minUtil ? maxUtil : oddUtil;
				}
				else {
					neededObj = (int) Math.floor(btree.singleLeafNodeCapacity * evenUtil);
					evenUtil += 0.01;
					evenUtil = evenUtil > maxUtil ? minUtil : evenUtil;
				}
				curNode++;
			}
			else {
				neededObj = (int) ((rand.nextDouble() * (maxUtil - minUtil) + minUtil) * btree.singleLeafNodeCapacity);
			}

			counter = 0;
			ArrayList<Object> bulk = new ArrayList<>();
			while (objects.hasNext()) {
				Object obj = objects.next();
				if (first == null) first = obj;
				last = obj;
				counter++;
				bulk.add(obj);
				if (counter >= neededObj) {
					break;
				}
			}
			insertLeaf(bulk);
			btree.numOfDataInLeafNodes += counter;
			btree.leafNodeNum++;
			btree.totalLeafNodeCapacity += btree.singleLeafNodeCapacity;
		}

		// at this point we have a getpath which contains the right flank.
		// write unsaved nodes
		for (int level = 0; level < treePath.size();) {
			BPlusTree.IndexEntry indexEntry = saveBulk(treePath.get(level).getKey(),
					treePath.get(level).getValue(),  isDuplicateEnabled);

			level++;

			if (level < treePath.size())
				insertBulk(indexEntry, level);
			else {
				tree.rootEntry = indexEntry;
				tree.rootDescriptor = tree.createKeyRange(tree.key(first), tree.key(last));
			}
		}
		adjustFlankPath(null, 0);
	}

	/**
	 * Added by LX
	 */
	public BPlusTreeBulkLoading(BPlusTree tree, Iterator objects, Function<BPlusTree.Node,Container> determineContainer, Predicate overflows, boolean order, double initCap){
		btree = tree;
		this.initCap = initCap;
		this.determineTreeContainer =  determineContainer;
		this.treeOverflows = overflows;
		this.treePath = new ArrayList<MapEntry<Object,BPlusTree.Node>>();
		this.bufferPath = new ArrayList<MapEntry<Object,BPlusTree.Node>>();
		boolean isDuplicateEnabled = btree.duplicate;
		Object first = null, last = null;

		// new code
		this.descending = order;

		// insert all objects
		while (objects.hasNext()) {
			Object obj = objects.next();
			if (first == null) first = obj;
			last = obj;
			/**
			 * Added by LX
			 */
			btree.numOfDataInLeafNodes += 1;

			insertBulk(obj, 0);
		}

		// Added by LX
//		int i = (int) Math.ceil(10000000 / (btree.singleNodeCapacity * initCap));
//		int l = (int) Math.ceil(80000000 / (btree.singleNodeCapacity * 0.7)); // assume reaches steady-state
//		btree.splitPerBatch = (int) Math.ceil((l - i) / 7000.0); // assume 7000 batches
//		btree.splitPerBatch = 20;
//		System.out.println("Averaged split per batch: " + btree.splitPerBatch);

		// at this point we have a getpath which contains the right flank.
		// write unsaved nodes
		for (int level = 0; level < treePath.size();) {
			BPlusTree.IndexEntry indexEntry = saveBulk(treePath.get(level).getKey(),
					treePath.get(level).getValue(),  isDuplicateEnabled);

			level++;

			if (level < treePath.size())
				insertBulk(indexEntry, level);
			else {
				tree.rootEntry = indexEntry;
				tree.rootDescriptor = tree.createKeyRange(tree.key(first), tree.key(last));
			}
		}
		adjustFlankPath(null, 0);

	}

	/**
	 * Added by LX
	 * Sound remedy
	 * @param tree
	 * @param objects
	 * @param determineContainer
	 * @param overflows
	 * @param order
	 * @param seed
	 * @param indexSize
	 */
	public BPlusTreeBulkLoading(BPlusTree tree, Iterator objects, Function<BPlusTree.Node,Container> determineContainer,
								Predicate overflows, boolean order, int seed, int indexSize) throws FileNotFoundException {
		btree = tree;
		this.determineTreeContainer =  determineContainer;
		this.treeOverflows = overflows;
		this.treePath = new ArrayList<MapEntry<Object,BPlusTree.Node>>();
		this.bufferPath = new ArrayList<MapEntry<Object,BPlusTree.Node>>();
		boolean isDuplicateEnabled = btree.duplicate;
		Object first = null, last = null;

		// new code
		this.descending = order;

		/** Sound Remedy */
		// First determine the value of j
		Random random = new Random(seed);
		initializaCDF();
		int remainingObjectSize = indexSize;
		int recordNum;
		int counter;

		boolean isLastIter = false, lastTwoPages = false;
		while (!isLastIter) {
			if (remainingObjectSize <= maxCap) {
				isLastIter = true;
			}
			else if (remainingObjectSize <= (3.0 * maxCap + 1) / 2) {
				lastTwoPages = true;
			}

			if (lastTwoPages) {
				recordNum = remainingObjectSize / 2;
			}
			else {
				recordNum = getOneLoadSize(random);
			}
			remainingObjectSize -= recordNum;
			counter = 0;
			ArrayList<Object> bulk = new ArrayList<>();
			while (objects.hasNext()) {
				Object obj = objects.next();
				if (first == null) first = obj;
				last = obj;
				counter++;
				bulk.add(obj);
				if (!isLastIter && counter >= recordNum) {
					break;
				}
			}
			insertLeaf(bulk);
			btree.numOfDataInLeafNodes += counter;
			btree.leafNodeNum++;
			btree.totalLeafNodeCapacity += btree.singleLeafNodeCapacity;
		}
		// at this point we have a getpath which contains the right flank.
		// write unsaved nodes
		for (int level = 0; level < treePath.size();) {
			BPlusTree.IndexEntry indexEntry = saveBulk(treePath.get(level).getKey(),
					treePath.get(level).getValue(),  isDuplicateEnabled);

			level++;

			if (level < treePath.size())
				insertBulk(indexEntry, level);
			else {
				tree.rootEntry = indexEntry;
				tree.rootDescriptor = tree.createKeyRange(tree.key(first), tree.key(last));
			}
		}
		adjustFlankPath(null, 0);

	}
   /**
    * Method corrects left flank while bulk loading the data objects in descending order mode
    * and right flank in ascending mode
    * @param newSep
    * @param level
    */ 
   protected void adjustFlankPath(Separator newSep, int level){
	   Container container =  btree.container();
	   Node flankNode = (Node)treePath.get(level).getValue();
	   Object flankEntry = treePath.get(level).getKey();
	   int secondIndex = 1; // index of the second entry in the current node
	   int firstIndex = 0; // index of the first entry in the current node
	   if (level < bufferPath.size()){
		   Node sibling =  (Node)bufferPath.get(level).getValue(); 
		   Object siblingEntry  = bufferPath.get(level).getKey();
		   Separator newSeparator = null; // 
		   if (level > 0  && newSep != null) {//below root node and have created new Separator with redistribute 
    			if (descending){ // left flank
    				((IndexEntry)flankNode.getFirst()).separator = (Separator)newSep.clone();//father node has more than one entry so we must update last node instead of left sibling.	
    			}else{ // right flank
        				if (flankNode.number() > 1){//more than one item in right node, therfore update
        					((IndexEntry)flankNode.getEntry(flankNode.number()-2)).separator = (Separator)newSep.clone();//father node has more than one entry so we must update last node instead of left sibling.
        				}
        				else{
        					((IndexEntry)sibling.getLast()).separator = (Separator)newSep.clone();
        				} 			
    			} // end adjust separator
    		}
		   //redistribute
			if (flankNode.underflows()){
   			int D = sibling.level() == 0 ? btree.D_LeafNode : btree.D_IndexNode; //
   			if (descending){
   				List newEntries = sibling.entries.subList(firstIndex, D); //
       			flankNode.entries.addAll(flankNode.number(), newEntries);//
       			newEntries.clear();
       			newSeparator =  (Separator)btree.separator(flankNode.getLast()).clone();
   			}else{
       			List newEntries = sibling.entries.subList(D, sibling.number()); //
       			flankNode.entries.addAll(0, newEntries);//
       			newEntries.clear();
       			newSeparator = (Separator)btree.separator(sibling.getLast()).clone();
   			}
   			//update nodes
   			container.update(flankEntry, flankNode);
   			container.update(siblingEntry, sibling);
			}// end redistribute
			//shared parent node without underflow
			if (newSeparator == null && newSep != null ){
					container.update(flankEntry, flankNode);
			}
			//Recursion
			adjustFlankPath(newSeparator, level + 1);
	   }else{
		   if ( newSep != null ){
			   if (descending){		   
	    				((IndexEntry)flankNode.getFirst()).separator = (Separator)newSep.clone();    	
			   }else{   
		   				((IndexEntry)flankNode.getEntry(flankNode.number()-2)).separator = (Separator)newSep.clone();
				  }
			   container.update(flankEntry, flankNode);
		   }
	   }
   }

	/**
	 * Added by LX
	 */
	protected void insertLeaf(ArrayList<Object> objs) {
		BPlusTree.Node node = null;
		Object id = null;
		if (treePath.size() > 0) {
			// Obtain the most recent leaf node
			node = treePath.get(0).getValue();
			id = treePath.get(0).getKey();
		}

		// Create a node to hold these new objects
		BPlusTree.Node newNode = (BPlusTree.Node) btree.createNode(0);
		Object newId = determineTreeContainer.invoke(newNode).reserve(new Constant(newNode));
		objs.forEach((n) -> insertIntoNode(newNode, n));

		if (node == null) {
			// We just created the first leaf node
			treePath.add(new MapEntry<Object, BPlusTree.Node>(newId, newNode));
		}
		else {
			// Add link to sibling leaf nodes and create their parent
			// newNode is now the most recent leaf node
			treePath.set(0, new MapEntry<Object,BPlusTree.Node>(newId, newNode));
			if (descending){
				newNode.nextNeighbor = (BPlusTree.IndexEntry)btree.createIndexEntry(1);
				newNode.nextNeighbor.initialize(id);
			}
			else{
				node.nextNeighbor = (BPlusTree.IndexEntry)btree.createIndexEntry(1);
				node.nextNeighbor.initialize(newId);
			}
			if (bufferPath.size() <= 0){
				bufferPath.add(new MapEntry<Object,BPlusTree.Node>(id, node));
			}
			else{
				bufferPath.set(0, new MapEntry<Object,BPlusTree.Node>(id, node));
			}
			insertBulk(saveBulk(id, node, btree.isDuplicatesEnabled()), 1);
		}
	}

	/**
     * Inserts an entry into the given level.
     * @param entry
     * @param level
     */
    protected  void insertBulk(Object entry, int level){
    	if (treePath.size() <= level) {

			BPlusTree.Node newNode = (BPlusTree.Node)btree.createNode(level);
			Object newId = determineTreeContainer.invoke(newNode).reserve(new Constant(newNode));

			treePath.add(new MapEntry<Object,BPlusTree.Node>(newId, newNode));
			/** Added by LX */
			if (level == 0) {
				btree.leafNodeNum += 1;
				btree.totalLeafNodeCapacity += btree.singleLeafNodeCapacity;
				newNode.nodeId =  btree.leafNodeNum;
			}
		}

		BPlusTree.Node node = treePath.get(level).getValue();
		Object id = treePath.get(level).getKey();

		insertIntoNode(node, entry );
		/** Added by LX */
		if (node.level == 0 && !btree.leafNodeToRange.containsKey(node.nodeId)) {
			btree.leafNodeToRange.put(node.nodeId, btree.getDescriptor.invoke(entry));
		}

		// Modified by LX:
		// To accomodate for the initial node capacity
		if (treeOverflows.invoke(node) || (node.level == 0 && node.number() >= initCap*btree.singleLeafNodeCapacity) ||
				(node.level> 0 && node.number() >= initCap * btree.singleIndexNodeCapacity)) {


			BPlusTree.Node newNode = (BPlusTree.Node)btree.createNode(node.level);
			Object newId = determineTreeContainer.invoke(newNode).reserve(new Constant(newNode));
			insertIntoNode(newNode, entry);
			treePath.set(level, new MapEntry<Object,BPlusTree.Node>(newId, newNode));

			Iterator entries;
			for (entries = node.entries(); !entries.next().equals(entry); );
			entries.remove();
			if (descending){
				newNode.nextNeighbor = (BPlusTree.IndexEntry)btree.createIndexEntry(level+1);
				newNode.nextNeighbor.initialize(id);
			}
			else{
				node.nextNeighbor = (BPlusTree.IndexEntry)btree.createIndexEntry(level+1);
				node.nextNeighbor.initialize(newId);
			}
			if (bufferPath.size() <= level){
				bufferPath.add(new MapEntry<Object,BPlusTree.Node>(id, node)); 	
			}
			else{
				bufferPath.set(level, new MapEntry<Object,BPlusTree.Node>(id, node)); 	
			}

			/** Added by LX */
			int idx = (int) (initCap * 100 - 50);
			if (node.level == 0) {
//				btree.leafNodeToRange.put(node.nodeId, ((HilbertRTree)btree).computeMBR(node.entries()));
				btree.leafNodeNum += 1;
				btree.totalLeafNodeCapacity += btree.singleLeafNodeCapacity;
				newNode.nodeId = btree.leafNodeNum;
				btree.leafNodeToRange.put(newNode.nodeId, btree.getDescriptor.invoke(entry));
			}

			insertBulk(saveBulk(id, node, btree.isDuplicatesEnabled()), level+1);
//			System.out.println(node.number() + "," + newNode.number());
		}
    }
    
    /**
     * Inserts entry depending on the order 
     * @param node
     * @param entry
     */
    protected  void insertIntoNode(BPlusTree.Node node, Object entry ){
    	if (descending){
			node.entries.add(0,entry);
		}
		else{
			node.entries.add(entry);
		}
    }
    /**
     * Saves a node of the tree to external memory.
     * @param id
     * @param node
     * @param isDuplicateEnabled
     * @return
     */
    protected  BPlusTree.IndexEntry saveBulk (Object id, BPlusTree.Node node, boolean isDuplicateEnabled){
    	Container container = determineTreeContainer.invoke(node);
		container.update(id, node);
		Separator sep = (Separator) btree.separator(node.getLast()).clone();
		return (BPlusTree.IndexEntry)((BPlusTree.IndexEntry)btree.createIndexEntry(node.level+1)).initialize(sep).initialize(container, id);
    }
	
}
