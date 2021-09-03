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
import java.util.Random;

import xxl.core.collections.containers.Container;
import xxl.core.functions.Function;
import xxl.core.predicates.Predicate;

/** This class provides functionality to bulk-load different types of trees.
 * The tree is created buttom-up.  
 */
public class SortBasedBulkLoading {

	/** The tree which is bulk-loaded.
	 */
	protected ORTree tree;
	
	/** A Function determining into which container a node is saved.
	 */
	protected Function determineContainer;	// Node -> Container
	
	/** A predicate determining if a node overflows (and i.e. a new node
	 * has to be allocated).
	 */
	protected Predicate overflows;
	
	/** An <tt>ArrayList</tt> containing the path of nodes in the tree where 
	 * insertions take place, i.e. the right edge of the tree.
	 */
	protected ArrayList path = new ArrayList();
	
	/** A flag indicating if <tt>tree</tt> is an <tt>MTree</tt>.
	 */
	protected boolean isMTree = false;

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
	protected boolean isRRstarTree = false; // Need to set the middle point of the node

	/** Creates a new <tt>SortBasedBulkLoading</tt> and bulk loads the given <tt>tree</tt>
	 * with the given <tt>objects</tt>.
	 *
	 * @param tree the empty <tt>ORtree</tt> to bulk-load
	 * @param objects sorted iterator of the objects to be loaded into the tree
	 * @param determineContainer a function determining into which container nodes are saved
	 * @param overflows a predicate determining when a node overflows
	 */
	public SortBasedBulkLoading (ORTree tree, Iterator objects, Function determineContainer, Predicate overflows) {
		this.tree = tree;
		this.determineContainer = determineContainer;
		this.overflows = overflows;
		this.isMTree = tree instanceof MTree;
		this.isRRstarTree = tree instanceof RevisedRTree; // LX

		// insert all objects
		if (isMTree) {
			while (objects.hasNext()) {
				insert(((MTree)tree).new LeafEntry(objects.next()), 0);
			}
		}
		else
			while (objects.hasNext()) {
				insert(objects.next(), 0);
				/**
				 * Added by LX
				 */
				tree.numOfDataInLeafNodes += 1;
			}

		// write unsaved nodes
		for (int level = 0; level<path.size();) {
			ORTree.IndexEntry indexEntry = save((ORTree.Node)path.get(level++));

			if (level<path.size())
				insert(indexEntry, level);
			else {
				tree.rootEntry = indexEntry;
				tree.rootDescriptor = indexEntry.descriptor();
			}
		}
	}

	/**
	 * Added by LX
	 */
	public SortBasedBulkLoading (ORTree tree, Iterator objects, Function determineContainer, Predicate overflows, double initCap) {
		this.tree = tree;
		this.initCap = initCap;
		this.determineContainer = determineContainer;
		this.overflows = overflows;
		this.isMTree = tree instanceof MTree;
		this.isRRstarTree = tree instanceof RevisedRTree; // LX

		// insert all objects
		if (isMTree) {
			while (objects.hasNext()) {
				insert(((MTree)tree).new LeafEntry(objects.next()), 0);
			}
		}
		else
			while (objects.hasNext()) {
				Object obj = objects.next();
				insert(obj, 0);
				/**
				 * Added by LX
				 */
				tree.numOfDataInLeafNodes += 1;
			}

		// write unsaved nodes
		for (int level = 0; level<path.size();) {
			ORTree.IndexEntry indexEntry = save((ORTree.Node)path.get(level++));

			if (level<path.size())
				insert(indexEntry, level);
			else {
				tree.rootEntry = indexEntry;
				tree.rootDescriptor = indexEntry.descriptor();
			}
		}
	}

	private void initializaCDF() {
		// This assumes the min capacity is 0.5 and the max capacity is 1
		maxCap = tree.singleLeafNodeCapacity;
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
	 * Added by LX
	 * Practical remedy: linear or random
	 */
	public SortBasedBulkLoading (ORTree tree, Iterator objects, Function determineContainer, Predicate overflows, int seed,
								 double initCap, String outDir, boolean isLinear) throws FileNotFoundException {

		this.tree = tree;
		this.determineContainer = determineContainer;
		this.overflows = overflows;
		this.isMTree = tree instanceof MTree;

		/** Practical remedy */
		// First determine the value of j
		Random rand = new Random(seed);
//		initializaCDF();
		int remainingObjectSize = indexSize;
		int recordNum;

		/**
		 * LX: add implementation of practical remedy
		 */
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
					neededObj = (int) Math.floor(tree.singleLeafNodeCapacity * oddUtil);
					oddUtil -= 0.01;
					oddUtil = oddUtil < minUtil ? maxUtil : oddUtil;
				}
				else {
					neededObj = (int) Math.floor(tree.singleLeafNodeCapacity * evenUtil);
					evenUtil += 0.01;
					evenUtil = evenUtil > maxUtil ? minUtil : evenUtil;
				}
				curNode++;
			}
			else {
				neededObj = (int) ((rand.nextDouble() * (maxUtil - minUtil) + minUtil) * tree.singleLeafNodeCapacity);
			}

			counter = 0;
			ArrayList<Object> bulk = new ArrayList<>();
			while (objects.hasNext()) {
				Object obj = objects.next();
				counter++;
				bulk.add(obj);
				if (counter >= neededObj) {
					break;
				}
			}
			insertLeaf(bulk);
			tree.numOfDataInLeafNodes += counter;
			tree.leafNodeNum++;
			tree.totalLeafNodeCapacity += tree.singleLeafNodeCapacity;
		}

//		while (objects.hasNext()) {
//			insert(objects.next(), 0);
//
//			tree.numOfDataInLeafNodes += 1;
//		}

		// write unsaved nodes
		for (int level = 0; level<path.size();) {
			ORTree.IndexEntry indexEntry = save((ORTree.Node)path.get(level++));

			if (level<path.size())
				insert(indexEntry, level);
			else {
				tree.rootEntry = indexEntry;
				tree.rootDescriptor = indexEntry.descriptor();
			}
		}
	}

	/**
	 * Added by LX
	 */
	public SortBasedBulkLoading (ORTree tree, Iterator objects, Function determineContainer, Predicate overflows, int seed, int indexSize, String outDir) throws FileNotFoundException {
		this.tree = tree;
		this.determineContainer = determineContainer;
		this.overflows = overflows;
		this.isMTree = tree instanceof MTree;

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
				counter++;
				bulk.add(obj);
				if (!isLastIter && counter >= recordNum) {
					break;
				}
			}
			insertLeaf(bulk);
			tree.numOfDataInLeafNodes += counter;
			tree.leafNodeNum++;
			tree.totalLeafNodeCapacity += tree.singleLeafNodeCapacity;

		}

		// write unsaved nodes
		for (int level = 0; level<path.size();) {
			ORTree.IndexEntry indexEntry = save((ORTree.Node)path.get(level++));

			if (level<path.size())
				insert(indexEntry, level);
			else {
				tree.rootEntry = indexEntry;
				tree.rootDescriptor = indexEntry.descriptor();
			}
		}
	}


	/** Creates a new <tt>SortBasedBulkLoading</tt> and bulk loads the given <tt>tree</tt>
	 * with the given <tt>objects</tt>. The predicate overflows of the tree is used to
	 * determine when a node overflows.
	 *
	 * @param tree the empty <tt>ORtree</tt> to bulk-load
	 * @param objects sorted iterator of the objects to be loaded into the tree
	 * @param determineContainer a function determining into which container nodes are saved
	 */
	public SortBasedBulkLoading (ORTree tree, Iterator objects, Function determineContainer) {
		this(tree, objects, determineContainer, tree.overflows);
	}

	/**
	 * Added by LX
	 */
	public SortBasedBulkLoading (ORTree tree, Iterator objects, Function determineContainer, double initCap) {
		this(tree, objects, determineContainer, tree.overflows, initCap);
	}

	/**
	 * Added by LX
	 */
	public SortBasedBulkLoading (ORTree tree, Iterator objects, Function determineContainer, int seed, int indexSize, String outDir) throws FileNotFoundException {
		this(tree, objects, determineContainer, tree.overflows, seed, indexSize, outDir);
	}

	/**
	 * Added by LX
	 */
	public SortBasedBulkLoading (ORTree tree, Iterator objects, Function determineContainer, int seed, double initCap, String outDir, boolean isLinear) throws FileNotFoundException {
		this(tree, objects, determineContainer, tree.overflows, seed, initCap, outDir, isLinear);
	}

	/**
	 * Added by LX
	 */
	protected void insertLeaf(ArrayList<Object> objs) {
		// Create a node to hold these new objects
		Tree.Node newNode = tree.createNode(0);
		objs.forEach((n) -> ((ORTree.Node) newNode).entries.add(n));

		ORTree.Node node = null;
		if (path.size() > 0) {
			// Obtain the most recent leaf node
			node = (ORTree.Node)path.get(0);
			path.set(0, newNode);
			insert(save(node), 1);
		}
		else {
			path.add(newNode);
		}
	}

	/** Inserts an entry into the given level.
	 * 
	 * @param entry the entry to add
	 * @param level the level into which <tt>entry</tt> is inserted
	 */
	protected void insert (Object entry, int level) {
		ORTree.Node node;

		if (path.size()<=level) {
			/** Added by LX */
			if (level == 0) {
				tree.leafNodeNum += 1;
				tree.totalLeafNodeCapacity += tree.singleLeafNodeCapacity;
			}

			path.add(tree.createNode(level));
		}
		node = (ORTree.Node)path.get(level);
		node.entries.add(entry);
		/** Added by LX */
		if (node.level == 0 && !tree.leafNodeToRange.containsKey(node.nodeId)) {
			tree.leafNodeToRange.put(node.nodeId, tree.getDescriptor.invoke(entry));
		}
		// Modified by LX:
		// To accomodate for the initial node capacity
		if (overflows.invoke(node) || (node.level == 0 && node.number()>=initCap*tree.singleLeafNodeCapacity) ||
				(node.level > 0 && node.number() >= initCap * tree.singleIndexNodeCapacity)) {

			ORTree.Node newNode = (ORTree.Node)tree.createNode(node.level);
			Iterator entries;

			newNode.entries.add(entry);
			path.set(level, newNode);
			for (entries = node.entries(); entries.next()!=entry;);
			entries.remove();
			/** Added by LX */
			if (level == 0) {
				tree.leafNodeNum += 1;
				tree.totalLeafNodeCapacity += tree.singleLeafNodeCapacity;
				newNode.nodeId = tree.leafNodeNum;
				tree.leafNodeToRange.put(newNode.nodeId, tree.getDescriptor.invoke(entry));
			}
			insert(save(node), level+1);
		}
	}
	
	/** Saves a node of the tree to external memory.
	 * 
	 * @param node the node to save
	 * @return an <tt>IndexEntry</tt> pointing to the saved node  
	 */
	protected ORTree.IndexEntry save (ORTree.Node node) {
		Container container = (Container)determineContainer.invoke(node);
		Object id = container.insert(node);
		Descriptor descriptor = tree.computeDescriptor(node.entries);
		// LX
		if (this.isRRstarTree) {
			// in case that the tree is too small to call nodeConverter.write()
			((RevisedRTree.Node) node).setMiddlePoint(descriptor);
		}
		return (ORTree.IndexEntry)((ORTree.IndexEntry)tree.createIndexEntry(node.level+1)).initialize(descriptor).initialize(container, id);
	}
}
