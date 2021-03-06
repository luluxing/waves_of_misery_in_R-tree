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

import java.io.PrintWriter;
import java.util.*;

import xxl.core.cursors.Cursors;
import xxl.core.cursors.mappers.Mapper;
import xxl.core.cursors.sources.Enumerator;
import xxl.core.functions.AbstractFunction;
import xxl.core.spatial.points.DoublePoint;
import xxl.core.spatial.rectangles.DoublePointRectangle;
import xxl.core.spatial.rectangles.Rectangle;

/** An <tt>RTree</tt> implementing the linear split-strategy proposed
 * in the original R-Tree paper.  
 * 
 * For a detailed discussion see Guttman, A.: "R-trees: a dynamic index structure
 * for spatial searching", Proc. ACM-SIGMOD Int. Conf. on Management of Data, 47-57, 1984.
 * 
 * @see Tree
 * @see ORTree
 * @see RTree
 */
public class LinearRTree extends RTree {
	private int count = 1;
	private int interval = 1;
	public String outDir;
	private Random rand = new Random(21);

	/* (non-Javadoc)
	 * @see xxl.core.indexStructures.Tree#createNode(int)
	 */
	public Tree.Node createNode (int level) {
		return new Node().initialize(level, new LinkedList());
	}

	/** A modification of {@link RTree.Node node} implementing the linear split-algorithm.
	 *
	 * @see RTree.Node 
	 */
	public class Node extends RTree.Node {
//		public int isBiggerChild = 2;
//		public double myArea;

		/** 
		 *	Performs an RTree split with complexity O(d*n).
		 *	First, the dimension is chosen by maximizing a normalized
		 *	separation distance that is computed for the rectangles
		 *	that have the greatest distance in a dimension.
		 *	Then, for the split dimension, the two most separated rectangles 
		 *	are the seeds for the distribution algorithm.
		 *	Then, each rectangle (one after the other) is added to the partition so
		 *	that the area enlargement of the MBRs of the partitions is minimized.
		 *	Furthermore, the algorithm has to make a distribution so that 
		 *	the number of objects in the new nodes is between the minimal and maximal
		 *	number.
		 *
		 * @param path the nodes already visited during this insert
		 * @return a <tt>SplitInfo</tt> containig all information needed about the split
		 */
		protected Tree.Node.SplitInfo split (final Stack path) {
			final Node node = (Node)node(path);
			final int dimensions = ((Rectangle)rootDescriptor()).dimensions();
			int number = node.number(), minNumberTmp = node.splitMinNumber(), maxNumberTmp = node.splitMaxNumber();
//			System.out.printf("node number=%d, splitmin=%d, splitmax=%d\n", number, minNumberTmp, maxNumberTmp);

			/**
			 * LX:
			 * Randomly choose minNumber and maxNumber in this range
			 */
			int minNumber, maxNumber;
			if (unequalRandomSplit) {
//				System.out.println("Here???");
				int idx = rand.nextInt(Math.abs(maxNumberTmp - minNumberTmp));
				if (minNumberTmp > maxNumberTmp) {
					maxNumber = maxNumberTmp + idx;
					minNumber = number - maxNumber;
				}
				else {
					minNumber = minNumberTmp + idx;
					maxNumber = number - minNumber;
				}
			}
			else {
				minNumber = minNumberTmp;
				maxNumber = maxNumberTmp;
			}

			Iterator seeds = new Mapper(
				new AbstractFunction() {
					Rectangle MBR = rectangle(indexEntry(path));

					public Object invoke (Object object) {
						final int dimension = ((Integer)object).intValue();
						Object e1 = Cursors.minima(node.entries(),
							new AbstractFunction() {
								public Object invoke (Object object) {
									return new Double(rectangle(object).getCorner(true).getValue(dimension));
								}
							}
						).getFirst();
						Object e2 = Cursors.maxima(node.entries(),
							new AbstractFunction() {
								public Object invoke (Object object) {
									return new Double(rectangle(object).getCorner(false).getValue(dimension));
								}
							}
						).getLast();
						Double normalizedSeparation = new Double(
							(rectangle(e2).getCorner(false).getValue(dimension)-rectangle(e1).getCorner(true).getValue(dimension))/
							(MBR.getCorner(true).getValue(dimension)-MBR.getCorner(false).getValue(dimension))
						);

						return new Object [] {normalizedSeparation, e1, e2};
					}
				}
			,new Enumerator(dimensions));
			Object [] seed = (Object[])Cursors.maxima(seeds,
				new AbstractFunction() {
					public Object invoke (Object seed) {
						return ((Object[])seed)[0];
					}
				}
			).getFirst();
			Rectangle oldNodesMBR = new DoublePointRectangle(rectangle(seed[1])), newNodesMBR = new DoublePointRectangle(rectangle(seed[2]));
			
			this.entries.add(seed[2]); // this is definitely in this partition!			
			// the following line has to be here, because if entries.remove is called
			// the following line has to be here, because if entries.remove is called
			// during iteration, it must be considered the case when seed[2] is the
			// last object in the iteration (not very nice for the if-condition).
			node.entries.remove(seed[2]); // object is already in the node's collection
			
			int distrNumber=number-2; // nodes that still have to be distributed
			double areaEnlargementDifference, areaDifference;
			Rectangle rectangle;
			Object entry;

			ArrayList<Double> ori_x = null, ori_y = null, fst_x = null, fst_y = null, scd_x = null, scd_y = null;

			// LX
			if (count % interval == 0) {
				ori_x = new ArrayList<>();
				ori_y = new ArrayList<>();
				fst_x = new ArrayList<>();
				fst_y = new ArrayList<>();
				scd_x = new ArrayList<>();
				scd_y = new ArrayList<>();
			}

			for (Iterator entries = node.entries(); entries.hasNext();) {
				entry = entries.next();
				// LX
				DoublePoint p = null;
				if (count % interval == 0 && node.level() == 0) {
					p = (DoublePoint) entry;
					ori_x.add(p.getValue(0));
					ori_y.add(p.getValue(1));
				}

				if (entry!=seed[1]) {
					// now, entry is not one of the seed objects
					rectangle = rectangle(entry);

					if (
						(node.number()>minNumber) && // still enough objects in the node?
						(
							node.number()==maxNumber || 
							distrNumber+number()<=minNumber ||
							(
								number-number()!=minNumber && (
									(areaEnlargementDifference =
										(Descriptors.union(newNodesMBR, rectangle).area()-newNodesMBR.area())-
										(Descriptors.union(oldNodesMBR, rectangle).area()-oldNodesMBR.area())
									)<0 ||
									areaEnlargementDifference==0 && (
										(areaDifference = newNodesMBR.area()-oldNodesMBR.area())<0 ||
										areaDifference==0 && number()<node.number()
									)
								)
							)
						)
					) {
						// putting objects into the new partition
						this.entries.add(entry);
						// removing it from the old one
						entries.remove();
						if (count % interval == 0 && node.level() == 0) {
							scd_x.add(p.getValue(0));
							scd_y.add(p.getValue(1));
						}
						// compute the new descriptor
						newNodesMBR.union(rectangle);
					}
					else
						// compute the new descriptor
						oldNodesMBR.union(rectangle);
					// element was distributed
					distrNumber--;
				}
			}
			((IndexEntry)indexEntry(path)).descriptor = oldNodesMBR;
			Tree.Node.SplitInfo spi = new SplitInfo(path).initialize(newNodesMBR);
			spi.splitOverlap = newNodesMBR.overlap(oldNodesMBR);// Added by LX
			spi.splitOverlapRatio = spi.splitOverlap / (Math.min(newNodesMBR.area(), oldNodesMBR.area()));
			spi.splitOverlapRatio2 = spi.splitOverlap / (Math.max(newNodesMBR.area(), oldNodesMBR.area()));
			spi.areaDiff = Math.abs(oldNodesMBR.area() - newNodesMBR.area()); // Added by LX
			spi.splitRatio = Math.min(newNodesMBR.area(), oldNodesMBR.area()) / Math.max(newNodesMBR.area(), oldNodesMBR.area());

			if (count % interval == 0 && node.level() == 0) {
				Iterator it = node.entries();
				while (it.hasNext()) {
					DoublePoint d = (DoublePoint) it.next();
					fst_x.add(d.getValue(0));
					fst_y.add(d.getValue(1));
				}
//				writePointToFile(ori_x, ori_y, outDir + "/original_leaf_" + count + ".txt");
//				writePointToFile(fst_x, fst_y, outDir + "/fst_leaf_" + count + ".txt");
//				writePointToFile(scd_x, scd_y, outDir + "/scd_leaf_" + count + ".txt");
			}

			count++;
			return spi;
		}

		private void writePointToFile(ArrayList<Double> x, ArrayList<Double> y, String fileName) {
			try {
				PrintWriter pw = new PrintWriter(fileName);
				for (int i = 0; i < x.size(); i++) {
					pw.printf("%f,%f\n", x.get(i), y.get(i));
				}
				pw.close();
			}
			catch (Exception e) {
				System.out.println("File not found");
			}

		}
	}
}
