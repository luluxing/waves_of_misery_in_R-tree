/* XXL: The eXtensible and fleXible Library for data processing

Copyright (C) 2000-2014 Prof. Dr. Bernhard Seeger
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
package xxl.core.indexStructures.rtrees;

import xxl.core.collections.containers.Container;
import xxl.core.collections.containers.CounterContainer;
import xxl.core.collections.containers.io.BlockFileContainer;
import xxl.core.collections.containers.io.BufferedContainer;
import xxl.core.collections.containers.io.ConverterContainer;
import xxl.core.collections.queues.Queue;
import xxl.core.collections.queues.io.BlockBasedQueue;
import xxl.core.collections.queues.io.QueueBuffer;
import xxl.core.cursors.Cursor;
import xxl.core.cursors.Cursors;
import xxl.core.cursors.mappers.Mapper;
import xxl.core.cursors.sorters.MergeSorter;
import xxl.core.cursors.sources.io.FileInputCursor;
import xxl.core.functions.AbstractFunction;
import xxl.core.functions.Function;
import xxl.core.functions.Functional.NullaryFunction;
import xxl.core.functions.Functional.UnaryFunction;
import xxl.core.functions.Identity;
import xxl.core.indexStructures.ORTree;
import xxl.core.indexStructures.RTree;
import xxl.core.indexStructures.SortBasedBulkLoading;
import xxl.core.indexStructures.Tree;
import xxl.core.indexStructures.rtrees.AbstractIterativeRtreeBulkloader.ProcessingType;
import xxl.core.indexStructures.rtrees.GenericPartitioner.CostFunctionArrayProcessor;
import xxl.core.indexStructures.rtrees.GenericPartitioner.DefaultArrayProcessor;
import xxl.core.io.LRUBuffer;
import xxl.core.io.converters.ConvertableConverter;
import xxl.core.io.converters.Converter;
import xxl.core.predicates.AbstractPredicate;
import xxl.core.spatial.SpaceFillingCurves;
import xxl.core.spatial.SpatialUtils;
import xxl.core.spatial.TestPlot;
import xxl.core.spatial.points.DoublePoint;
import xxl.core.spatial.rectangles.DoublePointRectangle;
import xxl.core.spatial.rectangles.Rectangles;
import xxl.core.util.Pair;

import java.io.*;
import java.util.*;


/**
 * This file is added by LX
 * To demo how to read and write DoublePoint as well as DoublePointRectangle
 *
 *
 */
public class WritePointOrQueryFile {

	private static int seed;

	private static void generateInputFile(boolean isPoint, int rowNumber, int dim, String fileName, boolean normalDist)
			throws IOException {
		System.out.println("Generate points for bulk-loading");
		double [] leftCorner ;
		double [] rightCorner ;

		Random random = new Random(seed);

		DataOutputStream dataOut = new DataOutputStream(new FileOutputStream(fileName));
		for (int j = 0; j < rowNumber; j++) {
			leftCorner = new double[dim];
			rightCorner = new double[dim];
			for (int i = 0; i < dim; i++) {
				if (normalDist) {
					leftCorner[i] = random.nextGaussian() * 0.5 / 3 + 0.5;
					rightCorner[i] = random.nextGaussian() * 0.5 / 3 + 0.5;
				} else {
					leftCorner[i] = random.nextDouble();
					rightCorner[i] = random.nextDouble();
				}
			}
			if (!isPoint) {
				DoublePointRectangle dpr = new DoublePointRectangle(leftCorner, rightCorner);
				dpr.write(dataOut);
			} else {
				DoublePoint dp = new DoublePoint(leftCorner);
				dp.write(dataOut);
			}
		}
	}

	private static void generateFixedSizeQueries(double area, int queryNum, String fileName)
			throws IOException {
		System.out.println("Generate uniform queries");
		double [] leftCorner;
		double [] rightCorner;

		Random random = new Random(seed);
		double width = Math.sqrt(area);
		DataOutputStream dataOut = new DataOutputStream(new FileOutputStream(fileName));
		for (int i = 0; i < queryNum; i++) {
			leftCorner = new double[2];
			rightCorner = new double[2];
			double x = random.nextDouble() * (1 - width) + width / 2;
			double y = random.nextDouble() * (1 - width) + width / 2;
			leftCorner[0] = x - width / 2;
			leftCorner[1] = y - width / 2;
			rightCorner[0] = x + width / 2;
			rightCorner[1] = y + width / 2;

			DoublePointRectangle dpr = new DoublePointRectangle(leftCorner, rightCorner);
			dpr.write(dataOut);
		}

	}

	/**
	 * Generate normal distributed queries around a random focal point
	 * This simulates high density queries in the center
	 * @param area
	 * @param queryNum
	 * @param fileName
	 */
	private static void generateNonUniformQueries(double area, int queryNum, String fileName)
			throws IOException {
		System.out.println("Generate non-uniform queries");
		// First randomly generate a focal point (a, b)
		Random random = new Random(seed);
		double[] focalPoint = {random.nextDouble(), random.nextDouble()};

		// Radius is min(a, b, 1-a, 1-b)
		double radius = Math.min(Math.min(focalPoint[0], focalPoint[1]), Math.min(1-focalPoint[0], 1-focalPoint[1]));
		double sd = radius / 3;

		double width = Math.sqrt(area);
		double [] leftCorner;
		double [] rightCorner;

		DataOutputStream dataOut = new DataOutputStream(new FileOutputStream(fileName));
		int realQueryNum = 0;
		while (realQueryNum < queryNum) {
			leftCorner = new double[2];
			rightCorner = new double[2];
			double x = random.nextGaussian() * sd + focalPoint[0];
			double y = random.nextGaussian() * sd + focalPoint[1];
			leftCorner[0] = x - width / 2;
			leftCorner[1] = y - width / 2;
			rightCorner[0] = x + width / 2;
			rightCorner[1] = y + width / 2;

			if (leftCorner[0] < 0 || rightCorner[0] < 0 || leftCorner[1] > 1 || rightCorner[1] > 1)
				continue;
			DoublePointRectangle dpr = new DoublePointRectangle(leftCorner, rightCorner);
			dpr.write(dataOut);
			realQueryNum++;
		}
	}

	private static void generateFromFile(String readFromFile, String outputFileName) throws FileNotFoundException {
		System.out.println("Generate points for bulk-loading from a csv file");
		DataOutputStream dataOut = new DataOutputStream(new FileOutputStream(outputFileName));
		try (BufferedReader br = new BufferedReader(new FileReader(readFromFile))) {
			String line;
			while ((line = br.readLine()) != null) {
				String[] values = line.split(",");
				DoublePoint dp = new DoublePoint(new double[]{Double.parseDouble(values[0]),
																Double.parseDouble(values[1])});
				dp.write(dataOut);
			}
		} catch (FileNotFoundException e) {
//			e.printStackTrace();
		} catch (IOException e) {
//			e.printStackTrace();
		} catch(NumberFormatException e) {
//			e.printStackTrace();
		}
	}

	public static void main(String[] args) throws IOException {

		if (args.length != 1) {
			System.out.println("Usage: generateFile.properties");
			return;
		}
		try (InputStream input = new FileInputStream(args[0])) {
			Properties prop = new Properties();
			prop.load(input);

			seed = Integer.parseInt(prop.getProperty("Rand_seed"));
			String outputFileName = prop.getProperty("Write_file_name");

			boolean isQuery = Boolean.parseBoolean(prop.getProperty("Generate_query"));
			if (isQuery) {
				int queryNum = Integer.parseInt(prop.getProperty("Query_number"));
				double rectArea = Double.parseDouble(prop.getProperty("Query_rect_area"));
				boolean uniformQuery = Boolean.parseBoolean(prop.getProperty("Uniform_query", "true"));
				if (uniformQuery)
					generateFixedSizeQueries(rectArea, queryNum, outputFileName);
				else
					generateNonUniformQueries(rectArea, queryNum, outputFileName);
			} else {
				boolean isPoint = Boolean.parseBoolean(prop.getProperty("Generate_point_file"));
				boolean normalDist = Boolean.parseBoolean(prop.getProperty("Normal_distribution"));
				int inputNumber = Integer.parseInt(prop.getProperty("Row_number"));
				int dim = Integer.parseInt(prop.getProperty("Dimension"));
				String readFromFile = prop.getProperty("Read_from_file", "");
				if (readFromFile.length() > 0) {
					generateFromFile(readFromFile, outputFileName);
				} else {
					generateInputFile(isPoint, inputNumber, dim, outputFileName, normalDist);
				}
			}

		}
		System.out.println("Done");
	}

}
