package automation;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.awt.BorderLayout;
import java.awt.Color;
import java.util.Arrays;
import java.util.InputMismatchException;
import java.util.LinkedList;
import java.util.Scanner;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JPanel;

import java.util.Queue;
import java.util.Random; 

@SuppressWarnings("serial")
public class PartitioningAndRouting extends ApplicationFrame  {
		
	static XYSeriesCollection dataset = new XYSeriesCollection( );  /*Plotting*/
	// keys of MST nodes
	static int [] key;
	// parents of MST nodes 
	static int [] parent;
	// vector of partitions 
	static Vector<Partition> v; 
	// cost of total MSTs of all graph 
	static int totalMST;
	// cost of total chain topology edges 
	static int totalCutset;
	// n is the number of elements (vertices) 
	static int n;
	// k is the number of partitions 
	static int k;
	// the  graph matrix 
	static int [][]Graph = new int[n][n];
	// to hold cutset edges
	static LinkedList<ExternalEdge> SecondLevelEdges = new LinkedList<ExternalEdge>();
	// class of objects to store the information of the partitions
	static class Partition
	{		
		// linked list to store vertices that belong to this partition
		LinkedList <Integer> myVertices = new LinkedList<Integer>(); 
		// a n by n matrix that stores the vertices and edges of the vertices in the partition
		int [][] B = new int[n][n];
		// marked partition i.e. connected to another partition or not
		boolean marked = false;
		// constructor
		public Partition(LinkedList<Integer> myVertices, int[][] b, boolean marked)
		{
			super();
			this.myVertices = myVertices;
			this.B = b;
		}
				/* getters and setters */
		public LinkedList<Integer> getMyVertices() {
			return myVertices;
		}
		public void setMyVertices(LinkedList<Integer> myVertices) {
			this.myVertices = myVertices;
		}
		public int[][] getB() {
			return B;
		}
		public void setB(int[][] b) {
			B = b;
		}
		public String toString() 
		{
	        return myVertices.toString() + "\n";
		}
		public boolean isMarked() {
			return marked;
		}
		public void setMarked(boolean marked) {
			this.marked = marked;
		}
	} // end Partition class
	// class of object Tuple < vertex , parent of vertex >
	public static class ExternalEdge 
	{ 		  
		public int v; 
		public  int x; 
		public ExternalEdge(int v, int x) 
		{ 
		   this.v = v; 
		   this.x = x; 
		} 
		public int getV() {
			return v;
		}
		public void setV(int v) {
			this.v = v;
		}
		public int getX() {
			return x;
		}
		public void setX(int x) {
			this.x = x;
		}
		@Override
		public String toString()
		{
			return "(" + v + "," + x + ")";
		}
	} // end ExternalEdge class
	// a method that reverses a given edge's vertices 
	public static ExternalEdge reverse(ExternalEdge e) 
	{
	  int v = e.getX();
	  int x = e.getV();
	  ExternalEdge e2 = new ExternalEdge(v,x);
	  return e2;
	}
	// a method that calculates and returns number of vertices in each partition
	public static int[] calcVertices()
	{
		// calculate number of vertices in each partition 
		int [] num = new int[k];
		// find the result of v / k
		int result = (int)Math.floor(n/k);
		// find the remainder of v / k
		int remainder = n%k;
		// fill out the number of vertices in partitions 
		for (int i=0; i<k; i++) 
		{
			num[i] = result;
		} // end for
		// fill out the remaining vertices 
		while(remainder!=0) 
		{
			for(int i=0; i<k; i++) 
			{
				num[i]++;
				remainder--;
				if(remainder == 0)
					break;
			} // end for
		} // end while 
		// print for testing purposes 
		//System.out.println("\nNumber of vertices in each partition: ");
		// print partitions 
//		for (int i=0; i<k; i++) 
//		{
//			System.out.println("Partition "+ (i+1) +" :"+ num[i] + " vertices.");
//		} // end for
		return num;
	}
	// a method that returns the partition of a given vertex 
	public static int partitionOfVertex(int vertex, Vector <Partition> partitions) 
	{
		LinkedList <Integer> myVertices = new LinkedList<Integer>(); 
		int partition = -1;
		for(int i=0; i<k; i++) 
		{
			myVertices = partitions.get(i).getMyVertices();
			if(myVertices.contains(vertex) == true) 
			{
				partition = i;
				break;
			}
		}
		return partition;
	}
	// a method that finds the initial solution for simulated annealing algorithm
	public static boolean InitialSolution()
	{
		totalCutset = 0;
		totalMST=0;
		SecondLevelEdges.clear();
		// calculate number of vertices in each partition 
		int [] num = calcVertices();
		// create vector of Objects(Partitions) where each node represents a partition 
        v = new Vector <Partition>(k); 
        // queue to store nodes
        Queue <Integer> q = new LinkedList<>(); 
        //add all nodes to queue
		for(int i=0; i<n; i++)
			q.add(i);
		// create a random object
        Random rand = new Random(System.currentTimeMillis());
        // traverse through num[z] number of vertices in each partition
    	int z = 0;
        // while queue is not empty 
        while(!q.isEmpty()) 
		{	
        	// create a new linked list of Integers for the current partition
    		LinkedList <Integer> myVertices = new LinkedList<Integer>(); 
    		// fill this partition with vertices randomly 
        	while(num[z] != 0) 
        	{
        		int randv;
        		do 
        		{
        			// 1- pick a random vertex [ 0 --> n ] that is in queue
        			randv = rand.nextInt(n); // hathi we changed it from 41 to n ( we discovered this bug on 27 nov at 4:30 pm in library)
        		}while(!q.contains(randv));
				// 2- place it in current partition's linked list
        		myVertices.add(randv);
	        	// 3- dequeue it 
				q.remove(randv);
	        	// 4- decrement number of vertices 
				num[z]--;
        	} // end while 
        	// print vertices in this partition
        	//System.out.println("Vertices of partition "+ (z+1) +" : " + myVertices);
        	// place vertices in matrix and store their edges 
        	int [][] A = new int[n][n];
        	for(int i=0;i<n;i++) // for every row
    		{
    			for(int j=0;j<n;j++) // for every column in this row 
    			{
    				if(i<j) // for upper triangle 
    				{
    					if(myVertices.contains(i) && myVertices.contains(j)) 
    						A[i][j] = Graph[i][j];
    					else 
    						A[i][j] = -1; // empty
    				} // end if
    				else
    					A[i][j] = 0; // lower triangle and diagonal not used
    			} // end inner for 
    		} // end outer for 
        	// display matrix of each partition
//        	System.out.println("Matrix of partition "+ (z+1) +" :");
//        	display_matrix(A);
        	boolean connected = BFS_Partition(A ,myVertices.peek(), myVertices);
        	if(connected) 
        	{
	        	// perform Prim algorithm 
	        	totalMST+=MST_Prim_Partition(A,myVertices.peek(), myVertices);
	        	// initialize a new partition
	        	// add myVertices and (null matrix for now)
	        	Partition p = new Partition(myVertices, A, false); 
	        	p.setMyVertices(myVertices);
	        	p.setB(A);
	        	// add the partition to the vector
	        	v.add(p);
	        	z++;
        	}
        	else return false;
        } // end while
        System.out.println("\nThe whole partitioned graph of the initial solution: \n" + v);
/*trademark*/
        
        // find the second level edges between the partitions 
        // will form a chain topology saved in the vecor "Cutset"
        // the size of the vecor will be = number of partitions -1
        	int [] Cutset = new int[k-1];
        	// initialize Cutset 
        	for(int i=0; i<(k-1); i++)
        		Cutset[i] = -1; // has no edges initially 
        	// traverse through partitions 
			for(int i=0; i<(k-1); i++) 
			{
				// store partition number i
				LinkedList <Integer> P1 = v.get(i).myVertices;
				// store a copy of partition i
				LinkedList <Integer> tempP1 = (LinkedList<Integer>) P1.clone();
				// store the partition number i+1
				LinkedList <Integer> P2 = v.get((i+1)).myVertices;
				// store a copy of partition i+1
				LinkedList <Integer> tempP2 = (LinkedList<Integer>) P2.clone();
				// a variable that keeps track of the minimum edge between partition i and partition i+1
				int minPartitionCut=40;
				// initialize the edge -1, -1 means initially empty 
		        ExternalEdge this_edge = new ExternalEdge(-1,-1);
		        // traverse through partition i until its empty 
				while(!tempP1.isEmpty()) 
				{
					// retreive and remove the first vertex in partition i
					int P1vertex = tempP1.poll();
			        // traverse through partition i+1 until its empty 
					while(!tempP2.isEmpty()) 
					{
						// retreive and remove the first vertex in partition i+1
						int P2vertex = tempP2.poll();
						// search for edge with min weight in upper triangle of Graph only 
						if(P1vertex < P2vertex) 
						{
							if( (Graph[P1vertex][P2vertex]) > 0 && (Graph[P1vertex][P2vertex]) <= minPartitionCut)  
							{
								this_edge.setV(P1vertex); 
								this_edge.setX(P2vertex); 
								minPartitionCut = Graph[P1vertex][P2vertex];
							}
						}
						else if (P1vertex > P2vertex) 
						{
							if((Graph[P2vertex][P1vertex]>0) && (Graph[P2vertex][P1vertex] < minPartitionCut))
							{
								this_edge.setV(P2vertex); 
								this_edge.setX(P1vertex);
								minPartitionCut = Graph[P2vertex][P1vertex];
							}
						}
						else
							continue;
					} // finish vertices of P2
				} // finish vertices of P1
				SecondLevelEdges.add(this_edge);
				Cutset[i] = minPartitionCut;
			} // finish visiting all partitions 
/*trademark*/
		for(int i=0; i<(k-1); i++) 
		{
			if(Cutset[i]==-1) 
			{
				System.out.println("The partitioned graph is disconnected.");
        		return false;
			} // end checking  if no slot is empty 
		} // end traversing through vector of second level edges 
		
		// print the edges 
		System.out.println("The second level edges of the initial solution are: "+ SecondLevelEdges);
		//System.out.println("The cutset of the initial solution is : ");
		for(int i=0; i<(k-1); i++) 
		{
			//System.out.println(Cutset[i] +" ");
			totalCutset+=Cutset[i];
		}
		System.out.println("The total cutset of the initial solution = "+ totalCutset);
		System.out.println("The total routing of the initial solution = "+ totalMST);
		System.out.println("The total cost of initial solution: "+Function(totalCutset, totalMST));


        return true; 
            	
	} // end InitialSolution 
		
	// a method that reads from the user the number of elements and the number of partitions	
	public static void read_inputs() 
	{
		Scanner input=new Scanner(System.in);
		// receive the number of elements (vertices)
		do {
		    try {
		        System.out.print("Please enter matrix size (must be positive nonzero integer): ");
		        n = input.nextInt(); // fetch user input
		    } catch (InputMismatchException e) { // catch type mismatch 
		        System.out.print(""); // handle the exception 
		    }
		    input.nextLine(); // clear buffer
		} while (n <= 0 || n!=(int)n ); // while input is invalid
		// receive the number of partitions
		do {
		    try {
		        System.out.print("Please enter number of partitions (must be integer at least 2 and at most "+ (n/2) +" ): ");
		        k = input.nextInt(); // fetch user input
		    } catch (InputMismatchException e) { // catch type mismatch 
		        System.out.print(""); // handle the exception 
		    }
		    input.nextLine(); // clear buffer
		} while (k < 2 || k!=(int)k || k>(n/2) ); // while input is invalid
	}// end read_inputs
	
	// a method that takes integer n and generates and returns a double nxn matrix of uniformly distributed random real numbers 
	public static int[][] generate_matrix(int n) 
	{
		int [] probabilities = new int[4];
		// For a matrix of size n x n, the number of elements in the upper triangle is n * (n - 1) / 2
		int upperElements = (int)n*(n-1)/2;
		int result = (int)Math.floor(upperElements/4);
		int remainder = upperElements%4;
		probabilities[0]= result;
		probabilities[1]= result;
		probabilities[2]= result;
		probabilities[3]= result;
		while(remainder!=0) 
		{
			for(int i=0; i<4; i++) 
			{
				probabilities[i]++;
				remainder--;
				if(remainder == 0)
					break;
			}
		}
		Random rand = new Random();
		//System.out.println("Creating a "+n+"x"+n+" matrix ...");
		int[][] squareMatrix = new int[n][n]; // declare nxn 2D array of type double  
		for(int i=0;i<n;i++) // for every row
		{
			for(int j=0;j<n;j++) // for every column in this row 
			{
				if(i<j) // for upper triangle 
				{
					boolean flag = true;
					do { // fill the matrix with edges having weights subject to these probabilities 
						int temp = rand.nextInt(61) ;
						//System.out.println(temp);
						if(temp == 0 || temp > 40) {
							if(probabilities[0] > 0)
							{
								squareMatrix[i][j]=0;
								probabilities[0]--;
								flag=false;
							}
						}
						else if(temp > 30) {
							if(probabilities[1] > 0)
							{
								squareMatrix[i][j]=temp;
								probabilities[1]--;
								flag=false;
							}
						}
						else if(temp > 10) {
							if(probabilities[2] > 0)
							{
								squareMatrix[i][j]=temp;
								probabilities[2]--;
								flag=false;
							}
						}
						else {
							if(probabilities[3] > 0)
							{
								squareMatrix[i][j]=temp;
								probabilities[3]--;
								flag=false;
							}
						}
					}while(flag); // end while 
				}
				else if(i==j)
					squareMatrix[i][j]=0;
//				else 
//					squareMatrix[i][j]=0;
			}// end inner for
		}// end outer for
	return squareMatrix;
	}// end generate_matrix
	
	// a method that takes an double nxn matrix and returns it with main diagonal set to zero
	public static int[][] set_main_diagonal_values(int[][]M) 
	{
		for(int i=0;i<M.length;i++) // for every row
		{
			for(int j=0;j<M.length;j++) // for every column in this row 
			{
				if(i==j) 
				{//the main diagonal = 0
					M[i][j]=0;
				}
			}//end inner for
		}//end outer for	
		return M;
	}// end set_main_diagonal_values
	
	// a method that takes double nxn matrix and prints it 
	public static void display_matrix(int [][]M)
	{		
		for(int i=0;i<M.length;i++) //for every row
		{
			for(int j=0;j<M.length;j++) //for every column in this row 
			{
				System.out.printf("%4d", M[i][j]); //set grid and print (f% -->six decimal places)
			}//end inner for
			System.out.println(); //new line	
		}//end outer for	
	}// end display_matrix
	
	// a method that receives a matrix that stores weights of edges (type double) and returns corresponding adjacency matrix
	public static int [][] make_adjacency_list(int [][] M) 
	{
		int adjacency_list[][] = new int[M.length][M.length]; // adjacency list implemented as an integer 2D matrix
     
		for(int i=0;i<M.length;i++) //for every row
		{
			for(int j=0;j<M.length;j++) //for every column in this row 
			{
				if ( M[i][j] > 0 ) 
				{
					adjacency_list[i][j]=1; // there is an edge
				}
				else 
				{
					adjacency_list[i][j]=0; // there is no edge
				}
			}// end inner for
		}// end outer for	
		return adjacency_list;
	} // end make_adjacency_list

	// method that performs BFS on a given adjacency list of a an undirected graph and a source node 
	// code implementation is based on CpE-465 graph algorithms hand out no.1
	public static boolean BFS (int [][]AL ,int s) 
	{
		// declare color, distance and parent arrays for nodes
		// each array index describes a feature of a node i.e. color[0] refers to the color of node 0 ,
		// and parent[5] points to the parent of node 5 and so on
		String []color = new String [AL.length];
		int []distance = new int [AL.length];
		int []parent = new int [AL.length];
		
		// 1- initialize colors, distances and parents for all nodes except source node
		for(int i=0;i<AL.length;i++) 
		{
			if(i==s)
				continue; // skip source initializing for now
			else 
			{
				color[i]="white"; // not discovered yet
				distance[i]=-1 ; // -1 refers to infinity 
				parent[i]=-1; // -1 refers to null parent 
			}	
		}
		// 2- initialize source node color, distance and parent
		color[s]="gray"; // discovered ( but not all of it's adjacent nodes yet)
		distance[s]=0; // set distance to root = 0
		parent[s]=-1; // -1 refers to null parent 
		
		Queue<Integer> q = new LinkedList<>(); // queue to store nodes
		q.add(s); // enqueue source node  
		
		// 3- perform BFS
		while(!q.isEmpty()) // while queue is not empty 
		{
			int head = q.peek(); // returns head of queue 
			for(int i=0;i<AL.length;i++)  // traverse through adjacency list 
			{
				if(AL[head][i] == 1 || AL[i][head] == 1 ) // there is a relationship 
				{
					if(color[i] == "white") // if not discovered yet
					{
						color[i] = "gray"; // discovered
						distance[i] = distance[head] + 1; // distance from root = parent distance + 1
						parent[i] = head; // set node parent 
						q.add(i); // enqueue node
					}
				}
					
			} // end for
			q.remove(); // dequeue head
			color[head] = "black"; // explored (discovered with all of it's adjacenct nodes) 
		} // end while
		for (int j=0;j<AL.length;j++) 
		{
			if (!color[j].equals("black")) // all nodes need to be black (explored) else graph is disconnected 
				{
				//System.out.println("Graph is disconnected");
				return false;
				// System.exit(0); // exit current program with 0 (indicating no error has occurred)
				} // end if 
		} // end for
		//System.out.println("Graph is connected");
		return true;
	} // end BFS 
	
	//a method that takes array of vertices keys and returns vertex number with maximum key
	public static int extract_max(int[] key, int [] Q) 
	{
		double maxKey = 0; //initially empty
		// vertex number with maximum key (initially empty)
		int u = -1;
		for(int i=0;i<Q.length;i++)
			if(Q[i]!=-1) 
			{
				u = Q[i];
				break;
			}
		for(int i=0;i<key.length;i++) 
		{
			if( Q[i] != -1 && key[i] >= maxKey) 
			{
				maxKey = key[i]; // keep track of maximum key
				u=Q[i]; // store vertex number with maximum key
			}
		}
		return u;
	} // end extract_max
	
	// a method that checks whether a vertex v belongs to a given set(array) or not
	public static boolean belong(int [] Q, int v) 
	{
		if(Q[v] == -1)
			return false;
		return true;
	} // end belong
	
	// a method that performs Prim algorithm to find the minimum spanning tree
	// code implementation is based on CpE-465 graph algorithms hand out no.1 MST-Prim
	public static void MST_Prim(int [][] AL, int[][] UTM, int r)
	{
		// an empty vector for comparison 
		int [] empty = new int [AL.length];
		// initialize empty vector to -1 to indicate empty slots 
		for (int i=0;i<empty.length;i++)
			empty[i]=-1;  
		// vector Q to hold the vertices of the graph
		int [] Q = new int [AL.length] ;
		// initialize vector Q to hold number of each vertex 
		for (int i=0;i<Q.length;i++)
			Q[i]=i;
		// to hold key of each vertex such that the index of the vector indicates the vertex number and this index stores the key of it
		key = new int [Q.length];
		// to hold parent of each vertex such that the index of vector indicates the vertex number and this index stores the parent of it 
		parent = new int [Q.length];
		// initialize  all keys to a very small number 
		for (int i=0;i<key.length;i++)
			key[i]=-999; 
		// key of root = 999 (max)
		key[r] = 999;
		// parent of root is NIL
		parent[r] = -1;
		
		// keep looping until the vector Q of vertices is empty
		while(!Arrays.equals(Q,empty)) 
		{
			// get the vertex with the maximum key from
			int u = extract_max(key,Q);
			// extract vertex with maximum key from Q
			Q[u]=-1;
			// for all nodes 
			for (int v=0;v<AL.length;v++)
			{
				// if node v is adjacent to node u
				if(AL[u][v] == 1 || AL[v][u] == 1 ) 
				{
					// if v belongs to Q, and edge (u,v) exists (greater than zero), and  is greater than key(v)
					if ( belong(Q,v) && ((UTM[u][v] > key[v] && UTM[u][v] > 0) || (UTM[v][u] > key[v] && UTM[v][u] > 0)) )
					{
						// make u the parent of v
						parent[v]=u;
						// update the key of v
						if (UTM[u][v] > 0 && UTM[u][v] > key[v])
							key[v]=UTM[u][v];
						else
							key[v]=UTM[v][u];
					} // end if
				} // end if
			} // end for
		} // end while
		// print the edges of the found MST
//		System.out.println("MST edges: ");
		// to calculate the overall weight of the MST
		int sum =0;
		for (int i=0;i<parent.length;i++) 
		{
//			System.out.println("( " + i +" , " + parent[i] + " ) , key = "+key[i]);
			if(i!=r)
				sum+=key[i];
		}
//		System.out.printf("MST cost = %d", sum);
	} // end MST_Prim
	
	// a method that performs Prim algorithm to find the MAXIMUM spanning tree
		// code implementation is based on CpE-465 graph algorithms hand out no.1 MST-Prim
		public static int MST_Prim_Partition( int[][] M, int r, LinkedList<Integer> myVertices)
		{
			Integer [] MyVertices = new Integer [myVertices.size()];
			MyVertices = myVertices.toArray(new Integer[0]);
			// an empty vector for comparison 
			int [] empty = new int [myVertices.size()];
			// initialize empty vector to -1 to indicate empty slots 
			for (int i=0;i<empty.length;i++)
				empty[i]=-1;  
			// vector Q to hold the vertices of the graph
			int [] Q = new int [myVertices.size()] ;
			// initialize vector Q to hold number of each vertex 
			for (int i=0;i<Q.length;i++)
				Q[i]=myVertices.get(i);
			// to hold key of each vertex such that the index of the linked list indicates the vertex number and this index stores the key of it
			int [] key = new int [myVertices.size()];
			// to hold parent of each vertex such that the index of linked list  indicates the vertex number and this index stores the parent of it 
			int [] parent = new int [myVertices.size()];
			// initialize  all keys to a very small number 
			for (int i=0;i<key.length;i++)
				key[i]=-999; 
			// key of root = 999 (max)
			key[myVertices.indexOf(r)] = 999;
			// parent of root is NIL
			parent[myVertices.indexOf(r)] = -1;
			
			// keep looping until the vector Q of vertices is empty
			while(!Arrays.equals(Q,empty)) 
			{
				// get the vertex with the maximum key from
				int u = extract_max(key,Q);
				// extract vertex with maximum key from Q
				Q[myVertices.indexOf(u)]=-1;
				// for all nodes 
				for (int v=0;v<M.length;v++)
				{
					// if node v is adjacent to node u
					if(M[u][v] > 0 || M[v][u] > 0) 
					{
						// if v belongs to Q, and edge (u,v) exists (greater than zero), and  is greater than key(v)
						if ( belong(Q,myVertices.indexOf(v)) && ((M[u][v] > key[myVertices.indexOf(v)] && M[u][v] > 0) || (M[v][u] > key[myVertices.indexOf(v)] && M[v][u] > 0)))
						{
							// make u the parent of v
							parent[myVertices.indexOf(v)]=u;
							// update the key of v
							if (M[u][v] > 0 && M[u][v] > key[myVertices.indexOf(v)])
								key[myVertices.indexOf(v)]=M[u][v];
							else
								key[myVertices.indexOf(v)]=M[v][u];
						} // end if
					} // end if
				} // end for
			} // end while
			// print the edges of the found MST
//			System.out.println("MST edges: ");
			// to calculate the overall weight of the MST
			int sum =0;
			for (int i=0;i<parent.length;i++) 
			{
//				System.out.println("( " + MyVertices[i] +" , " + parent[i] + " ) , key = "+key[i]);
				if(i!=myVertices.indexOf(r))
					sum+=key[i];
			}
//			System.out.printf("MST cost = %d\n", sum);
			return sum;
		} // end MST_Prim_Partition
		// method that performs BFS on a given partition and a source node 
		// code implementation is based on CpE-465 graph algorithms hand out no.1
		public static boolean BFS_Partition (int [][]AL ,int r, LinkedList<Integer> myVertices) 
		{
			// declare color, distance and parent arrays for nodes
			// each array index describes a feature of a node i.e. color[0] refers to the color of node 0 ,
			// and parent[5] points to the parent of node 5 and so on
			String []color = new String [myVertices.size()];
			int []distance = new int [myVertices.size()];
			int []parent = new int [myVertices.size()];
			int s = r;
			// 1- initialize colors, distances and parents for all nodes except source node
			for(int i=0;i<myVertices.size();i++) 
			{
				if( myVertices.get(i) == s)
					continue; // skip source initializing for now
				else 
				{
					color[i]="white"; // not discovered yet
					distance[i]=-1 ; // -1 refers to infinity 
					parent[i]=-1; // -1 refers to null parent 
				}	
			}
			// 2- initialize source node color, distance and parent
			color[myVertices.indexOf(s)]="gray"; // discovered ( but not all of it's adjacent nodes yet)
			distance[myVertices.indexOf(s)]=0; // set distance to root = 0
			parent[myVertices.indexOf(s)]=-1; // -1 refers to null parent 
			
			Queue<Integer> q = new LinkedList<>(); // queue to store nodes
			q.add(s); // enqueue source node  
			
			// 3- perform BFS
			while(!q.isEmpty()) // while queue is not empty 
			{
				int head = q.peek(); // returns head of queue 
				for(int i=0;i<AL.length;i++)  // traverse through adjacency list 
				{
					if(AL[head][i] > 0 || AL[i][head] > 0 ) // there is a relationship 
					{
						if(color[myVertices.indexOf(i)] == "white") // if not discovered yet
						{
							color[myVertices.indexOf(i)] = "gray"; // discovered
							distance[myVertices.indexOf(i)] = distance[myVertices.indexOf(head)] + 1; // distance from root = parent distance + 1
							parent[myVertices.indexOf(i)] = head; // set node parent 
							q.add(i); // enqueue node
						}
					}
						
				} // end for
				q.remove(); // dequeue head
				color[myVertices.indexOf(head)] = "black"; // explored (discovered with all of it's adjacent nodes) 
			} // end while
			
			for (int j=0 ; j < myVertices.size() ; j++) 
			{ 
				if ( ( (color[j].equals("white")==true) || (color[j].equals("grey")==true) ) &&  myVertices.contains(j)==true) // all nodes need to be black (explored) else graph is disconnected 
					{
						//System.out.println("The partition is disconnected");
						return false;
					} // end if 
			} // end for
			//System.out.println("The partition is connected");
			return true;
		} // end BFS 
		
		//Method that computes and returns probability of acceptance --> "Acceptance Function" 
		public static double ProbabilityFunc (double w, double t )
		{
			double prob;
			prob=Math.exp(-w/t);
			return prob;
		}//End probabilityFunc  
		//Method that computes and returns new temperature --> " Cooling Function" 
		public static double TemperatureFunc(double t)
		{
			double newT;
			newT=t*0.9; //given function 
			return newT;
		}//End TemperatureFunc
		//Method that performs a function on a given decimal number z --> "Objective Function" 
		public static double Function(int chainEdges, int MST)
		{
			double F;
			F = (double)chainEdges + ((double)1/(double)MST); // calculate cost function 
			return F;
		} //end Function
		//Method that flips a bit (from 0 to 1 and vice versa) in the given string x at a given index i
		public static boolean NeighborFunction()
		{
			totalMST=0;
			totalCutset=0;
			SecondLevelEdges.clear();
			v.clear();
			// calculate number of vertices in each partition 
			int [] num = calcVertices();
			// create vector of Objects(Partitions) where each node represents a partition 
	        v = new Vector <Partition>(k); 
	        // queue to store nodes
	        Queue <Integer> q = new LinkedList<>(); 
	        //add all nodes to queue
			for(int i=0; i<n; i++)
				q.add(i);
			// create a random object
	        Random rand = new Random(System.currentTimeMillis());
	        // traverse through num[z] number of vertices in each partition
	    	int z = 0;
	        // while queue is not empty 
	        while(!q.isEmpty()) 
			{	
	        	// create a new linked list of Integers for the current partition
	    		LinkedList <Integer> myVertices = new LinkedList<Integer>(); 
	    		// fill this partition with vertices randomly 
	        	while(num[z] != 0) 
	        	{
	        		int randv;
	        		do 
	        		{
	        			// 1- pick a random vertex [ 0 --> n ] that is in queue
	        			randv = rand.nextInt(n); //we changed it from 41 to n ( we discovered this bug on 27 nov at 4:30 pm in library)
	        		}while(!q.contains(randv));
					// 2- place it in current partition's linked list
	        		myVertices.add(randv);
		        	// 3- dequeue it 
					q.remove(randv);
		        	// 4- decrement number of vertices 
					num[z]--;
	        	} // end while 
	        	// place vertices in matrix and store their edges 
	        	int [][] A = new int[n][n];
	        	for(int i=0;i<n;i++) // for every row
	    		{
	    			for(int j=0;j<n;j++) // for every column in this row 
	    			{
	    				if(i<j) // for upper triangle 
	    				{
	    					if(myVertices.contains(i) && myVertices.contains(j)) 
	    						A[i][j] = Graph[i][j];
	    					else 
	    						A[i][j] = -1; // empty
	    				} // end if
	    				else
	    					A[i][j] = 0; // lower triangle and diagonal not used
	    			} // end inner for 
	    		} // end outer for 
	        	// display matrix of each partition
	        	boolean connected = BFS_Partition(A ,myVertices.peek(), myVertices);
	        	if(connected) 
	        	{
		        	// perform Prim algorithm 
		        	totalMST+=MST_Prim_Partition(A,myVertices.peek(), myVertices);
		        	// initialize a new partition
		        	// add myVertices and (null matrix for now)
		        	Partition p = new Partition(myVertices, A, false); 
		        	p.setMyVertices(myVertices);
		        	p.setB(A);
		        	// add the partition to the vector
		        	v.add(p);
		        	z++;
	        	}
	        	else return false;
	        } // end while
	        // check if partitioned graph is connected 

	        // find the second level edges between the partitions 
	        // will form a chain topology saved in the vecor "Cutset"
	        // the size of the vecor will be = number of partitions -1
	        	int [] Cutset = new int[k-1];
	        	// initialize Cutset 
	        	for(int i=0; i<(k-1); i++)
	        		Cutset[i] = -1; // has no edges initially 
	        	// traverse through partitions 
				for(int i=0; i<(k-1); i++) 
				{
					// store partition number i
					LinkedList <Integer> P1 = v.get(i).myVertices;
					// store a copy of partition i
					LinkedList <Integer> tempP1 = (LinkedList<Integer>) P1.clone();
					// store the partition number i+1
					LinkedList <Integer> P2 = v.get((i+1)).myVertices;
					// store a copy of partition i+1
					LinkedList <Integer> tempP2 = (LinkedList<Integer>) P2.clone();
					// a variable that keeps track of the minimum edge between partition i and partition i+1
					int minPartitionCut=40;
					// initialize the edge -1, -1 means initially empty 
			        ExternalEdge this_edge = new ExternalEdge(-1,-1);
			        // traverse through partition i until its empty 
					while(!tempP1.isEmpty()) 
					{
						// retreive and remove the first vertex in partition i
						int P1vertex = tempP1.poll();
				        // traverse through partition i+1 until its empty 
						while(!tempP2.isEmpty()) 
						{
							// retreive and remove the first vertex in partition i+1
							int P2vertex = tempP2.poll();
							// search for edge with min weight in upper triangle of Graph only 
							if(P1vertex < P2vertex) 
							{
								if( (Graph[P1vertex][P2vertex]) > 0 && (Graph[P1vertex][P2vertex]) <= minPartitionCut)  
								{
									this_edge.setV(P1vertex); 
									this_edge.setX(P2vertex); 
									minPartitionCut = Graph[P1vertex][P2vertex];
								}
							}
							else if (P1vertex > P2vertex) 
							{
								if((Graph[P2vertex][P1vertex]>0) && (Graph[P2vertex][P1vertex] < minPartitionCut))
								{
									this_edge.setV(P2vertex); 
									this_edge.setX(P1vertex);
									minPartitionCut = Graph[P2vertex][P1vertex];
								}
							}
							else
								continue;
						} // finish vertices of P2
					} // finish vertices of P1
					SecondLevelEdges.add(this_edge);
					Cutset[i] = minPartitionCut;
				} // finish visiting all partitions
				for(int i=0; i<Cutset.length; i++)
					if(Cutset[i] ==-1)
						return false;
				for(int i=0; i<Cutset.length; i++)
					totalCutset+=Cutset[i];
				return true;
		}//end NeighborFunction 
		
		//Method that performs simulated annealing algorithm 	
		@SuppressWarnings("unchecked")
		public static void  SimAnnealing ( double t )
		{
			System.out.println("Start Simulated Annealing...");
			int BestMST=totalMST;
			int BestCutset=totalCutset;
			LinkedList<ExternalEdge> BestEE = (LinkedList<ExternalEdge>) SecondLevelEdges.clone();
			double Ebest=Function(totalCutset, totalMST);; //best objective function result stored here
			Vector <Partition> BestSol = (Vector<Partition>) v.clone();  //global maximum solution
			Vector <Partition> currentSol = new Vector<Partition>(k);  // to store the currrent solution x
			Vector <Partition> neighborSol = new Vector<Partition>(k);  // to store  the neighbor solution x'
			currentSol=(Vector<Partition>) v.clone();
			neighborSol=(Vector<Partition>) v.clone();
			double Temperature=t;
			
			boolean autoSort = false;
			boolean allowDuplicateXValues = true;
			XYSeries series = new XYSeries("Partitions = " + k + " Nodes = "+ n , autoSort, allowDuplicateXValues);
			double E1 = Function(totalCutset, totalMST);
			while(Temperature>0.00001) //Stopping criteria 
			{			 
				System.out.printf("Temperature : %3f", Temperature);
				System.out.printf("%7s", "");
				System.out.printf("Cost = %6f" , E1 );
				System.out.println();
			    series.add( -1*Temperature, Ebest);
			    //dataset.addValue( Ebest , "Temperature" , String.valueOf(Temperature) ); // hatha dakhel el SA loop
//				System.out.println("___________________________________________________\n");
//				System.out.printf("T = %f\n", Temperature);
//				System.out.println("___________________________________________________\n");
//			
				int i=0; // set the for loop controller i to zero for each temperature 
				for(i=0;i<100;i++)	
				{
					
				E1= Function(totalCutset, totalMST); // applies objective function on current solution
				Vector <Partition> temp0 = (Vector<Partition>) currentSol.clone(); //temporary store the current Solution 	
				int tempTotalCutset = totalCutset;
				int tempTotalMST = totalMST;
				LinkedList<ExternalEdge> tempEE = (LinkedList<ExternalEdge>) SecondLevelEdges.clone();
				boolean FeasibleNeighborFound = false;
				do {
					FeasibleNeighborFound= NeighborFunction(); //generate a random neighbor
				}while(FeasibleNeighborFound == false);
				double E2= Function(totalCutset, totalMST); // applies objective function on neighbor solution
				
				Vector <Partition> temp1 = (Vector<Partition>) v.clone(); //temporary store the neighbor Solution 

//				System.out.println("Iteration ("+ i +")");
//				System.out.println("Solution= ");
//				System.out.println(currentSol);
//				System.out.println("Total Cost= " + E2);
				double E_delta= E2-E1; 
				//this is to store and safeguard the global maximum and its result since SA is memoryless
				if(Ebest>E2)
				{
					Ebest=E2;
					BestSol= (Vector<Partition>) temp1.clone();
					BestMST=totalMST;
					BestCutset=totalCutset;
					BestEE= (LinkedList<ExternalEdge>) SecondLevelEdges.clone();
				}
				if(Ebest>E1)
				{
					Ebest=E1;
					BestSol= (Vector<Partition>) temp0.clone();
					BestMST=tempTotalMST;
					BestCutset=tempTotalCutset;
					BestEE=(LinkedList<ExternalEdge>) tempEE.clone();
										
					totalMST = tempTotalMST;
					totalCutset = tempTotalCutset;
					SecondLevelEdges = (LinkedList<ExternalEdge>) tempEE.clone();
				}
				// fe problem bel second level edges of the neighbour solution 
				//this is the critera of deciding wether to accept or reject the neighbor solution
				if (E_delta<0 || ProbabilityFunc(E_delta,Temperature)<0.9)
					if(E1>E2)
						currentSol=(Vector<Partition>) temp1.clone();	 //accepted 
				// else keep the current solution 

				}// end for loop
				Temperature=TemperatureFunc(Temperature);	// cooling function call 
			}//end while 
			System.out.println("The Minimum Solution Found So Far Is : ");
			System.out.println(BestSol);
			System.out.println("With External Edge(s):" + BestEE);
			System.out.println("\nWith Cutset = " + BestCutset);
			System.out.println("And Total Routing Cost = " + BestMST);
			System.out.printf("And Total Cost = %6f", Ebest);
			dataset.addSeries(series);
		}//End SimAnnealing 
		
		////	PLOTTING 	/////

		   public PartitioningAndRouting( String applicationTitle , String chartTitle ) {
		        super("Temperature vs Cost Function For "+ n + "Vertices and "+k+" Partitions");

			   JPanel chartPanel = createChartPanel();
		        add(chartPanel, BorderLayout.CENTER);
		        
		        setSize(640, 480);
		        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		        setLocationRelativeTo(null);
		   }
		   private JPanel createChartPanel() {
			    String chartTitle = "Temperature vs Cost Function For "+ n + " Vertices and "+k+" Partitions";
			    String xAxisLabel = "Temperature";
			    String yAxisLabel = "Cost Function";
			 			 
			    JFreeChart chart = ChartFactory.createXYLineChart(chartTitle,
			            xAxisLabel, yAxisLabel, dataset, PlotOrientation.VERTICAL, false, false, false);
			 
			    return new ChartPanel(chart);
			}


	public static void main(String[] args) 
	{
//		for(int i=0; i<3;i++) 
//		{
			boolean startSA = false;
			// read user input
			read_inputs();
			// then generate nxn matrix and display it 
			Graph = generate_matrix(n);
	//		System.out.println("Initial Solution:");
			display_matrix(Graph); 
			// call adjacency_list to create adjacency matrix called AL from M 
			int [][] AL = make_adjacency_list(Graph); 
			// call BFS to perform BFS: arguments --> adjacency matrix AL , source node (2 or any node)
			if (BFS(AL,2)) 		
			{
				// generate the initial solution 
				do {
					startSA = InitialSolution();
				}while(startSA == false);
				SimAnnealing(10000);
		}
//			
//			do {
//				startSA = InitialSolution();
//			}while(startSA == false);
//			SimAnnealing(750);
//			
//			do {
//				startSA = InitialSolution();
//			}while(startSA == false);
//			SimAnnealing(750);
//			
			/*Plotting*/
			
		    PartitioningAndRouting chart = new PartitioningAndRouting(
		    	         "Temperature Vs Cost Function" ,
		    	         "Temperature Vs Cost Function");

		    	      chart.pack( );
		    	      RefineryUtilities.centerFrameOnScreen( chart );
		    	      chart.setVisible( true );
		    	      chart.setBackground(Color.BLACK);
       
		//} // end if
	} // end main
} // end class

