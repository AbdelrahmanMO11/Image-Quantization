using System;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Windows.Forms;
using System.Drawing.Imaging;
using System.Diagnostics;
using System.Linq;
///Algorithms Project
///Intelligent Scissors
///

namespace ImageQuantization
{
    /// <summary>
    /// Holds the pixel color in 3 byte values: red, green and blue
    /// </summary>
    /// 
   
    public struct RGBPixel
    {
        public byte red, green, blue;
        public RGBPixel(byte red,byte green,byte blue)
        {
            this.red = red;
            this.green = green;
            this.blue = blue;
        }
        
    }
    public struct RGBPixelD
    {
        public double red, green, blue;
    }
    public struct node
    {
        public RGBPixel parent, current;
        public double distance;
        public int parentIndex,clusterIndex,childIndex;
        
        public node(RGBPixel color1, RGBPixel color2, double distance)
        {
            this.parent = color1;
            this.current = color2;
            this.distance = distance;
            this.parentIndex = -1;
            clusterIndex = -1;
            childIndex = -1;
           
        }
    }

    /// <summary>
    /// Library of static functions that deal with images
    /// </summary>
    public class ImageOperations
    {
        public static int[] cluster;
        public static int[,] distinctIndex;
        /// <summary>
        /// Open an image and load it into 2D array of colors (size: Height x Width)
        /// </summary>
        /// <param name="ImagePath">Image file path</param>
        /// <returns>2D array of colors</returns>
        public static RGBPixel[,] OpenImage(string ImagePath)
        {
            Bitmap original_bm = new Bitmap(ImagePath);
            int Height = original_bm.Height;
            int Width = original_bm.Width;

            RGBPixel[,] Buffer = new RGBPixel[Height, Width];

            unsafe
            {
                BitmapData bmd = original_bm.LockBits(new Rectangle(0, 0, Width, Height), ImageLockMode.ReadWrite, original_bm.PixelFormat);
                int x, y;
                int nWidth = 0;
                bool Format32 = false;
                bool Format24 = false;
                bool Format8 = false;

                if (original_bm.PixelFormat == PixelFormat.Format24bppRgb)
                {
                    Format24 = true;
                    nWidth = Width * 3;
                }
                else if (original_bm.PixelFormat == PixelFormat.Format32bppArgb || original_bm.PixelFormat == PixelFormat.Format32bppRgb || original_bm.PixelFormat == PixelFormat.Format32bppPArgb)
                {
                    Format32 = true;
                    nWidth = Width * 4;
                }
                else if (original_bm.PixelFormat == PixelFormat.Format8bppIndexed)
                {
                    Format8 = true;
                    nWidth = Width;
                }
                int nOffset = bmd.Stride - nWidth;
                byte* p = (byte*)bmd.Scan0;
                for (y = 0; y < Height; y++)
                {
                    for (x = 0; x < Width; x++)
                    {
                        if (Format8)
                        {
                            Buffer[y, x].red = Buffer[y, x].green = Buffer[y, x].blue = p[0];
                            p++;
                        }
                        else
                        {
                            Buffer[y, x].red = p[2];
                            Buffer[y, x].green = p[1];
                            Buffer[y, x].blue = p[0];
                            if (Format24) p += 3;
                            else if (Format32) p += 4;
                        }
                    }
                    p += nOffset;
                }
                original_bm.UnlockBits(bmd);
            }

            return Buffer;
        }
        
        /// <summary>
        /// Get the height of the image 
        /// </summary>
        /// <param name="ImageMatrix">2D array that contains the image</param>
        /// <returns>Image Height</returns>
        public static int GetHeight(RGBPixel[,] ImageMatrix)
        {
            return ImageMatrix.GetLength(0);
        }

        /// <summary>
        /// Get the width of the image 
        /// </summary>
        /// <param name="ImageMatrix">2D array that contains the image</param>
        /// <returns>Image Width</returns>
        public static int GetWidth(RGBPixel[,] ImageMatrix)
        {
            return ImageMatrix.GetLength(1);
        }

        /// <summary>
        /// Display the given image on the given PictureBox object
        /// </summary>
        /// <param name="ImageMatrix">2D array that contains the image</param>
        /// <param name="PicBox">PictureBox object to display the image on it</param>
        public static void DisplayImage(RGBPixel[,] ImageMatrix, PictureBox PicBox)
        {
            // Create Image:
            //==============
            int Height = ImageMatrix.GetLength(0);
            int Width = ImageMatrix.GetLength(1);

            Bitmap ImageBMP = new Bitmap(Width, Height, PixelFormat.Format24bppRgb);

            unsafe
            {
                BitmapData bmd = ImageBMP.LockBits(new Rectangle(0, 0, Width, Height), ImageLockMode.ReadWrite, ImageBMP.PixelFormat);
                int nWidth = 0;
                nWidth = Width * 3;
                int nOffset = bmd.Stride - nWidth;
                byte* p = (byte*)bmd.Scan0;
                for (int i = 0; i < Height; i++)
                {
                    for (int j = 0; j < Width; j++)
                    {
                        p[2] = ImageMatrix[i, j].red;
                        p[1] = ImageMatrix[i, j].green;
                        p[0] = ImageMatrix[i, j].blue;
                        p += 3;
                    }

                    p += nOffset;
                }
                ImageBMP.UnlockBits(bmd);
            }
            PicBox.Image = ImageBMP;
        }


       /// <summary>
       /// Apply Gaussian smoothing filter to enhance the edge detection 
       /// </summary>
       /// <param name="ImageMatrix">Colored image matrix</param>
       /// <param name="filterSize">Gaussian mask size</param>
       /// <param name="sigma">Gaussian sigma</param>
       /// <returns>smoothed color image</returns>
        public static RGBPixel[,] GaussianFilter1D(RGBPixel[,] ImageMatrix, int filterSize, double sigma)
        {
            int Height = GetHeight(ImageMatrix);
            int Width = GetWidth(ImageMatrix);

            RGBPixelD[,] VerFiltered = new RGBPixelD[Height, Width];
            RGBPixel[,] Filtered = new RGBPixel[Height, Width];

           
            // Create Filter in Spatial Domain:
            //=================================
            //make the filter ODD size
            if (filterSize % 2 == 0) filterSize++;

            double[] Filter = new double[filterSize];

            //Compute Filter in Spatial Domain :
            //==================================
            double Sum1 = 0;
            int HalfSize = filterSize / 2;
            for (int y = -HalfSize; y <= HalfSize; y++)
            {
                //Filter[y+HalfSize] = (1.0 / (Math.Sqrt(2 * 22.0/7.0) * Segma)) * Math.Exp(-(double)(y*y) / (double)(2 * Segma * Segma)) ;
                Filter[y + HalfSize] = Math.Exp(-(double)(y * y) / (double)(2 * sigma * sigma));
                Sum1 += Filter[y + HalfSize];
            }
            for (int y = -HalfSize; y <= HalfSize; y++)
            {
                Filter[y + HalfSize] /= Sum1;
            }

            //Filter Original Image Vertically:
            //=================================
            int ii, jj;
            RGBPixelD Sum;
            RGBPixel Item1;
            RGBPixelD Item2;

            for (int j = 0; j < Width; j++)
                for (int i = 0; i < Height; i++)
                {
                    Sum.red = 0;
                    Sum.green = 0;
                    Sum.blue = 0;
                    for (int y = -HalfSize; y <= HalfSize; y++)
                    {
                        ii = i + y;
                        if (ii >= 0 && ii < Height)
                        {
                            Item1 = ImageMatrix[ii, j];
                            Sum.red += Filter[y + HalfSize] * Item1.red;
                            Sum.green += Filter[y + HalfSize] * Item1.green;
                            Sum.blue += Filter[y + HalfSize] * Item1.blue;
                        }
                    }
                    VerFiltered[i, j] = Sum;
                }

            //Filter Resulting Image Horizontally:
            //===================================
            for (int i = 0; i < Height; i++)
                for (int j = 0; j < Width; j++)
                {
                    Sum.red = 0;
                    Sum.green = 0;
                    Sum.blue = 0;
                    for (int x = -HalfSize; x <= HalfSize; x++)
                    {
                        jj = j + x;
                        if (jj >= 0 && jj < Width)
                        {
                            Item2 = VerFiltered[i, jj];
                            Sum.red += Filter[x + HalfSize] * Item2.red;
                            Sum.green += Filter[x + HalfSize] * Item2.green;
                            Sum.blue += Filter[x + HalfSize] * Item2.blue;
                        }
                    }
                    Filtered[i, j].red = (byte)Sum.red;
                    Filtered[i, j].green = (byte)Sum.green;
                    Filtered[i, j].blue = (byte)Sum.blue;
                }

            return Filtered;
        }


        /// <summary>
        /// Get the distinct Colors in the Image
        /// </summary>
        /// <param name="image">Colored image matrix</param>
        /// <returns>list of the distinct colors</returns>

        public static RGBPixel[] colorsIndices;
        public static List<RGBPixel> distinctColor(RGBPixel[,] image) //O(N^2)
        {
            /////////////////// O(1)
            int height = GetHeight(image); //O(1)
            int width = GetWidth(image); //O(1)
            distinctIndex = new int[height, width]; //O(1)
            List<RGBPixel> finalColors = new List<RGBPixel>(); //O(1)
            bool[] colors = new bool[256 * 256 * 256]; //O(1)
            colorsIndices = new RGBPixel[256 * 256 * 256]; //O(1)

            for (int i = 0; i < height; i++)        //o(n^2)
            {
                for (int j = 0; j < width; j++)     //O(N)
                {   ///////// all the next are O(1) //////////////////////

                    int red = image[i, j].red; //O(1)

                    int green = image[i, j].green;//O(1)

                    int blue = image[i, j].blue;//O(1)

                    if (colors[(red * 256 * 256) + (green * 256) + blue] == false) //O(1)
                    {
                        colors[(red * 256 * 256) + (green * 256) + blue] = true; //O(1)
                        finalColors.Add(image[i, j]); //O(1) 
                        colorsIndices[(red * 256 * 256) + (green * 256) + blue] = image[i, j]; //O(1)
                    }
                }
            }
            ///////////////////         total O(N^2)         /////////////////////////////////

            MainForm.distinctColors = finalColors; //O(1)

            return finalColors; //O(1)
        }

        /// <summary>
        /// Find the eculidean distance between two colors
        /// </summary>
        /// <param name="color1">the first color</param>
        /// <param name="color2">the second color</param>
        /// <returns>The eculidean distance calculated</returns>
        public static double getDistance(RGBPixel color1, RGBPixel color2) //O(1) 
        {
            ////// O(1)
            return Math.Sqrt( ((color1.red - color2.red)* (color1.red - color2.red)) + ((color1.green - color2.green) * (color1.green - color2.green)) + ((color1.blue - color2.blue) * (color1.blue - color2.blue)));
        }

        /// <summary>
        /// MST Make or prepare a graph for an Image by getting the distinct colors of it and calculating the distance between each one 
        /// and record this distance
        /// </summary>
        /// <param name="image">RGB Image matrix</param>
        /// <returns>a full detailed graph of the distances as cost of the edge and the distinct colors as vertexes</returns>
        public static node[] mst(RGBPixel[,] image) // O(V^2)
        {
            //timer
           // Stopwatch timer = new Stopwatch();
            //-timer.Start();
            /////////////////////////////////////////////////////////////////

            List<RGBPixel> dColors = distinctColor(image); //O(1)
            cluster = new int[dColors.Count]; //O(1)
            node[] graph2 = new node[dColors.Count];//O(1)
            double mstCost = 0; //O(1)
            int minIndex = 0,//O(1)
            verticesCount = dColors.Count; //O(1)

            bool[] visited = new bool[verticesCount]; //O(1)
            for (int i = 0; i < verticesCount; i++) //O(D)
            {
                graph2[i].current = dColors[i]; //O(1)
                graph2[i].distance = double.MaxValue; //O(1)
            }
            graph2[0].distance = 0;//O(1)


            for (int counter = 0; counter < verticesCount - 1; counter++) //O(n^2)
            {
                int newMin = minIndex; //O(1)
                visited[newMin] = true;//O(1)
                double mindistance = double.MaxValue; //O(1)

                ////////////////////////////////// Get Min Distance and the index of the node satisfying it//////////////////////////////////
                for (int i = 0; i < verticesCount; i++)//O(n)

                {

                    if (visited[i] == false)//O(1)
                    {
                        double distance = getDistance(graph2[i].current, graph2[newMin].current); //O(1)

                        if (distance < graph2[i].distance)//O(1)
                        {
                            graph2[i].distance = distance;//O(1)
                            graph2[i].parent = graph2[newMin].current;//O(1)
                            graph2[i].parentIndex = newMin;//O(1)

                        }
                        ////////////////////////////////// assign min distance and min index for MST//////////////////////////////////       
                        if (graph2[i].distance < mindistance)//O(1)
                        {
                            mindistance = graph2[i].distance;//O(1)
                            minIndex = i;//O(1)

                        }
                    }
                }



                /////////////////////////////////////////////////////////////
                //

            }
            for (int i = 1; i < verticesCount; i++) //O(n)
            {
                graph2[i].childIndex = i; //O(1)
            }
            ////////////////////////////////// calculate MST Cost//////////////////////////////////
            for (int i = 1; i < verticesCount; i++) //O(n)
            {
                mstCost += graph2[i].distance; //O(1)
            }
            MainForm.MSTCostt = mstCost;
            
            // final result of distinct colors for the form
            MainForm.dColorsCounter = dColors.Count;
            //////////////////////////////////////////////////////////////
            return graph2;
        }



        public static List<int>[] adjlist;
        public static char[] BFSchars;
        public static double reds;
        public static double greens;
        public static double blues;
        public static int av;

        public static RGBPixel[] clustering(node[] MSTgraph, int K , List<RGBPixel> distincitColors)// o(k*d)
        {
           

            node[] graph = MSTgraph;
           
            
            int index = 0;
            double maxDistance = -1;

            int[] clustersCounter  = new int[K];
            adjlist = new List<int>[distincitColors.Count];
            for(int i = 0 ; i <distincitColors.Count; i++)
                adjlist[i] = new List<int>();
           
            for (int i = K - 1; i > 0; i--)
            {
                maxDistance = -1;
                for (int j = 0; j < graph.Length; j++)
                {
                    if (graph[j].distance > maxDistance)
                    {
                        maxDistance = graph[j].distance;
                        index = j;
                    }
                }
                graph[index].distance = -1;
                graph[index].parent = graph[index].current;
                graph[index].parentIndex = -1;
                
                clustersCounter[i]++;

            }

           
            for (int i = 0; i < graph.Length; i++)
            {
                if (graph[i].distance != -1)
                {
                    adjlist[graph[i].parentIndex].Add(graph[i].childIndex);
                    adjlist[graph[i].childIndex].Add(graph[i].parentIndex);
                }
            }
              
            //BFS
            BFSchars = new char[distincitColors.Count];
            RGBPixel[] pallette = new RGBPixel[K];
            int currentClusterIndex = 0;
            for (int i = 0;i <distincitColors.Count; i++)
            {
                BFSchars[i] = 'w';
            }
            for (int i = 0; i < distincitColors.Count; i++)
            {   if (BFSchars[i] == 'w')
                {
                    reds = 0.0;
                    blues = 0.0;
                    greens = 0.0;
                    av = 0;
                    BFS(i, distincitColors,currentClusterIndex);
                  
                    RGBPixel color = new RGBPixel((byte)reds , (byte) greens , (byte) blues);
                    
                    pallette[currentClusterIndex++] = color;

                }
            }
           
            return pallette;
        }

        public static void BFS(int adjIndex,List<RGBPixel> distincitColors , int currentIndex)
        {
            
            Queue<int> queue = new Queue<int>();
            queue.Enqueue(adjIndex);
            int k = 0;

            while (queue.Count > 0)
            {
                int i = queue.Dequeue();
                BFSchars[i] = 'g';
                reds += distincitColors[i].red;
                blues += distincitColors[i].blue;
                greens += distincitColors[i].green;
                av++;
                cluster[i] = currentIndex;
                //if(BFSchars[i] == 'w')
                //graph[i].clusterIndex =k;
                foreach (int j in adjlist[i])
                {
                    if(BFSchars[j] == 'w')
                    {
                        queue.Enqueue(j);
                    }
                }
                BFSchars[i] = 'b';
                k++;
            }
            reds = reds / av;
            blues = blues / av;
            greens = greens / av;
        }

        //O(E^2)
        // First Bonus
        public static int kDetect(node[] graph) // O(E^2)
        {
            int k = 0; //O(1)
            double average = 0; //O(1)

            node[] graphTmp = graph; //O(1)

            double totalSTD = double.MaxValue, newTotalSTD = 0; //O(1)
            List<double> distances = new List<double>(); //O(1)

            //Adding all Distances to list
            for (int i = 0; i < graphTmp.Length; i++) //O(n)
            {

                distances.Add(graphTmp[i].distance); //O(1)
            }

            //Infinte loop untill finish with worst case O(E^2)
            for (int count = 0; count > -1; count++)
            {
                double maxSTD = double.MinValue; //O(1)
                int maxIndex = 0; //O(1)

                totalSTD = newTotalSTD; //O(1)
                average = distances.Sum() / distances.Count; //O(1)

                //getting first total Std for the loop
                if (count == 0) //O(1)
                {
                    for (int i = 0; i < distances.Count; i++)  //O(n)
                    {
                        totalSTD += ((distances[i] - average) * (distances[i] - average));  //O(1)
                    }
                    totalSTD = Math.Sqrt(totalSTD / (distances.Count - 1));  //O(1)
                }


                for (int i = 0; i < distances.Count; i++) //O(n)
                {
                    double curSTD = Math.Sqrt(((distances[i] - average) * (distances[i] - average)) / (distances.Count - 1)); //O(1)

                    if (curSTD > maxSTD) //O(1)
                    {
                        maxSTD = curSTD; //O(1)
                        maxIndex = i; //O(1)
                    }

                }


                distances.RemoveAt(maxIndex);  //O(1)
                k++; //O(1)
                average = distances.Sum() / distances.Count; //O(1)
                newTotalSTD = 0; //O(1)
                for (int i = 0; i < distances.Count; i++)  //O(n)
                {
                    newTotalSTD += ((distances[i] - average) * (distances[i] - average)); //O(1)
                }
                newTotalSTD = Math.Sqrt(newTotalSTD / (distances.Count - 1)); //O(1)

                if (Math.Abs(totalSTD - newTotalSTD) <= 0.0001 || distances.Count == 0) //O(1)
                {
                    if (distances.Count != 0)  //O(1)
                    {
                        k++; //O(1)
                    }
                    break;
                }

            }


            MessageBox.Show("number of Clusters = " + k);

            return k;


        }

        public static void colorsIndexing(List<RGBPixel> dColors, RGBPixel[] ColorPallette)//O(D)
        {
            for (int i = 0; i < dColors.Count; i++) // O(D) 
            {

                int red = dColors[i].red; // o 1
                int green = dColors[i].green;// 0 1 
                int blue = dColors[i].blue; // o 1 
                colorsIndices[(red * 256 * 256) + (green * 256) + blue] = ColorPallette[cluster[i]];// o 1 
            }
        }

        public static RGBPixel[,] colorReplacement(RGBPixel[,] image,List<RGBPixel>dColors,RGBPixel[] ColorPallette) //O(N^2)
        {   int width = GetWidth(image); // 1 
            int height = GetHeight(image); // 1 
            colorsIndexing(dColors,ColorPallette);// 1 
            for(int i = 0; i < height; i++)// O(N) ==> O(N^2)
            {
                for (int j = 0; j < width; j++)//O(N)
                {
                    int red = image[i, j].red;// 1 
                    int green = image[i, j].green;// 1 
                    int blue = image[i, j].blue;// 1 

                    image[i, j] = colorsIndices[(red * 256 * 256) + (green * 256) + blue];// 1 

                }
            }

            return image;
        }













        //Second Bonus
        public static RGBPixel[] clusteringBonus(node[] MSTgraph,int k, List<RGBPixel> distincitColors, RGBPixel[,]image)//O(E*N^2)
        {

            double Mean = 0, std = 0;  //O(1)
            int Q = 1; //O(1)

            for (int i = 0; i < MSTgraph.Length; i++)
            {
                Mean += MSTgraph[i].distance; //O(1)
            }
            Mean /= MSTgraph.Length; //O(1)
            for (int i = 0; i < MSTgraph.Length; i++)
            {
                std += ((MSTgraph[i].distance - Mean) * (MSTgraph[i].distance - Mean)); //O(1)
            }
            std = Math.Sqrt(std / (MSTgraph.Length - 1)); //O(1)

            for (int i = 0; i < MSTgraph.Length; i++)
            {
                if (MSTgraph[i].distance > Mean + std)
                {
                    MSTgraph[i].distance = -1; //O(1)
                    MSTgraph[i].parent = MSTgraph[i].current; //O(1)
                    MSTgraph[i].parentIndex = -1; //O(1)
                    Q++; //O(1)
                }
            }
            int[] clustersCounter = new int[k]; //O(1)
            adjlist = new List<int>[distincitColors.Count]; //O(1)
            for (int i = 0; i < distincitColors.Count; i++) //O(n)
                adjlist[i] = new List<int>(); //O(1)

            RGBPixel[] pallette = new RGBPixel[Q + k];
            if (Q <= k)
            {
                if (Q < k)
                {
                    double maxDistance = -1; //O(1)
                    int index = 0; //O(1)
                    for (int i = k - Q; i > 0; i--)
                    {
                        maxDistance = -1; //O(1)
                        for (int j = 0; j < MSTgraph.Length; j++)
                        {
                            if (MSTgraph[j].distance > maxDistance) //O(1)
                            {
                                maxDistance = MSTgraph[j].distance; //O(1)
                                index = j; //O(1)
                            }
                        }
                        MSTgraph[index].distance = -1; //O(1)
                        MSTgraph[index].parent = MSTgraph[index].current; //O(1)
                        MSTgraph[index].parentIndex = -1; //O(1)
                        clustersCounter[i]++; //O(1)
                        Q++; //O(1)
                    }
                }
                for (int i = 0; i < MSTgraph.Length; i++) //O(n)
                {
                    if (MSTgraph[i].distance != -1) //O(1)
                    {
                        adjlist[MSTgraph[i].parentIndex].Add(MSTgraph[i].childIndex); //O(1)
                        adjlist[MSTgraph[i].childIndex].Add(MSTgraph[i].parentIndex); //O(1)
                    }
                }

                //BFS
                BFSchars = new char[distincitColors.Count]; //O(1)

                int currentClusterIndex = 0; //O(1)
                for (int i = 0; i < distincitColors.Count; i++) //O(n)
                {
                    BFSchars[i] = 'w'; //O(1)
                }
                for (int i = 0; i < distincitColors.Count; i++)
                {
                    if (BFSchars[i] == 'w') //O(1)
                    {
                        reds = 0.0; //O(1)
                        blues = 0.0; //O(1)
                        greens = 0.0; //O(1)
                        av = 0; //O(1)
                        BFS(i, distincitColors, currentClusterIndex); //O(V+E)
                        RGBPixel color = new RGBPixel((byte)reds, (byte)greens, (byte)blues);
                        pallette[currentClusterIndex++] = color; //O(1)

                    }
                }
                return pallette;

            }
            else
            {


                for (int i = 0; i < MSTgraph.Length; i++)
                {
                    if (MSTgraph[i].distance != -1) //O(1)
                    {
                        adjlist[MSTgraph[i].parentIndex].Add(MSTgraph[i].childIndex); //O(1)
                        adjlist[MSTgraph[i].childIndex].Add(MSTgraph[i].parentIndex); //O(1)
                    }
                }

                //BFS
                BFSchars = new char[distincitColors.Count]; //O(1)
                int currentClusterIndex = 0; //O(1)
                for (int i = 0; i < distincitColors.Count; i++)
                {
                    BFSchars[i] = 'w'; //O(1)
                }
                for (int i = 0; i < distincitColors.Count; i++)
                {
                    if (BFSchars[i] == 'w') //O(1)
                    {
                        reds = 0.0; //O(1)
                        blues = 0.0; //O(1)
                        greens = 0.0; //O(1)
                        av = 0; //O(1)
                        BFS(i, distincitColors, currentClusterIndex);
                        RGBPixel color = new RGBPixel((byte)reds, (byte)greens, (byte)blues);
                        pallette[currentClusterIndex++] = color; //O(1)

                    }
                }
                RGBPixel[,] recImage = colorReplacement(MainForm.ImageMatrix, distincitColors, pallette);
                List<RGBPixel> dColors = distinctColor(recImage); //O(N^2)
                node[] recGraph = mst(recImage); //O(V^2)
                clusteringBonus(recGraph, k, dColors, recImage);
            }
            return pallette;
        }

    }
}
