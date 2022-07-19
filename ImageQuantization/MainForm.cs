using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Text;
using System.Windows.Forms;
using System.Linq;
using System.Diagnostics;

namespace ImageQuantization
{
    public partial class MainForm : Form
    {
        public MainForm()
        {
            InitializeComponent();
        }

        public static RGBPixel[,] ImageMatrix;

        node[] graph;
        public static List<RGBPixel> distinctColors;
        public static int dColorsCounter;
        RGBPixel[] ColorPallette;
        int KClusters=0;
        public static double MSTCostt;
       Stopwatch timer = new Stopwatch();
     
        private void btnOpen_Click(object sender, EventArgs e)
        {
           
            OpenFileDialog openFileDialog1 = new OpenFileDialog();
            if (openFileDialog1.ShowDialog() == DialogResult.OK)
            {
                //Open the browsed image and display it
                string OpenedFilePath = openFileDialog1.FileName;
                ImageMatrix = ImageOperations.OpenImage(OpenedFilePath);
                ImageOperations.DisplayImage(ImageMatrix, pictureBox1);
            }
            txtWidth.Text = ImageOperations.GetWidth(ImageMatrix).ToString();
            txtHeight.Text = ImageOperations.GetHeight(ImageMatrix).ToString();
        }

        private void btnGaussSmooth_Click(object sender, EventArgs e)
        {
            double sigma = double.Parse(txtGaussSigma.Text);
            int maskSize = (int)nudMaskSize.Value ;
            ImageMatrix = ImageOperations.colorReplacement(ImageMatrix, distinctColors, ColorPallette);
            ImageMatrix = ImageOperations.GaussianFilter1D(ImageMatrix, maskSize, sigma);
            timer.Start();
            ImageOperations.DisplayImage(ImageMatrix, pictureBox2);
            timer.Stop();
            this.textBox1.Text = timer.Elapsed.ToString();


        }

        private void button1_Click(object sender, EventArgs e)
        {
            //timer1.Start();
            timer.Restart();
            graph = ImageOperations.mst(ImageMatrix);
            timer.Stop();
            MessageBox.Show("mst cost = " + MSTCostt.ToString());
            //timer1.Enabled = false;
            this.textBox2.Text = dColorsCounter.ToString();
            List<node> newGraph = graph.ToList();
            newGraph.RemoveAt(0);
            graph = newGraph.ToArray();
            this.textBox1.Text = timer.Elapsed.ToString();

        }

        private void textBox1_TextChanged(object sender, EventArgs e)
        {

        }

        private void button2_Click(object sender, EventArgs e)
        {
            //timer1.Enabled = true;
            timer.Start();
            ColorPallette = ImageOperations.clustering(graph, Int32.Parse(K.Text),distinctColors);
            timer.Stop();
            MessageBox.Show("clustered");
            this.textBox1.Text = timer.Elapsed.ToString();
        }

        private void Kdetect_Click(object sender, EventArgs e)
        {
            //timer1.Enabled = true;
            timer.Start();
            KClusters = ImageOperations.kDetect(graph);
            timer.Stop();
            this.textBox1.Text = timer.Elapsed.ToString();
        }

        private void button4_Click(object sender, EventArgs e)
        {
            timer1.Enabled = true;
            ColorPallette = ImageOperations.clusteringBonus(graph, Int32.Parse(K.Text), distinctColors,ImageMatrix);
            timer1.Stop();
            this.textBox1.Text = timer.Elapsed.ToString();
            MessageBox.Show("clustered");
        }

        private void button5_Click(object sender, EventArgs e)
        {
            SaveFileDialog sf = new SaveFileDialog();
            sf.Filter = "jpg(.jpg)|.jpg";
            if (sf.ShowDialog() == DialogResult.OK)
            {
                pictureBox2.Image.Save(sf.FileName);
            }
        }

        private void pictureBox2_Click(object sender, EventArgs e)
        {

        }

        private void textBox1_TextChanged_1(object sender, EventArgs e)
        {
         
        }

        private void timer1_Tick(object sender, EventArgs e)
        {
         
        }

        private void label7_Click(object sender, EventArgs e)
        {

        }
    }
}