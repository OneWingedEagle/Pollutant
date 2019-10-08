package components;


import java.awt.Color;
import java.awt.Component;
import java.awt.ComponentOrientation;
import java.awt.Dimension;
import java.awt.FileDialog;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Frame;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

import main.Main;
import math.util;

public class GUI extends JFrame implements ActionListener{


	
	public		JPanel		panel,pp1;
	public JTextArea dataArea=new JTextArea(), iccgArea=new JTextArea();
	public  TextField tfD,tfV,tfX0,tfXn,tfT0,tfTn,tfxDiv,tftDiv;
	public Label lbD,lbV,lbX0,lbXn,lbT0,lbTn,lbxDiv,lbtDiv;
	public  Button Run,bTerminate;
	
	public static int screenWidth,screenHeight;

		
	 
	 public GUI(String path) {

		
				Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
				screenWidth = (int)(screenSize.getWidth());
				screenHeight = (int)(screenSize.getHeight());
		
				int d2=(int)(100.*GUI.screenWidth/1200);
				

			panel = new JPanel(new FlowLayout(0,10,10));
			panel.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
			getContentPane().add(panel);
			setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			setTitle(" FEM Analysis : "+path);
		
			int width = (int)(.6*screenWidth);
			int height = (int)(.9*screenHeight);
			
			setSize(width,height);
			setLocation((int)(.3*screenWidth),(int)(.05*screenHeight));
	
		 
		    
	 //================================================================ redirecting console to text area
			
			int ww1=(int)(.9*width);
			int hh1=(int)(.7*height);
		
			dataArea.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
			dataArea.setEditable(false);;
			dataArea.setFont(new Font("Arial", 0, (int)(10.0*screenWidth/1200)));
		
			dataArea.setBorder(BorderFactory.createLineBorder(Color.blue,1));
			   JScrollPane scrollPane = new JScrollPane(dataArea);
			   scrollPane.setBorder(BorderFactory.createCompoundBorder(
						BorderFactory.createTitledBorder("Progress"),
						BorderFactory.createEmptyBorder(10,5,5,5)));
		  scrollPane.setPreferredSize(new Dimension(ww1,hh1));
		  
		  iccgArea.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
		  iccgArea.setEditable(false);;
			iccgArea.setFont(new Font("Arial", 0, (int)(10.0*screenWidth/1200)));
			
			iccgArea.setBorder(BorderFactory.createLineBorder(Color.blue,1));
			   JScrollPane scrollPane2 = new JScrollPane(iccgArea);
			   scrollPane2.setBorder(BorderFactory.createCompoundBorder(
						BorderFactory.createTitledBorder("Parameters"),
						BorderFactory.createEmptyBorder(10,5,5,5)));
		  scrollPane2.setPreferredSize(new Dimension((int)(.3*width),hh1));

	//=================================================================== configuring the panel		

				Label lbMeshFile=new  Label("Load Mesh"   , Label.RIGHT);
				Label lbDataFile=new  Label("Load Data"   , Label.RIGHT);
				
				int d1=(int)(80.*GUI.screenWidth/1200);
				int h1=(int)(20.*GUI.screenWidth/1200);
				lbMeshFile.setPreferredSize(new Dimension(d1,h1));
				lbDataFile.setPreferredSize(new Dimension(d1,h1));
	
				
				lbD=new Label("  D [m2/s]");
				
				lbV=new Label("   V [m/s]");
				lbX0=new Label("   X0 [m]");
				lbXn=new Label("   Xn [m]");
				
				lbT0=new Label("   T0 [days]");
				lbTn=new Label("   Tn [days]");
				
				lbxDiv=new Label("   X-Div ");
				lbtDiv=new Label("   t_Div ");
				
				tfD=new TextField("1e-3");
				tfD.setPreferredSize(new Dimension(d2/3,h1));
	
				tfV=new TextField("1e-4");
				tfV.setPreferredSize(new Dimension(d2/3,h1));
				
				tfX0=new TextField("0.0");
				tfX0.setPreferredSize(new Dimension(d2/3,h1));

				tfXn=new TextField("1000.0");
				tfXn.setPreferredSize(new Dimension(d2/3,h1));

				
				tfT0=new TextField("0.0");
				tfT0.setPreferredSize(new Dimension(d2/3,h1));

				
				tfTn=new TextField("10.");
				tfTn.setPreferredSize(new Dimension(d2/3,h1));

				
				tfxDiv=new TextField("100");
				tfxDiv.setPreferredSize(new Dimension(d2/3,h1));

				tftDiv=new TextField("10");
				tftDiv.setPreferredSize(new Dimension(d2/3,h1));


				panel.setComponentOrientation(ComponentOrientation.LEFT_TO_RIGHT);
	
				Label lbsaveTo=new Label("Save output to :", Label.RIGHT);
				lbsaveTo.setFont(new Font("Arial", 1, 16));
				
				 Run=new Button("  Run  ");
				Run.setBackground(Color.GREEN);
				Run.setPreferredSize(new Dimension(d2/3,h1));
				 
				Label empty1=new Label();
				empty1.setPreferredSize(new Dimension(30,3));
				Label empty2=new Label();
				empty2.setPreferredSize(new Dimension(30,3));

				JPanel filesPanel1 = new JPanel(new GridLayout(1,7,5,5));
				
				filesPanel1.add(lbX0);
				filesPanel1.add(tfX0);
				filesPanel1.add(lbXn);
				filesPanel1.add(tfXn);
				
				filesPanel1.add(lbT0);
				filesPanel1.add(tfT0);
				
				JPanel filesPanel2 = new JPanel(new GridLayout(1,8,5,5));
				filesPanel2.add(lbTn);
				filesPanel2.add(tfTn);
				
				filesPanel2.add(lbD);
				filesPanel2.add(tfD);
				
				filesPanel2.add(lbV);
				filesPanel2.add(tfV);

				filesPanel2.add(empty1);
				filesPanel2.add(empty2);
				filesPanel2.add(Run);

				
				JPanel iterPanel = new JPanel(new FlowLayout(0,10,10));

				
				JPanel filesPanel = new JPanel(new GridLayout(2,1,10,10));
				filesPanel.add(filesPanel1);
				filesPanel.add(filesPanel2);
				JPanel textPanel = new JPanel(new GridLayout(2,1,10,10));
				textPanel.add(scrollPane);
				//textPanel.add(scrollPane2);
				
				panel.add(filesPanel);
				panel.add(iterPanel);
				panel.add(textPanel);
			
				
				 //======================  file paths	
						
						
						String meshFile1= "";
						String dataFile1= "";
						
				
					this.update(null);
					//====================================
				
			

	}
	 
		public static void main(String[] args){

			new Main();
		}

		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			
		}
	 


}




