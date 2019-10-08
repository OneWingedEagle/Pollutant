package main;
import math.*;
import static java.lang.Math.*;
import io.Console;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.Scanner;

import javax.swing.JOptionPane;

import java.io.FileDescriptor;
import java.io.FileOutputStream;

import components.GUI;


	public class Main implements ActionListener{

	public GUI gui;
	PollutionTransport pt;

	private String path = System.getProperty("user.dir");
	public  boolean console=true;

	public Main()
	{		

		this.gui=new GUI(this.path);
		this.gui.Run.addActionListener(this);


		if(console){
			Console.redirectOutput(this.gui.dataArea);
		}
		
		this.gui.setVisible(true);

	}	
	
	public static void main(String[] args){
		
		
		new Main();
	}

		
	public void runGUI(){
		

		pt=new PollutionTransport();
		
		pt.x0=Double.parseDouble(gui.tfX0.getText());
		pt.xn=Double.parseDouble(gui.tfXn.getText());
	
		
		pt.nT1=Integer.parseInt(gui.tftDiv.getText());
		pt.xDiv=Integer.parseInt(gui.tfxDiv.getText());
		
		
		pt.t0=86400*Double.parseDouble(gui.tfT0.getText());
		pt.tn=86400*Double.parseDouble(gui.tfTn.getText());
		
		pt.D=Double.parseDouble(gui.tfD.getText());
		pt.V=Double.parseDouble(gui.tfV.getText());
		
		pt.Run();
		
	}



	

	@Override
	public void actionPerformed(ActionEvent e)
	{	

		 if(e.getSource()==this.gui.Run){
			this.gui.Run.setBackground(Color.gray);
			this.gui.Run.setEnabled(false);

			runGUI();
			this.gui.Run.setBackground(Color.green);
			this.gui.Run.setEnabled(true);


		}
	

	
	}








}

