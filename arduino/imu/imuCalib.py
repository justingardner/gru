#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 21:34:56 2024

@author: justin
"""

import serial
import serial.tools.list_ports

####################
# class definition
####################
class imuCalib:
    ####################    
    # init function
    ####################
    def __init__(self, baudRate=9600, stepsPerTurn=1600):
        
        # set internal variables
        self.baudRate = baudRate
        self.stepsPerTurn = 1600
        
        # open up serial port
        self.portName = self.selectSerialPort()
        self.s = self.openSerialPort(self.portName)
        
        # set status
        self.sendCommand(b'i')
        self.imuStatus = False
        self.sendCommand(b'c')
        self.colorStatus = False
        
        # setRotation
        self.setRotation();
    ####################    
    # setRotation
    ####################
    def setRotation(self, angle=90, velocity=360):
        # caluclate numSteps
        numSteps = abs((angle / 360) * self.stepsPerTurn);
        # calucluate stepTime (in microseconds)
        stepTime = round(1000000 / ((velocity / 360) * self.stepsPerTurn))
        # keep set values
        self.angle = 360*(numSteps/self.stepsPerTurn)
        self.velocity = 360/((stepTime * self.stepsPerTurn)/(1000000));
        
        if (angle > 0):
            self.sendCommandString(f"r{numSteps},{stepTime};")
        else:
            self.sendCommandString(f"l{numSteps},{stepTime};")

    ####################    
    # go
    ####################
    def go(self):
        # send go command
        self.sendCommand(b'G')   
        # return angle and achieved velocity 
        # FIX - this should get the time value back from arduino
        # and use that to calculate
        return (self.angle,self.velocity)             

    ####################    
    # imu
    ####################
    def imu(self):
        if self.imuStatus:
            # send imu command to turn off IMU display
            self.sendCommand(b'i')            
            self.imuStatus = False
        else:
            # send imu command to turn on IMU display
            self.sendCommand(b'I')
            self.imuStatus = True
            
    ####################    
    # color
    ####################
    def color(self):
        if self.colorStatus:
            # send color command to turn off color display
            self.sendCommand(b'c')
            self.colorStatus = False
        else:
            # send command command to turn on color display
            self.sendCommand(b'C')
            self.colorStatus = True
         
    ####################    
    # selectSerialPort
    ####################
    def selectSerialPort(self):    
        # get ports
        ports = serial.tools.list_ports.comports()
        if not ports:
            print("(imuCalib:selectSerialPort) No serial ports available.")
            return None

        print("Available serial ports:")
        for i, port in enumerate(ports):
            print(f"{i+1}: {port}")

        selection = input("Enter the number corresponding to the port you want to open (or 'q' to quit): ")
        if selection.lower() == 'q':
            return None

        try:
            port_index = int(selection) - 1
            selected_port = ports[port_index].device
            return selected_port
        except (ValueError, IndexError):
            print("Invalid selection. Please enter a valid number.")
            return self.selectSerialPort()

    ######################
    # Open serial port
    #######################
    def openSerialPort(self,port):
        try:
            ser = serial.Serial(port,baudrate=9600)
            print(f"Serial port {port} opened successfully.")
            return ser
        except serial.SerialException as e:
            print(f"Failed to open serial port {port}: {e}")
            return None
    
    #######################
    # Send command over serial port
    #######################
    def sendCommand(self, command):
        try:
            # Send command
            self.s.write(command)
            print(f"(imuCalb:sendCommand) '{command}' sent successfully.")
        except serial.SerialException as e:
            print(f"(imuCalb:sendCommand) Failed to send command '{command}': {e}")

    #######################
    # Send command over serial port
    #######################
    def sendCommandString(self, command):
        try:
            # Send command
            self.s.write(command.encode())
            print(f"(imuCalb:sendCommandString) Command '{command}' sent successfully.")
        except serial.SerialException as e:
            print(f"(imuCalb:sendCommandString) Failed to send command '{command}': {e}")
            
#portName=selectSerialPort()
#p=openSerialPort(portName)
#sendCommand(p,b'I')

