#*********************************************************************
#
# File:          Makefile 
#
# Project:       McPhase
#
# Description:   Makefile for the java module
#
# Author(s):     M. Rotter
#
# Original:	 August 2001
#
#**********************************************************************

all : pack 

compile : IniConfigurator.java IniFile.java IniItem.java IniKey.java
	javac ./*.java

pack    : 
	jar cvf mcphas.jar ./*.class 

clean  :
	rm -vf ./*.class
        rm -vf mcphas.ja
