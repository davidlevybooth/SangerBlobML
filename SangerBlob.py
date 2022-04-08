#!/usr/bin/env python

#---------------------------------------------------------------------------
# SangerBlob.py
# Usage: run SangerBlob.py to Compile Trace data from Sanger Sequences 
# Author: David Levy-Booth
#---------------------------------------------------------------------------


#Generic imports
import os, sys, requests, smtplib, time
from pathlib import Path
from datetime import datetime, date, timedelta

#BioPython imports
from Bio import SearchIO
from Bio import SeqIO
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from xml.etree.ElementTree import XML, fromstring

#Biopython traces
from collections import defaultdict
from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np

#Emailing imports
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase 
from email import encoders

# Datatables
import html
import pandas as pd
from pretty_html_table import build_table


#Config
#--------------------------------------------------------------------------
run_dir1 = "/ubc/smounts/sbic/data/SangerResults/"
run_dir2 = "/ubc/smounts/sbic/data/SangerResults/rf10701-10799"
log_dir  = "logs/SangerBlast"
nres     = 1 #Max number of BLAST hits to return
ethresh  = 0.000001 #BLAST e-value threshold
wait     = 1 #sleep time after BLAST request


#Functions
#--------------------------------------------------------------------------
def run_finder(run_dir): 
    #Capture results directories and pass to list
    for run in os.listdir(run_dir):
        resultsdirs = [f for f in os.listdir(run_dir) if os.path.isdir(os.path.join(run_dir, f))]
        return resultsdirs


def results_finder(resdir, run_dir):
    #Capture results files and pass to list
    resultsfiles = [f for f in os.listdir(os.path.join(run_dir, resdir)) if os.path.isfile(os.path.join(run_dir, resdir, f))]
    resultsfiles = [f for f in resultsfiles if f.endswith(".ab1")]
    return resultsfiles

def extractIntensity(record):
    '''
    Extract the intensity of 'G' bases
    From the raw abi data
    '''
    #print(record.annotations["abif_raw"]["S/N%1"])
    intensities = record.annotations["abif_raw"]["S/N%1"]
    Gspot = intensities[0]
    return f"G({Gspot})"


def extractTrace(record):
    '''
    From Biopython Seq record, 
    extract the Sanger trace.
    '''
    #PBAS1 P2RL1 P2AM1 P1RL1 P1AM1
    channels = ["DATA9", "DATA10", "DATA11", "DATA12"]
    trace = defaultdict(list)
    for c in channels:
        trace[c] = record.annotations["abif_raw"][c]
    return trace


def plotTrace(record):
    # Make an example plot with two subplots...
    q = list(record.annotations["abif_raw"]["PCON2"])
    trace = extractTrace(record)
    rate  = round(len(trace["DATA9"])/len(record.seq))
    xend  = 100 if len(record.seq) > 100 else len(record.seq) 

    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1)
    ax1.bar(range(len(q[0:xend])), q[0:xend], color="grey")
    ax1.axvline(x = 60, color = 'black')
    ax1.axvline(x = 90, color = 'black')
    ax1.set_ylim([0, 65])

    ax2 = fig.add_subplot(2,1,2)
    ax2.plot(trace["DATA9"][0:round(xend*rate):rate], color="black")
    ax2.plot(trace["DATA10"][0:round(xend*rate):rate], color="green")
    ax2.plot(trace["DATA11"][0:round(xend*rate):rate], color="red")
    ax2.plot(trace["DATA12"][0:round(xend*rate):rate], color="blue")
    ax2.axvline(x = 60, color = 'grey')
    ax2.axvline(x = 90, color = 'grey')

    # Save the full figure...
    #fig.savefig(os.path.join(run_dir, log_dir, f"{resdir}.{record.name}.png"))
    #plt.clf() 


def extractTraceData(resfile, run_dir):
    '''
    Find files that include the word "CONTROL" in the 
    ID field, and return the sequences as a dict
    '''
    handle = os.path.join(run_dir, resdir, resfile)
    for record in SeqIO.parse(handle, "abi"):
            q      = list(record.annotations["abif_raw"]["PCON2"])
            trace  = extractTrace(record)
            rate   = round(len(trace["DATA9"])/len(record.seq))
            xstart = 60
            xend   = 100 if len(record.seq) > 100 else len(record.seq)
            
            RecordTrace = {
            "Name": record.name,
            "SeqLen": len(record.seq),
            "Q": tuple(q[xstart:xend]),
            "DATA9": trace["DATA9"][round(xstart*rate):round(xend*rate):rate], 
            "DATA10": trace["DATA10"][round(xstart*rate):round(xend*rate):rate], 
            "DATA11": trace["DATA11"][round(xstart*rate):round(xend*rate):rate], 
            "DATA12": trace["DATA12"][round(xstart*rate):round(xend*rate):rate]}

            return RecordTrace 

           

#Main implementation -- compile trace data
#--------------------------------------------------------------
traces = []

#Capture results directories and pass to list
resultsdir1 = run_finder(run_dir1)
resultsdir2 = run_finder(run_dir2)

#New Data
for resdir in resultsdir1:    
    #Extract results files in directory
    resultsfiles = results_finder(resdir, run_dir1)    

    #Extract traces from sequences as dict
    for resfile in resultsfiles:
        RecordTrace = extractTraceData(resfile, run_dir1)
        RecordTrace["Run"] = resdir
        traces.append(RecordTrace)

#Old Data
for resdir in resultsdir2:    
    #Extract results files in directory
    resultsfiles = results_finder(resdir, run_dir2)    

    #Extract traces from sequences as dict
    for resfile in resultsfiles:
        RecordTrace = extractTraceData(resfile, run_dir2)
        RecordTrace["Run"] = resdir
        traces.append(RecordTrace)


df = pd.DataFrame(traces)
df.to_csv(os.path.join(run_dir1, 'traces.csv'),
    index = False)

# df2 = pd.read_csv(os.path.join(run_dir, 'traces.csv')) 
# print(df2)