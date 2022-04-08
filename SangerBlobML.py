#!/usr/bin/env python

#---------------------------------------------------------------------------
# SangerBlobML.py
# Usage: run SangerBlobML.py to assess  
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

#tensorflow imports
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
tf.__version__


#Config
#--------------------------------------------------------------------------
run_dir = "/ubc/smounts/sbic/data/SangerResults/"
log_dir  = "logs/SangerBlob"

#Load tensorflow model 
model = keras.models.load_model(os.path.join(run_dir, 'DyeBlobModelV1'))
print(model.summary())

#Channels for dye blob prediction
channels = ["DATA9", "DATA10", "DATA11", "DATA12"]
labels = ["", "dye blob"]

#Functions
#--------------------------------------------------------------------------
def dictToArray(dict1):
    result = dict1.items()
    # Convert object to a list
    data = list(result)
    # Convert list to an array
    return np.array(data)

def subsetDict(keys, dict1):
    dict2 = {key: dict1[key] for key in keys}
    return dict2

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
            "ID": record.id,
            "SeqLen": len(record.seq),
            "Q": tuple(q[xstart:xend]),
            "DATA9": trace["DATA9"][round(xstart*rate):round(xend*rate):rate], 
            "DATA10": trace["DATA10"][round(xstart*rate):round(xend*rate):rate], 
            "DATA11": trace["DATA11"][round(xstart*rate):round(xend*rate):rate], 
            "DATA12": trace["DATA12"][round(xstart*rate):round(xend*rate):rate]}

            return RecordTrace 


#send email
def emailControls(resdir, controls): 
    subject = f"{resdir} SangerBlobML results" 
    sender = "SangerBlobML"
    sender_address = f'{sender}@sbcprodredhat.sequencing.ubc.ca'
    smtp_server = 'localhost'
    #names, emails = get_contacts('/home/sbcuser/mail/contacts.txt') # Read contacts
    #names = ["David"]
    emails = ["dlevyboo@mail.ubc.ca"]
    names = ["SBC Members"]
    #emails = ["sanger.submissions@ubc.ca"]  
    
    #controltable = controls.to_html(escape=False)
    controltable = build_table(controls, 'blue_light', font_size = 'small', escape = False)
     
    # Set up the SMTP server
    server = None
    try:
        server = smtplib.SMTP(smtp_server)    

        # For each contact, send the email:
        for name, email in zip(names, emails):

            message = f'''
            <html><body>
            <p>Hello {name},</p>
            <p>Here are the SangerBlobML results for:</p>
            <h3>run: {resdir}</h3>
            <p></p>
            <p><b>SangerBlobML uses AI to identify dye blobs in Sanger sequencing traces</b></p>
            <div id="controls" style>
            {controltable}
            </div>
            </br>
            {datetime.now().strftime('%Y-%m-%d - %H:%M:%S')}
            </body></html>
            '''
            msg = MIMEMultipart() # Create a message    

            # Setup the parameters of the message
            msg['From'] = sender_address
            msg['To'] = email
            msg['Subject'] = subject    

            # Add in the message body
            msg.attach(MIMEText(message, 'html'))    

            # Send the message via the server set up earlier.
            server.send_message(msg)
            del msg
    except Exception as e:
        print(datetime.now().strftime('%Y-%m-%d - %H:%M:%S'))
        print(e)
    finally:
        # Terminate the SMTP session and close the connection
        if server:
            server.quit()
   


#Main implementation -- apply SangerBlob.py
#--------------------------------------------------------------
if __name__ == '__main__':
    print("Starting SangerBlobML")
    
    #Capture results directories and pass to list
    resultsdirs = run_finder(run_dir)
    
    #Capture results files and pass to list
    for resdir in resultsdirs:    
        '''
        Results Directory | Primary level for emails
        '''
        log_path = os.path.join(run_dir, log_dir, f"{resdir}.csv")
        if not os.path.isfile(log_path): 
            #Extract results files in directory
            resultsfiles = results_finder(resdir, run_dir)    

            #Extract traces from sequences as array
            traces = []

            for resfile in resultsfiles:
                RecordTrace = extractTraceData(resfile, run_dir)
                RecordTrace["Run"] = resdir
                '''
                Isolate data channels for model input. 
                Should be a list of arrays
                ''' 
                record = subsetDict(["Run", "Name", "ID", "SeqLen"], RecordTrace)
                trace = subsetDict(["Q"] + channels, RecordTrace)
                test_features = []
                traceArray = np.array([[[float(i) for i in trace['Q']], 
                     [float(i) for i in trace['DATA9']], 
                     [float(i) for i in trace['DATA10']], 
                     [float(i) for i in trace['DATA11']],
                     [float(i) for i in trace['DATA12']]]])
                print(RecordTrace["Name"])

                try:
                    #Use logistic regression result for probablilty
                    p = model.predict(traceArray)
                    probability = [round(float(x), 2) for x in p][0]
                    prediction  = int([round(float(x)) for x in p][0])
                    prediction  = labels[prediction]

                    #Add to sequence dict
                    record["DyeBlob"] = f"<b>{prediction}</b>"
                    record["DyeBlobProb"] = f"{probability}%"
                except:
                    print("No prediction")
                    record["DyeBlob"] = "no prediction"
                    record["DyeBlobProb"] = ""
                traces.append(record)


            if traces:
                df = pd.DataFrame(traces)
                df = df.sort_values(by="Name")
                print(f"Writing SangerBlobML log {log_path} on {datetime.now().strftime('%Y-%m-%d - %H:%M:%S')}")
                df.to_csv(log_path)
                print(f"Emailing SangerBlobML for run {resdir} on {datetime.now().strftime('%Y-%m-%d - %H:%M:%S')}")
                emailControls(resdir, df)
                print("")



