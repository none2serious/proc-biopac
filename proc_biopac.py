#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 16:03:05 2017

@author: stan
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import bioread as br
import os
impot argparse
import sys
import matplotlib.gridspec as gridspec
import os.path as path

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def gracefulExit():
    parser.print_help()
    exit(2)

def new_physio_struct2(tasklist):
    physio_data = {}
    physio_data['pulse'] = []
    physio_data['resp'] = []
    subject_data = {}
    subject_data['subject_id'] = ""
    for task in tasklist:
        subject_data[task] = physio_data
    return(subject_data)

def new_physio_struct(tasklist):
    physio_data = {}
    physio_data['subject_id'] = ""
    for task in tasklist:
        physio_data[task] = []
    return(physio_data)

def lfr(X, order = 3, cutoff_freq = 0.01):
    #Identify low frequency drift and remove 
    B, A = signal.butter(order, cutoff_freq, output='ba') 
    # get the low-freq component
    tempf = signal.filtfilt(B,A, X)
    #remove it from the raw data
    Xnew = X - tempf
    return(Xnew, tempf)

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def proc_acq(acq, target_sampling_rate=50):
    trigger, pulse, resp = get_channels(acq)
    idx = np.where(trigger)
    samples_per_second = int(acq.samples_per_second)
    pulse = pulse[idx]
    resp = resp[idx]
    # resample data to match target sampling rate
    if samples_per_second != 50:
        scale = int(samples_per_second / target_sampling_rate)
        pulse = signal.decimate(pulse,scale,zero_phase=True)
        resp = signal.decimate(resp,scale,zero_phase=True)
    # de-spike data.
    d_pulse = signal.medfilt(pulse,3)
    d_resp = signal.medfilt(resp,3)
#    Low-frequency filter for pulse data ; resp belt is so loose on these kids 
#    That there is really only drift on the peak values. 
    d_pulse, _ = lfr(d_pulse)
    # find peaks for resp and pulse signals, compute # of peaks to peaks/min
    #  N.B. May need to adjust window width and preprocessing for some Ss/populations.
    hr_peaks = signal.find_peaks_cwt(d_pulse, np.arange(1,35))
    hr = int(len(hr_peaks) / ((len(d_pulse)/target_sampling_rate)/60.0))
    rr_peaks = signal.find_peaks_cwt(moving_average(d_resp,50), np.arange(1,70))
    rr = int(len(rr_peaks) / ((len(d_resp)/target_sampling_rate)/60.0))
        
    return(d_pulse, d_resp, hr, rr, hr_peaks, rr_peaks)
    
def get_channels(acq):
    trig = np.nan
    pulse = np.nan
    resp = np.nan
    pulse_str = u'PPG100C'
    trigger_str = u'Digital input'
    resp_str = u'Custom, DA100C'
    for k in range(len(acq.channels)):
        cname = acq.channels[k].name
        if cname == "TRIGGER" or cname == trigger_str:
            trig = acq.channels[k].data
        elif cname == "PULSE" or cname == pulse_str:
            pulse = acq.channels[k].data
        elif cname == "RESP" or cname == resp_str:
            resp = acq.channels[k].data                
    return(trig, pulse, resp)

def plot_subject_struct(physio_data, output_dir='~/'):
    fig = plt.figure(num=None, figsize=(12, 12), dpi=72, facecolor='w', edgecolor='k')
    fig.suptitle(physio_data['subject_id'], fontsize=16)
    gs = gridspec.GridSpec(4, 2)
    for k in range(4):
        ax = fig.add_subplot(gs[k,0])
        rate = physio_data[tasklist[k]][0]
        if rate is not np.nan:
            try:
                dat = physio_data[tasklist[k]][1]
                ticdat = np.zeros(len(dat))
                dat = dat[0:500]
                dat = dat - dat.mean()
                idx = physio_data[tasklist[k]][2]
                ticdat[idx] = dat.max()
                ticdat = ticdat[0:500]
                ax.plot(dat)
                ax.plot(ticdat,'r-')
                ax.set_title(tasklist[k]+", Pulse="+str(rate))
            except:
                ax.set_title(tasklist[k]+", Pulse=Unknown")
                ax.text(0.25,0.4, "ERROR OCCURRED DURING PROCESSING")
        else:
            ax.set_title(tasklist[k]+", Pulse=Unknown")
            ax.text(0.25,0.4, "FILE NOT FOUND,\nOR PROCESSED WITH ERROR")
            
        ax = fig.add_subplot(gs[k,1])
        rate = physio_data[tasklist[k]][3]
        if rate is not np.nan:
            dat = physio_data[tasklist[k]][4]
            ticdat = np.zeros(len(dat))
            dat = dat[0:500]
            dat = dat - dat.mean()
            idx = physio_data[tasklist[k]][5]
            ticdat[idx] = dat.max()
            ticdat = ticdat[0:500]
            ax.plot(dat)
            ax.plot(ticdat,'r-')
            ax.set_title(tasklist[k]+", Resp= "+str(rate))
        else:
            ax.set_title(tasklist[k]+", Resp= Unknown")
            ax.text(0.25,0.4, "FILE NOT FOUND,\nOR PROCESSED WITH ERROR")
    gs.update(wspace=0.1, hspace=0.4)
    figname = physio_data['subject_id'] + "_physio.png"
    outfile = os.path.join(output_dir,figname)
    fig.savefig(outfile)

def save_physio_csv(physio_data, output_dir, tasklist):
    output_csv = os.path.join(output_dir, 'physio_values.csv')
    if os.path.isfile(output_csv) is False:
        hdrstr = "subject_id"
        for task in tasklist:
            hdrstr = hdrstr+','+task+'_hr,'+task+'_rr'
        csv_out = open(output_csv, 'w')
        csv_out.write(hdrstr+'\n')
        csv_out.close()
    outstr = physio_data['subject_id']
    for task in tasklist:
        outstr = outstr+','+str(physio_data[task][0]) + ',' + str(physio_data[task][3])

def save_physio_tsv(physio_data, output_dir, tasklist):
    subid = physio_data['subject_id']
    tsv_name = {}
    tsv_name['rest1'] = subid+'_task-rest_run-01_physio'
    tsv_name['rest2'] = subid+'_task-rest_run-02_physio'
    tsv_name['face1'] = subid+'_task-faces_run-01_physio'
    tsv_name['face2'] = subid+'_task-faces_run-02_physio'    
    
    physio_json = '''{
        "SamplingFrequency": 50.0,
        "StartTime": 0.00,
        "Columns": ["cardiac", "respiratory"]
        }'''
    for task in tasklist:
        print(task)
        hr = physio_data[task][0]
        rr = physio_data[task][3]
        if hr is not np.nan and rr is not np.nan:
            outfile = os.path.join(output_dir,subid,"func",tsv_name[task]+'.tsv')
            k = []
            k.append(physio_data[task][1][:])
            k.append(physio_data[task][4][:])
            p = np.array(k)
            np.savetxt(outfile, p.transpose(), delimiter='\t')
            jsonpath = os.path.join(output_dir, tsv_name[task]+'.json')
            jsonfile = open(jsonpath,'w')
            jsonfile.write(physio_json)
            jsonfile.close()
        
parser = MyParser(prog="proc_biopac")
parser.add_argument('-pd','--project_directory', help="The project folder containing sourceData, etc. (Required)")
parser.add_argument('-out','--output_directory', help="The project folder containing sourceData, etc.",default='BIDS')
parser.add_argument("-qc", "--qc_directory", help="Destination folder for processed HR and resp data", default='QA')
parser.add_argument("-sl", "--subject_list", help="File name listing subjects to process", default='subject_list.txt')
parser.add_argument("-tl", "--task_list", help="File name containing tasks Ss performed while BIOPAC was running", default='physio_tasklist.txt')
args = parser.parse_args()

if (args.project_directory is None or args.output_directory is None):
    gracefulExit()
    
sourceDir = args.project_directory
outDir = args.output_directory
sublist = [fil.strip() for fil in open(args.subject_list)]
tasklist = [fil.strip() for fil in open(args.task_list)]

for sub in sublist:
    physio_data = new_physio_struct(tasklist)
    physio_data['subject_id'] = sub
    for task in tasklist:
        print(task)
        acqName = sub+"_"+task+".acq"
        infile = path.join(sourceDir,sub,'originals',"01+physio", acqName)
        if os.path.isfile(infile):
            print("Processing: \n\t"+infile)
            data = br.read_file(infile)
            if len(data.channels[0].data) > 1000:
                pulse, resp, hr, rr, hr_idx, rr_idx = proc_acq(data)
                physio_data[task].append(hr)
                physio_data[task].append(pulse)
                physio_data[task].append(hr_idx)
                physio_data[task].append(rr)
                physio_data[task].append(resp)
                physio_data[task].append(rr_idx)
            else:
                print("File Did not Load or was corrupt. \n\t"+acqName)
                for k in np.arange(0,6):
                    physio_data[task].append(np.nan)
        else:
            print("File Not Found: \n\t"+acqName)
            for k in np.arange(0,6):
                physio_data[task].append(np.nan)

    out_dir = os.path.join(sourceDir, sub,'QC','physio')

    if (os.path.exists(out_dir) is False):
        os.makedirs(out_dir)

    plot_subject_struct(physio_data, out_dir)
    
    if (outDir == "BIDS"):
        physio_out = os.path.join(sourceDir,outDir)
    else:
        physio_out = os.path.abspath(outDir)
        
    save_physio_tsv(physio_data, physio_out, tasklist)
