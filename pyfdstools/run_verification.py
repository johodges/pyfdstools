# -*- coding: utf-8 -*-
"""
Created on Mon Sep 15 08:11:18 2025

@author: jhodges
"""

import os, subprocess, platform, argparse, sys

def executeModel(inputs,showTerminal=False,printOutput=False):
    filename = inputs[0]
    systemDir = inputs[1]
    executable = inputs[2]
    path = inputs[3]
    mpiprocesses = '%d'%(inputs[4])
    console = []
    logname = os.path.basename(filename)+"_log.err"
    #cmd = ["start","cmd","/k",executable, filename, '>&', 'log.err']
    #cmd = ["start","cmd","/k",executable, filename]
    if platform.system() == 'Windows':
        if showTerminal:
            cmd = ["start","cmd","/k", 'mpiexec', '-np', mpiprocesses, executable, filename]
        else:
            cmd = ['mpiexec', '-np', mpiprocesses, executable, filename] #, '2>&1', logname]
    else:
        cmd = ['mpiexec', '-np', mpiprocesses, executable, filename, '>&', logname]
    
    my_env = os.environ.copy()
    if path is not False:
        my_env['PATH']=path
    if showTerminal:
        stdout = subprocess.PIPE
        stderr = subprocess.PIPE
        shell = True
    else:
        stdout=subprocess.DEVNULL
        stderr=subprocess.STDOUT
        shell = True
    #print(path, my_env['PATH'], shell)
    #print(cmd)
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=shell, cwd=systemDir, env=my_env)
    
    if showTerminal:
        p.wait()
    else:
        while p.stdout.readable():
            line = p.stdout.readline()
            console.append(line)
            if not line:
                break
            if printOutput:
                print(line.strip().decode('utf-8'))
    return console

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='Run verification cases',
                                     description='Example routine to run a verification case',
                                     epilog='Help goes here')
    parser.add_argument('--executable', help='path to fds executable to run', type=str, default="C:\\firemodels\\fds_master\\Build\\impi_intel_win\\fds_impi_intel_win.exe")
    parser.add_argument('--fdssource', help='path to firemodels directory', type=str, default="C:\\firemodels\\fds_master")
    parser.add_argument('--path', help='environment variable PATH for run within subprocess', type=str, default="C:\\Program Files (x86)\\Intel\\oneAPI\\mpi\\latest\\env\\..\\opt\\mpi\\libfabric\\bin;C:\\Program Files (x86)\\Intel\\oneAPI\\mpi\\latest\\env\\..\\bin;C:\\firemodels\\fds_master\\Build\\impi_intel_win")
    parser.add_argument('--cases', help='cases to run', type=int, nargs='*', default=[89])
    
    args = parser.parse_args()
    
    executable = args.executable
    path = args.path
    fdssource = os.path.abspath(args.fdssource)
    
    verification_bash = os.path.join(fdssource,'Verification','FDS_Cases.sh')
    
    with open(verification_bash, 'r') as f:
        verification_txt = f.readlines()
    
    if args.cases != None:
        cases = args.cases
        for i in range(0, len(cases)):
            txt = verification_txt[cases[i]-1]
            indir = txt.split('-d ')[1].split(' ')[0]
            filename = txt.split(' ')[-1]
            if '-p' in txt:
                mpiprocesses = int(txt.split('-p ')[1].split(' ')[0])
            else:
                mpiprocesses = 1
            txt = txt.replace('$QFDS','')
            systemDir = os.path.join(fdssource,'Verification',indir) + os.sep
            inputs = [filename, systemDir, executable, path, mpiprocesses]
            console = executeModel(inputs,showTerminal=False, printOutput=True)
            #print(mpiprocesses, indir, filename)
    