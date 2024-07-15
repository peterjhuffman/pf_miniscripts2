# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 15:48:50 2024

@author: HUFFMP
"""

import os
import shutil
from datetime import datetime
from time import sleep as slp
import glob

def list_files_in_directory(directory):
    return [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

def list_folders(directory_path):
    folders = []
    for item in os.listdir(directory_path):
        item_path = os.path.join(directory_path, item)
        if os.path.isdir(item_path):
            folders.append(item)
    return folders

def folderslicer(fraws):
    fraws = [f for f in fraws if len(f.split('_'))>4]
    fraws = [f for f in fraws if str(f).lower().count('wash')==0]
    
    numcaps = list(set(['_'.join(f.split('_')[:min(4, len(f.split('_')))]) for f in fraws]))

    foldproj = []
    numcounter = [len([raw for raw in fraws if raw.count(cap)>0]) for cap in numcaps]

    for cap,count in zip(numcaps,numcounter):
        if count>3:
            stocklimit=4
            bencher = fraws[[f.count(cap) for f in fraws].index(1)]
            while len([raw for raw in fraws if raw.count('_'.join(bencher.split('_')[:min(stocklimit, len(bencher.split('_')))]))>0])==count:
                stocklimit+=1
            foldproj.append('_'.join(bencher.split('_')[:min(stocklimit-1, len(bencher.split('_')))]))
    return foldproj


def file_in_dir(file, directory):
    for root,dirs,files in os.walk(directory):
        for name in files:
            if name == file:
                return True
    return False


# '_'.join(bencher.split('_')[:min(stocklimit, len(bencher.split('_')))])

# Use the function
# directory = '\\\\amer.pfizer.com\\pfizerfiles\\Research\\LAJ\\oru-tcb\\ChemBio\\HUFFMP\\misc\\meth'  # replace with your directory
# files = list_files_in_directory(directory)
# # print(files)


# for file in files:

    # specify file path


months = {1:'January',
          2:'February',
          3:'March',
          4:'April',
          5:'May',
          6:'June',
          7:'July',
          8:'August',
          9:'September',
          10:'October',
          11:'November',
          12:'December'}

input('Press any key to start.')
mach = input('Press E for eclipse, A for ascend:').lower()
switch=False
if mach=='e':
    machine = 'Eclipse'
elif mach=='a':
    machine='Ascend'
else:
    switch=True

while switch:
    mach = input('Press E for eclipse, A for ascend:').lower()
    switch=False
    if mach=='e':
        machine = 'Eclipse'
    elif mach=='a':
        machine='Ascend'
    else:
        switch=True



pull = True

while pull:
    cpath = 'D:\\'
    ds  = list_files_in_directory(cpath)
    search_dir = cpath

    # Get a list of all files (excluding directories and symlinks)
    files = list(filter(os.path.isfile, glob.glob(search_dir + "/*")))

    # Sort the files based on modification time
    files.sort(key=lambda x: os.path.getmtime(x))

    c = [str(cnew.split('\\')[-1]) for cnew in files if str(cnew)[-4:]=='.raw'][:-1]
    print(c)
    # path2 = 
    lf = list_folders('\\\\cifs.lajvfs100.pfizer.com\\oru-tcb\\ChemBio\\MS_RawFiles\\')

    for file in c:
        file = file[1:]
        if file.split('.')[-1]=='raw':
            if str(file).lower().count('wash')==0:
                print(file)
                ts = os.path.getctime(cpath+file)
                dt = datetime.fromtimestamp(ts)

                target_path = f'\\\\cifs.lajvfs100.pfizer.com\\oru-tcb\\ChemBio\\MS_RawFiles\\{dt.year}_{machine}_runs\\{dt.month}_{months[dt.month]}_{dt.year}\\'

                source_file = cpath+f'{file}'
                destination_folder = target_path
                if not os.path.exists(target_path):
                    os.makedirs(target_path)

                # get the base name of the source file
                file_name = os.path.basename(source_file)

                # construct the destination file path
                destination_file = os.path.join(destination_folder, file_name)

                # check if the file exists at the source and not at the destination
                if os.path.isfile(source_file) and not file_in_dir(file, destination_folder):
                    try:
                        # use shutil.move() to move the file
                        shutil.copy(source_file, destination_folder)
                        print(f"Copying {file} ===================================== >>>>>>>>>>>>>>>>>")
                        pass
                    except Exception as e:
                        print(f"Error occurred while moving file: {e}")
                else:
                    if not os.path.isfile(source_file):
                        print(f"The file {source_file} does not exist.")
                    if os.path.isfile(destination_file):
                        print(f"A file with the same name already exists at the destination {destination_file}.")




    path1 = '\\\\cifs.lajvfs100.pfizer.com\\oru-tcb\\ChemBio\\MS_RawFiles\\'
    lf = list_folders(path1)

    for fold in lf:
        lfquad = list_folders(path1+fold)
        for nfold in lfquad:
            print(f"{fold} - {nfold}")
            fstack = list_files_in_directory(f"{path1}{fold}\\{nfold}\\")
            fraws = [f for f in fstack if str(f).lower().split('.')[-1]=='raw']

            fraws_bench = sum([(str(f).lower().count('hela')+str(f).lower().count('bsa')) for f in fraws])>0
            # print(f"{fold},{nfold}   ",fraws_bench)

            if fraws_bench:
                destination_folder = path1+fold+f"\\{nfold}\\standards"
                if not os.path.exists(destination_folder):
                    os.makedirs(destination_folder)

                for f in [f for f in fraws if (str(f).lower().count('hela')+str(f).lower().count('bsa'))>0]:
                    print(f)
                    source_file = f"{path1}{fold}\\{nfold}\\{f}"
                    destination_file = f"{destination_folder}\\{f}"
                    # check if the file exists at the source and not at the destination
                    if os.path.isfile(source_file) and not os.path.isfile(destination_file):
                        try:
                            # use shutil.move() to move the file
                            shutil.move(source_file, destination_folder)
                            # print(f"{source_file} --> {destination_folder}")
                            pass
                        except Exception as e:
                            print(f"Error occurred while moving file: {e}")
                    else:
                        if not os.path.isfile(source_file):
                            print(f"The file {source_file} does not exist.")
                        if os.path.isfile(destination_file):
                            print(f"A file with the same name already exists at the destination {destination_file}.")


            fstack = list_files_in_directory(f"{path1}{fold}\\{nfold}\\")
            fraws = [f for f in fstack if str(f).lower().split('.')[-1]=='raw']

            if len(fraws)>3:
                new_folds = folderslicer(fraws)

                for foldi in new_folds:
                    destination_folder = path1+fold+f"\\{nfold}\\{foldi}"
                    # print(destination_folder)
                    if not os.path.exists(destination_folder):
                        os.makedirs(destination_folder)

            for fraw in fraws:
                ffold = list_folders(f"{path1}{fold}\\{nfold}\\")
                ffold.sort(key=len)
                ffold.reverse()
                
                for ff in ffold:
                    if fraw.count(ff)>0:
                        destination_folder = f"{path1}{fold}\\{nfold}\\{ff}\\"
                        source_file = f"{path1}{fold}\\{nfold}\\{fraw}"
                        destination_file = f"{path1}{fold}\\{nfold}\\{ff}\\{fraw}"
                        # check if the file exists at the source and not at the destination
                        if os.path.isfile(source_file) and not os.path.isfile(destination_file):
                            try:
                                # use shutil.move() to move the file
                                shutil.move(source_file, destination_folder)
                                print(fraw)
                                # pass
                                # print(f"{source_file} --> {destination_folder}")
                                pass
                            except Exception as e:
                                print(f"Error occurred while moving file: {e}")
                        else:
                            if not os.path.isfile(source_file):
                                print(f"The file {source_file} does not exist.")
                            if os.path.isfile(destination_file):
                                print(f"A file with the same name already exists at the destination {destination_file}.")
                        break
    
    print('Sleeping.')
    slp(60*60)
    # pull = False
input('Done.')
























