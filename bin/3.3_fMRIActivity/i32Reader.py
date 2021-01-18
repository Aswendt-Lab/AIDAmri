"""
Created on 10/08/2017

@author: Niklas Pallast
Neuroimaging & Neuroengineering
Department of Neurology
University Hospital Cologne

"""



import sys,os
import numpy as np
import peakdet as pk
import scipy.signal as sc



# noinspection PyTypeChecker
def getI32(file_name,numberOfSlices,numberOfAllRepitionsParTable):


    #fileSize = fileInfo.bytes/4
    fid = open(file_name)

    fileTable   = np.fromfile(fid, dtype=np.float32)


    #fpRawTime	= fileTable[0:len(fileTable):4]
    fpRawResp	= fileTable[1:len(fileTable):4]
    fpRawTrig	= fileTable[2:len(fileTable):4]
    fpRawCard	= fileTable[3:len(fileTable):4]

    # Process Respiration Data
    # Merge 10 Values
    N = 10
    RespBLC = np.convolve(fpRawResp, np.ones((N,)) / N, mode='same')

    # Evalute Baseline Shift
    RespBLC = RespBLC - np.median(RespBLC, axis=None)

    # derivative of respiration
    kernel = [1, 0, -1]
    RespDeriv = np.convolve(RespBLC, kernel, mode='same')
    #RespDeriv = np.gradient(RespBLC)
    # peak detection
    pksRespMax,pksResMin = pk.peakdet(RespBLC*20, delta=1)
    print('avg. Respiration Rate: '+str(len(pksRespMax) / (len(RespBLC) / 60000))+' 1/min')

    # Process Cardiac Data
    fs = 1000 # Sampling Frequency(1 kHz)
    lowcut = 2.5
    highcut = 10.0
    nyq = 0.5 * fs                  # Nyquist Frequency (Hz)
    Wp = [lowcut/nyq, highcut/nyq]  # Passband Frequencies (Normalised 2.5 - 10 Hz)
    Ws= [0.1/nyq, 35/nyq]           # Stopband Frequencies (Normalised)
    Rs = 40                         # Sroppband Ripple (dB)
    N = 3                           # Filter order
    b, a = sc.cheby2(N, Rs,Ws, btype='bandpass')
    filtCardBLC = sc.filtfilt(b, a, fpRawCard)


    N = 10
    CardBLC = np.convolve(filtCardBLC, np.ones((N,)) / N, mode='same')

    # Evalute Baseline Shift
    CardBLC = CardBLC - np.median(CardBLC, axis=None)

    # derivative of respiration
    kernel = [1, 0, -1]
    CardDeriv = np.convolve(CardBLC, kernel, mode='same')

    # peak detection
    pksCardpMax, pksCardMin = pk.peakdet(CardBLC * 20, delta=1)
    print('avg. Card Rate: ' + str(len(pksCardpMax) / (len(CardBLC) / 60000)) + ' 1/min')

    # if the trigger max is not equal 1 but higher
    if max(fpRawTrig) != 1.0:
        fpRawTrig = fpRawTrig-(max(fpRawTrig)-1)

    # find missing trigger and replace 1 by 0
    idx_missedTrigger = np.where(np.diff(fpRawTrig,2)==2)[0]+1
    if len(idx_missedTrigger)>0:
        fpRawTrig[idx_missedTrigger+1] =0

    triggerDataPoints = np.argwhere(fpRawTrig == 0)
    numberOfTiggers = len(triggerDataPoints)
    numberOfRepitions = numberOfTiggers / (numberOfSlices * 2)
    print('Number of Repetitions: ' + str(numberOfRepitions))

    # if two dataset in a single i32 file
    if numberOfTiggers >= ((numberOfAllRepitionsParTable + 5) * numberOfSlices * 2)*2:
        triggerDataPoints = triggerDataPoints[:int(numberOfTiggers/2)]# if more than two dataset in a single i32 file

    old_numberOfTriggers = numberOfTiggers
    # corrected number of triggers
    numberOfTiggers = len(triggerDataPoints)

    # if some wrong triggers in i32 file
    if numberOfTiggers > (numberOfAllRepitionsParTable + 5) * numberOfSlices * 2:
        wrongAmountOfTriggers = numberOfTiggers-(numberOfAllRepitionsParTable + 5) * numberOfSlices * 2
        fpRawTrig_cut = fpRawTrig[wrongAmountOfTriggers*100-1::]
        triggerDataPoints = np.argwhere(fpRawTrig_cut == 0)
        numberOfTiggers = len(triggerDataPoints)
        numberOfRepitions = numberOfTiggers / (numberOfSlices * 2)
        print('Number of Repetitions: ' + str(numberOfRepitions))



    triggerDataPoints_1st = triggerDataPoints[numberOfSlices * 5 * 2 : numberOfTiggers:2,0]
    triggerDataPoints_2nd = triggerDataPoints[numberOfSlices * 5 * 2 +1: numberOfTiggers:2,0]
    usedTriggerAmount = ((numberOfAllRepitionsParTable+5)*numberOfSlices*2-5*2*numberOfSlices)/2

    if not len(triggerDataPoints_1st) == len(triggerDataPoints_2nd):
        print('Miss one Trigger in file_name in %s', file_name)
        if len(triggerDataPoints_1st) < usedTriggerAmount:
            if len(triggerDataPoints_2nd) == usedTriggerAmount:
                triggerDataPoints_1st = triggerDataPoints_2nd
            else:
                sys.exit('Trigger does not relate to any slice or rep. Time')

    if len(RespBLC) == len(CardBLC):
        i32Table = np.zeros([len(RespBLC),4])
    else:
        sys.exit('Respiration and Cardiac Data do not have the same length!')

    i32Table[:,0] = RespBLC
    i32Table[:,1] = RespDeriv
    i32Table[:,2] = CardBLC
    i32Table[:,3] = CardDeriv

    return triggerDataPoints_1st,i32Table

if __name__ == "__main__":
    import argparse


    parser = argparse.ArgumentParser(description='I32 Reader')

    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-i','--input', help='Path to input file',required=True)

    args = parser.parse_args()


    if args.input is not None and args.input is not None:
        input = args.input
    if not os.path.exists(input):
        sys.exit("Error: '%s' is not an existing directory or file %s is not in directory." % (input, args.file,))

    result = getI32(input,16,100)
