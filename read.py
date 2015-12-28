# coding: UTF-8
import scipy.io


def read_mat_file(args):
    print "reading .mat file : " + args[1]
    x = scipy.io.loadmat(args[1])

    print "Properties of .mat data"
    for index, value in enumerate(x):
        if index == 2:
            sampleNames = x[value]
            print sampleNames[0][0]
        elif index == 3:
            cancerType = x[value]
            print "cancerType is \"" + cancerType[0] + "\""
        elif index == 4:
            types = value
            print 'types variation is \"%d\" ' % len(types)
        elif index == 5:
            strandTypes = value
            print strandTypes
        elif index == 6:
            subtypes = value
            print subtypes
        elif index == 8:
            originalGenomes = x[value]
            print value
        else:
            continue
    return originalGenomes
