#!/usr/bin/env python
class Vb:

    def __init__(self, data):
        self.data = data

    def show(self):
        print self.data

data = 'aaa'
a = Vb(data)
a.show()

labels = []
genomes = []
nums = []
motos = []
tikans = []

for line in open('Nik_Zainal_2012.mutationPositionFormat.txt', 'r'):
    itemList = line[:-1].split('\t')
    for index, value in enumerate(itemList):
        if index == 0:
            labels.append(value)
        if index == 1:
            genomes.append(value)
        if index == 2:
            nums.append(value)
        if index == 3:
            motos.append(value)
        if index == 4:
            tikans.append(value)

label_set = set(labels)
genome_set = set(genomes)
num_set = set(nums)
moto_set = set(motos)
tikan_set = set(tikans)


def replace2point(value, target_list, target_set, linenum):
    for point, label in enumerate(target_set):
        if value == label:
            labels[linenum] = point
        else:
            print point
    return target_list

for linenum, line in enumerate(open('Nik_Zainal_2012.mutationPositionFormat.txt', 'r')):
    itemList = line[:-1].split('\t')
    for index, value in enumerate(itemList):
        if index == 0:
            for point, label in enumerate(label_set):
                if value == label:
                    labels[linenum] = point
                    break
                else:
                    continue
        elif index == 1:
            for point, label in enumerate(genome_set):
                if value == label:
                    genomes[linenum] = point
                    break
                else:
                    continue
        elif index == 2:
            continue
            """
            for point, label in enumerate(num_set):
                if value == label:
                    nums[linenum] = point
                    break
                else:
                    continue
            """
        elif index == 3:
            for point, label in enumerate(moto_set):
                if value == label:
                    motos[linenum] = point
                    break
                else:
                    continue
        elif index == 4:
            for point, label in enumerate(tikan_set):
                if value == label:
                    tikans[linenum] = point
                    break
                else:
                    continue
        else:
            print "whats!?"
