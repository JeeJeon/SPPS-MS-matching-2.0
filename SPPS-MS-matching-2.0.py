"""Solid phase peptide synthesis
already know the adding sequence of amino acid
calculate all the possible mass of the peptides
on the basis of first AA was added on resin well
Version: 2.0
Author: Jee
Notice: The sequence was inputted line by line;If there is Boc/Fmoc/OMe protection,
just put them in the new line. for Boc-Phe:
Boc
Phe
Date: 2019-09-04"""

from itertools import combinations

aa=[]
dict=[]

for line in open('dict.txt'):
    temp1=line.strip('\n')
    temp2=temp1.split(', ')
    dict.append(temp2)

cf_CHNOS_s=dict[1::2]

aa=dict[0::2]

for i in range(0, len(cf_CHNOS_s)):
    for k in range(5):
        cf_CHNOS_s[i].append(float(cf_CHNOS_s[i][k]))
    del(cf_CHNOS_s[i][0:5])
for i in range(0, len(aa)):
    aa.append(aa[i][0])
del(aa[0:len(cf_CHNOS_s)])



pep=[]
pep_mass=0
num_C=0
num_H=0
num_O=0
num_N=0
num_S=0
for line in open('peptide.txt'):
    data=line.strip()
    pep.append(data)

output=open("result_peptide.txt","w")

for i in range(len(pep)):
    num_C+= cf_CHNOS_s[aa.index(pep[i])][0]
    num_H+= cf_CHNOS_s[aa.index(pep[i])][1]
    num_N+= cf_CHNOS_s[aa.index(pep[i])][2]
    num_O+= cf_CHNOS_s[aa.index(pep[i])][3]
    num_S+= cf_CHNOS_s[aa.index(pep[i])][4]
pep_mass=num_C*12.0+(num_H-2*len(pep)+2)*1.007825+(num_O-len(pep)+1)\
         *15.994915+num_N*14.003074+num_S*31.972072
print('%s\n'%str(pep))
print('%d\n'%len(pep))
print('C%dH%dN%dO%dS%d\n'%(num_C,num_H-2*len(pep)+2,num_N,num_O-len(pep)+1,num_S))
print('%.4f\n'%pep_mass)
output.write('Sequences,mass\n')
output.write('%s\n%.4f\n'%(str(pep),pep_mass))

pep_var=[]
for i in range(len(pep)-1):
    pep_var+=list(combinations(pep,i+1))


for i in range(len(pep_var)):
    num_C_v = 0
    num_H_v = 0
    num_O_v = 0
    num_N_v = 0
    num_S_v = 0
    for j in range(len(pep_var[i])):
        num_C_v+= cf_CHNOS_s[aa.index(pep_var[i][j])][0]
        num_H_v+= cf_CHNOS_s[aa.index(pep_var[i][j])][1]
        num_N_v+= cf_CHNOS_s[aa.index(pep_var[i][j])][2]
        num_O_v+= cf_CHNOS_s[aa.index(pep_var[i][j])][3]
        num_S_v+= cf_CHNOS_s[aa.index(pep_var[i][j])][4]


    pep_mass_v=num_C_v*12.0+(num_H_v-2*len(pep_var[i])+2)*1.007825+(num_O_v-len(pep_var[i])+1)\
           *15.994915+num_N_v*14.003074+num_S_v*31.972072
    output.write('%s\n%.4f'%(str(pep_var[i]),pep_mass_v))
    output.write('\n')

output.close()


"""input the sequences of peptide in peptide.txt
run 'all mass of peptides.py' to generate result_peptide.txt
which listed all of the possibilities of combinations
convert mass data into the mass_spectrum.txt
run this program to find out all the matches
Version: 0.1
Author: Jee
Date: 2019-09-05"""
from Read_txt_data import loadDatadet
from time import time

pep_mass=[]
data=[]
mass=[]
mass_match=[[]]
atom=['H','Na','K']
mass_1=[1.00728,22.98922,38.96316]
solvent=['ACN']
mass_sol=[41.0265]

for line in open('result_peptide.txt'):
    data=line.strip()
    pep_mass.append(data)


for i in range(2,len(pep_mass),2):          
    mass.append(float(pep_mass[i]))

for file_num in range(1,10):
    start=time()
    try:
       output=open("result_mass_scan_"+str(file_num)+".txt","w")
       output.write('Sequence:\tcal_mass:\tfound_mass:\tIntensity:\tdelt:\ttype:\n')
       mass_spec=loadDatadet('mass_'+str(file_num)+'.txt',3)  

       for j in range(len(mass)):
           if j==0:
               delt = 0.08  
               min_Intensity = 5 
           else:
               delt = 0.08  
               min_Intensity = 10  


           for k in range(len(mass_spec)):
                if mass_spec[k][2]>=min_Intensity:         
                    for num in range(1,5):
                        for atom1 in mass_1:
                            for num2 in range(1,4):
                                if num==1:
                                    cal1 = num2*mass[j] + atom1
                                    if cal1 <= (mass_spec[k][0] + delt) and cal1 >= \
                                    (mass_spec[k][0] - delt):
                                        output.write('%s\t' % pep_mass[2 * j + 1])
                                        output.write('%.4f\t' % cal1)
                                        output.write('%.4f\t' % mass_spec[k][0])
                                        output.write('%.2f\t' % mass_spec[k][2])
                                        output.write('%.4f\t' % (mass_spec[k][0] - cal1))
                                        if num2==1:
                                            output.write('[M+%s]+\n' % (str(atom[mass_1.index(atom1)])))
                                        else:
                                            output.write('[%dM+%s]+\n' % (num2,str(atom[mass_1.index(atom1)])))
                                    for sol in mass_sol:
                                        cal1s=num2*mass[j]+atom1+sol
                                        if cal1s <= mass_spec[k][0] + delt and cal1s >= \
                                        mass_spec[k][0] - delt:
                                            output.write('%s\t' % pep_mass[2 * j + 1])
                                            output.write('%.4f\t' % cal1s)
                                            output.write('%.4f\t' % mass_spec[k][0])
                                            output.write('%.2f\t' % mass_spec[k][2])
                                            output.write('%.4f\t' % (mass_spec[k][0] - cal1s))
                                            if num2==1:
                                                output.write('[M+%s+%s]+\n'%(str(atom[mass_1.index(atom1)]),\
                                              str(solvent[mass_sol.index(sol)])))
                                            else:
                                                output.write('[%dM+%s+%s]+\n' % (num2, str(atom[mass_1.index(atom1)]), \
                                                                                 str(solvent[mass_sol.index(sol)])))

                                else:
                                    if num==num2:
                                        continue
                                    cal2 = (num2*mass[j] + num*atom1)/num
                                    if cal2 <= (mass_spec[k][0] + delt) and cal2 >= \
                                        (mass_spec[k][0] - delt):
                                        output.write('%s\t' % pep_mass[2 * j + 1])
                                        output.write('%.4f\t' % cal2)
                                        output.write('%.4f\t' % mass_spec[k][0])
                                        output.write('%.2f\t' % mass_spec[k][2])
                                        output.write('%.4f\t' % (mass_spec[k][0] - cal2))
                                        if num2==1:
                                            output.write(
                                                '[M+%d%s]%d+\n' % (num, str(atom[mass_1.index(atom1)]), num))
                                        else:
                                            output.write('[%dM+%d%s]%d+\n' % (num2,num,str(atom[mass_1.index(atom1)]),num))
                                    for sol in mass_sol:
                                        cal2s = (num2 * mass[j] + num*atom1 + sol)/num
                                        if cal2s <= mass_spec[k][0] + delt and cal2s >= \
                                            mass_spec[k][0] - delt:
                                            output.write('%s\t' % pep_mass[2 * j + 1])
                                            output.write('%.4f\t' % cal2s)
                                            output.write('%.4f\t' % mass_spec[k][0])
                                            output.write('%.2f\t' % mass_spec[k][2])
                                            output.write('%.4f\t' % (mass_spec[k][0] - cal2s))
                                            if num2==1:
                                                output.write(
                                                    '[M+%d%s+%s]%d+\n' % (num, str(atom[mass_1.index(atom1)]), \
                                                                            str(solvent[mass_sol.index(sol)]), num))
                                            else:
                                                output.write('[%dM+%d%s+%s]%d+\n' % (num2,num,str(atom[mass_1.index(atom1)]), \
                                                                       str(solvent[mass_sol.index(sol)]),num))
                    continue
       end = time()
       print('mass_%d\nScanning accomplished!!! \n %.2fs' % (file_num, (end - start)))
    except FileNotFoundError:
        print('Can not find file:\nmass_%d!'%file_num)
    except LookupError:
        print('Unknow coding:\nmass_%d!'%file_num)
    except UnicodeDecodeError:
        print('Error in reading file:\nmass_%d!'%file_num)
    finally:
        output.close()
