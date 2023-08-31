import csv

import matplotlib.pyplot as plt
import numpy
import matplotlib
matplotlib.use("Qt5Agg")
import numpy as np

Patterson_2010_4c_file = "/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Published/Patterson_2010_4C.csv"
Patterson_2010_1g_file = "/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Published/Patterson_2010_1G.csv"
Tanaka_2012_S3a_PSD_file = "/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Published/Tanaka_2012_S3A_PSLM.csv"
Tanaka_2012_S3a_Dend_file = "/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Published/Tanaka_2012_S3A_NON_PSLM.csv"
Tanaka_2012_exo_data = "/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Published/Tanaka_2012_exo_NPLSM.csv"

Graves_2021_glua1 = "/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Published/Graves_2021_5DGreen.csv"
Graves_2021_dsred = "/Users/surbhitwagle/Desktop/Surbhit/Work/PhD/2020/PhD/MPIBR/PhD-Project/Experimental_collab/Published/Graves_2021_5DRed.csv"


def ReadCSVFull(filename):
    """
   function to read CSV data files:
       arguments : filename : path of the .csv file
   """
    data_mat = []
    with open(filename) as csvDataFile:
        # read file as csv file
        csvReader = csv.reader(csvDataFile)
        # loop over rows
        num_rows = 0
        for row in csvReader:
            # add cell [0] to list of dates
            if row != []:
                num_rows += 1
                row_data = [float(x) for x in row]
                data_mat.append(row_data)
        assert num_rows%2 == 0
        data_mat = np.asarray(data_mat)
        # print(data_mat)
        final_data = np.zeros((num_rows//2,3))

        for i in range(num_rows//2):
            final_data[i,0] = data_mat[i,0]
            final_data[i,1] = data_mat[i,1]
            final_data[i,2]=  np.abs(data_mat[(i+num_rows//2),1] -final_data[i,1])
            # print(pos,mean,std)
            # final_data[i] = [pos,mean,sem]
        # print(final_data)
        return final_data

def ReadCSVTime(filename):
    """
   function to read CSV data files:
       arguments : filename : path of the .csv file
   """
    data_mat = []
    with open(filename,encoding='utf-8-sig') as csvDataFile:
        # read file as csv file
        csvReader = csv.reader(csvDataFile)
        # loop over rows
        num_rows = 0
        for row in csvReader:
            # add cell [0] to list of dates
            if row != []:
                num_rows += 1
                # breakpoint()
                row_data = [float(x) for x in row]
                data_mat.append(row_data)
        data_mat = np.asarray(data_mat)
        return data_mat
def patterson_data_1g():
    data = ReadCSVFull(Patterson_2010_1g_file)
    return data
def patterson_data_4c():
    data = ReadCSVFull(Patterson_2010_4c_file)
    return data

def tanaka_data_s3():
    psd_data = ReadCSVFull(Tanaka_2012_S3a_PSD_file)
    dend_data = ReadCSVFull(Tanaka_2012_S3a_Dend_file)
    return psd_data, dend_data

def graves_data_5d():
    dsred_data = ReadCSVTime(Graves_2021_dsred)
    glua1_data = ReadCSVTime(Graves_2021_glua1)
    return glua1_data,dsred_data

def get_tanaka_exo_data():
    exo_data = ReadCSVFull(Tanaka_2012_exo_data)
    return exo_data
# glua1_data,dsred_data = graves_data_5d()
# ng1d = glua1_data.T/glua1_data[:,0]
# mean_g1 = ng1d.mean(axis=1)
# tps = np.array([-10,-5,0,10,20,30,40,50,60])
# sem_g1 = ng1d.std(axis=1)/np.sqrt(glua1_data.shape[0])
# dsredd = dsred_data.T/dsred_data[:,0]
# mean_dsred = dsredd.mean(axis=1)
# sem_dsred = dsredd.std(axis=1)/np.sqrt(dsred_data.shape[0])
# # breakpoint()
# for i in range(glua1_data.shape[0]):
#     plt.plot(tps,ng1d[:,i],'k',alpha=0.2)
# plt.errorbar(tps,mean_g1,sem_g1,color="green",label="SEP-gluA1")
# # plt.errorbar(tps,mean_dsred,sem_dsred,color="red",label="dsRed")
# plt.hlines(y=1,xmin=tps[0],xmax=tps[-1])
# plt.legend()
# plt.xlabel("Time relative to induction (min)")
# plt.ylabel("Normalized signal (a.u.)")
# plt.xticks(tps)
# plt.show()


# exo_data = get_tanaka_exo_data()
# breakpoint()