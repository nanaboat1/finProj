## name : castle group 
## date : 12/13/2022
## class : cs-473-1 


# performs data analysis of localized strands of RNA. 
import os 
from strandAnalyze import StrandAnalysis 
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd 
import openpyxl 


class RNAAlysis : 

    def __init__( self )  : 

        # get the file-paths of all the localized files.
        self.localized_names = os.scandir('localization/')
        self.analysisresult = {}
        self.graphdata= [ ]
        self.graphdataExcel = [ ]
    
    def Analysis( self ) : 

        ''' performs time analysis on different RNA strands of varying lengths. '''

        for filepath in self.localized_names : 

            localizedAnalysis : StrandAnalysis = StrandAnalysis('localization/'+filepath.name)
 
            self.analysisresult[filepath.name] = localizedAnalysis.seqAnalysisLCS


        # Debugging purposes. 
        # Atestanalysis : StrandAnalysis = StrandAnalysis('localization/out_Axon.csv')

        # Bt  = StrandAnalysis('localization/out_Plasma membrane vesicle.csv')

        # self.analysisresult['out_Axon.csv'] = Atestanalysis.seqAnalysisLCS

        # self.analysisresult['out_Plasma membrane'] = Bt.seqAnalysisLCS

        
    
    def dataPrepAnalysis ( self ) : 

        # when we get the analysis result, it has strand length and strand name, then we find the avg of each localization. 
        for localeName, resultObj in self.analysisresult.items() : 

            for idx, results in resultObj.items() : 
                self.graphdata.append( [results[2], results[1]] ) # append the length of strand and its corresponding time to lcs

        # sort the data based on lowest time to highest time. 
        self.graphdata.sort( key= lambda x:x[1])

        


    def graphAnalysis ( self ) : 

        inputX = [ x[1] for x in self.graphdata ]
        inputY = [ y[0] for y in self.graphdata ]

        # write data to excel. 
        dataframe = pd.DataFrame ( self.graphdata, columns=['Length of RNA Strand', 'Longest Subsequence Dynamic Programming'] )

        try :

            dataframe.to_excel('castle_LCS_Analysis2.xlsx')
        
        except : 
            print(' dont write to excel now')
        
    

        X, Y = np.array( inputX ), np.array( inputY )

        # computing line of best-fit. 
        S, C = np.polyfit( X, Y, 1) 

        # display graph
        plt.scatter( X, Y, color='purple') 

        plt.title(" Empirical Analysis of Longest Common Subsequence ( Dynamic Programming Approach ) ", loc='center')
        plt.xlabel(" Algorithm execution time (s) ")
        plt.ylabel( " Length of input String N (#) ")

        # add line of best fit to plot. to have a concept of the algorithm.
        plt.plot(X, S*X+C, color='blue', linestyle='--', linewidth=2 )
        plt.show()




def TDD( ) : 

    a = RNAAlysis() 

    a.Analysis()

    a.dataPrepAnalysis()

    a.graphAnalysis()

    print('ok')


TDD() 




