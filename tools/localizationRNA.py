## castle 
import pandas as pd 

class LocalizedRNA : 

    def __init__( self, file_path : str  ) : 

        ''' represents parsed localized RNAs in an object abstracted form '''
        self.filepath = file_path 
        self.RNA_strands = None 
        self.localizationName = None # try and get localization name

        self.__readFile() # processes RNA strand from file.  

    
    def __readFile( self ) -> None : 

        # open's the file path and starts processing the files. 
        dataframe = pd.read_csv( self.filepath, usecols=['Gene_ID', 'Sequence'])

        # store sequence data in memory for processing by using a dictionary
        self.RNA_strands = dataframe.to_dict('index')


def TDD() : 

    a = LocalizedRNA('localization/out_Axon.csv') 

    print( 'ok' )

# TDD()