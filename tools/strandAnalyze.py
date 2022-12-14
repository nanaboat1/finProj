## name : castle group 
## date : 12/13/2022
## class : cs-473-1 

import time 
from localizationRNA import LocalizedRNA

class StrandAnalysis : 

    ''' information about the strand after being processed by the various algorithms'''

    def __init__(self, strandsFile : str ) : 

        self.strand : LocalizedRNA  = LocalizedRNA( strandsFile )
        self.seqAnalysisLCS = { } 
        self.seqAnalysisLCSTrie = { }

        # run the longest common subsequence algorithm
        self.sequenceManagerLCS( ) 


    def sequenceManagerLCS( self ) : 

        ''' decides which two sequences to pass as args to lCS for processing '''
        strandlength : int = len ( self.strand.RNA_strands.keys() ) 
        
        
        if ( strandlength % 2 == 0 ) : # parity-check.

            for i in range( 2, strandlength ) : 

                # time-start 
                startTime = time.time()

                # check-s the longest common subsequence of the two strands.
                subseq = self.longestCommonSubsequence( self.strand.RNA_strands[i-2]['Sequence'], self.strand.RNA_strands[i-1]['Sequence'] ) 

                # time-end
                endTime = time.time()

                # length of longest input string.
                N , M = len( self.strand.RNA_strands[i-2]['Sequence'] ) , len ( self.strand.RNA_strands[i-1]['Sequence'] )

                inputLength : int = N*M if N > M else M*N

                execTime =  ( endTime - startTime ) # repr it in millisecs
     
                self.seqAnalysisLCS[i-1] = [ subseq, round( execTime, 2 ), inputLength ]

        else :

            for i in range( 1, strandlength ) : 

                # time-start
                startTime = time.time()

                subseq = self.longestCommonSubsequence( self.strand.RNA_strands[i-1]['Sequence'], self.strand.RNA_strands[i]['Sequence'])

                # time-end
                endTime = time.time() 

                # length of longest input string.
                N , M = len( self.strand.RNA_strands[i-1]['Sequence'] ) , len ( self.strand.RNA_strands[i]['Sequence'] )

                inputLength : int = N*M if N > M else M*N

                execTime = ( endTime - startTime ) 

                self.seqAnalysisLCS[i] = [ subseq, round( execTime, 2 ), inputLength ]


    def longestCommonSubsequence( self, strandA : str, strandB : str ) : 

        ''' a dynamic programming approach to find the lcs of two rna strands '''

        # time-complexity is O ( M N )

        N, M  = len ( strandA ), len ( strandB )

        dp = [[0 for j in range( M + 1 )] for i in range( N + 1 ) ]


        for i in range (1,  N + 1) : 

            for j in range (1, M + 1) : 

                if strandA[i-1] == strandB[j-1] : 

                    dp[i][j] = dp[i-1][j-1] + 1

                else : 
                    dp[i][j] = max( dp[i-1][j], dp[i][j-1] )
        

        return dp[N][M]



def TDD() : 

    a = StrandAnalysis() 

    b = a.longestCommonSubsequence('XMJYAUZ', 'MZJAWXU') 

    print ( b) 


# TDD()

