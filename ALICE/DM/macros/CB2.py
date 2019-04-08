#! /usr/bin/env python

class CB2_pdf(MASS) :
    """Double sided Crystal Ball function with both left and rigth sides
    It appears to be very powerful and is used for many LHCb papers to describe
    B-hadron mass signals, especially for B->J/psi X final states 
    
    Note:
    - Similar to CrystalBall_pdf and unlike the original definition,
    the parameters 'n' here are shifted by 1: n <- |n| + 1
    - Typical value of parameters alpha for ``physical'' peaks is 1.5<alpha<3.0,
    - For large alpha (e.g. alpha>3), there is no sensitivity for n;
    similarly in the limit of large n, sensitivity for alpha is minimal
    """
    def __init__ ( self              ,
                   name              ,
                   xvar              , 
                   mean      = None  ,
                   sigma     = None  ,
                   alphaL    = None  ,
                   alphaR    = None  ,
                   nL        = None  ,
                   nR        = None  ) : 
        
        #
        ## initialize the base
        # 
        MASS.__init__  ( self , name , xvar , mean , sigma )
        #
        ## treat the specific parameters
        #
        self.__aL    = self.make_var ( alphaL                  ,
                                 "aL_%s"          % name ,
                                 "#alpha_{L}(%s)" % name , alphaL    , 2.0 ,  0.01 ,  5 )
        self.__nL    = self.make_var ( nL                      ,                     
                                 "nL_%s"          % name ,
                                 "n_{L}(%s)"      % name , nL        , 1   , 1.e-8 , 50 )
        self.__aR    = self.make_var ( alphaR ,
                                 "aR_%s"          % name ,
                                 "#alpha_{R}(%s)" % name , alphaR    , 2.0 , 0.01  ,  5 )
        self.__nR    = self.make_var ( nR                      ,
                                 "nR_%s"          % name ,
                                 "n_{R}(%s)"      % name , nR        , 1   , 1.e-8 , 50 )
        
        self.pdf = cpp.Analysis.Models.CrystalBallDS(
            "cb2_"       + name ,
            "CB_{2}(%s)" % name ,
            self.xvar    ,
            self.mean    ,
            self.sigma   ,
            self.aL      ,
            self.nL      ,
            self.aR      ,
            self.nR      )

        ## save the configuration
        self.config = {
            'name'   : self.name  ,
            'xvar'   : self.xvar  ,
            'mean'   : self.mean  ,
            'sigma'  : self.sigma ,
            'alphaL' : self.aL    ,
            'alphaR' : self.aR    ,
            'nL'     : self.nL    ,
            'nR'     : self.nR    ,
            }

    @property
    def aL ( self ) :
        """(left) Alpha-parameter for Crystal Ball tail"""
        return self.__aL
    @aL.setter
    def aL ( self, value ) :
        value = float ( value )
        assert 0.01 < value < 5 , "alpha_L must be between 0.01 and 5" 
        self.__aL.setVal ( value )

    @property
    def nL ( self ) :
        """(left) N-parameter for Crystal Ball tail"""
        return self.__nL
    @nL.setter
    def nL ( self, value ) :
        value = float ( value )
        assert 1.e-6 < value < 40 , "N_L must be between 1.e-6 and 40" 
        self.__nL.setVal ( value )

    @property
    def aR ( self ) :
        """(right) Alpha-parameter for Crystal Ball tail"""
        return self.__aR
    @aR.setter
    def aR ( self, value ) :
        value = float ( value )
        assert 0.01 < value < 5 , "alpha_R must be between 0.01 and 5" 
        self.__aR.setVal ( value )

    @property
    def nR ( self ) :
        """(right) N-parameter for Crystal Ball tail"""
        return self.__nR
    @nR.setter
    def nR ( self, value ) :
        value = float ( value )
        assert 1.e-6 < value < 40 , "N_R must be between 1.e-6 and 40" 
        self.__nR.setVal ( value )
        
models.append ( CB2_pdf )    

