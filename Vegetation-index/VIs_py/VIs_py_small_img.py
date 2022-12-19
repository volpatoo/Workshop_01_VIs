# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 12:08:52 2022

@author: Leonardo Volpato
@email: volpato1@msu.edu
"""

## Check RAM available
from psutil import virtual_memory
ram_gb = virtual_memory().total / 1e9
print('Your runtime has {:.2f} gigabytes of available RAM\n'.format(ram_gb))

if ram_gb < 20:
  print('You should consider increasing the memory RAM before to start')
else:
  print('You are using a high-RAM runtime!')


# Import Libraries 

import numpy as np
import numpy.ma as ma
#from matplotlib import pyplot as plt
import rasterio
import geopandas as gpd
import rasterstats as rs
#from IPython.display import display
import pandas as pd
import os
import fiona
import rasterio.mask
#import matplotlib.pyplot as plt
import warnings
import rasterio
import rasterio.features
warnings.filterwarnings('ignore') #don't display warnings

# Class implemented to calculus the index
# Vegetation indices contribuetions from: 
    # [João Gustavo A. Amorim] (https://github.com/TheAlgorithms/Python/blob/master/digital_image_processing/index_calculation.py)

class IndexCalculation:
    """
    # Class Summary
            This algorithm consists in calculating vegetation indices, these
        indices can be used for precision agriculture for example (or remote
        sensing). There are functions to define the data and to calculate the
        implemented indices.
    # Vegetation index
        https://en.wikipedia.org/wiki/Vegetation_Index
        A Vegetation Index (VI) is a spectral transformation of two or more bands
        designed to enhance the contribution of vegetation properties and allow
        reliable spatial and temporal inter-comparisons of terrestrial
        photosynthetic activity and canopy structural variations
    # Information about channels (Wavelength range for each)
        * nir - near-infrared
            https://www.malvernpanalytical.com/br/products/technology/near-infrared-spectroscopy
            Wavelength Range 700 nm to 2500 nm
        * Red Edge
            https://en.wikipedia.org/wiki/Red_edge
            Wavelength Range 680 nm to 730 nm
        * red
            https://en.wikipedia.org/wiki/Color
            Wavelength Range 635 nm to 700 nm
        * blue
            https://en.wikipedia.org/wiki/Color
            Wavelength Range 450 nm to 490 nm
        * green
            https://en.wikipedia.org/wiki/Color
            Wavelength Range 520 nm to 560 nm
    # Implemented index list
            #"abbreviationOfIndexName" -- list of channels used
            #"ARVI2"            --  red, nir
            #"CCCI"             --  red, redEdge, nir
            #"CVI"              --  red, green, nir
            #"GLI"              --  red, green, blue
            #"NDVI"             --  red, nir
            #"BNDVI"            --  blue, nir
            #"redEdgeNDVI"      --  red, redEdge
            #"GNDVI"            --  green, nir
            #"GBNDVI"           --  green, blue, nir
            #"GRNDVI"           --  red, green, nir
            #"RBNDVI"           --  red, blue, nir
            #"PNDVI"            --  red, green, blue, nir
            #"ATSAVI"           --  red, nir
            #"BWDRVI"           --  blue, nir
            #"CIgreen"          --  green, nir
            #"CIrededge"        --  redEdge, nir
            #"CI"               --  red, blue
            #"CTVI"             --  red, nir
            #"GDVI"             --  green, nir
            #"EVI"              --  red, blue, nir
            #"GEMI"             --  red, nir
            #"GOSAVI"           --  green, nir
            #"GSAVI"            --  green, nir
            #"HUE"              --  red, green, blue
            #"IVI"              --  red, nir
            #"IPVI"             --  red, nir
            #"I"                --  red, green, blue
            #"RVI"              --  red, nir
            #"MRVI"             --  red, nir
            #"MSAVI"            --  red, nir
            #"NormG"            --  red, green, nir
            #"NormNIR"          --  red, green, nir
            #"NormR"            --  red, green, nir
            #"NGRDI"            --  red, green
            #"RI"               --  red, green
            #"IF"               --  red, green, blue
            #"DVI"              --  red, nir
            #"TVI"              --  red, nir
            #"NDRE"             --  redEdge, nir
            #"BI"               --  red, green, blue
            #"BIM"              --  red, green, blue
            #"HI"               --  red, green, blue
            #"SI"               --  red, blue
            #"VARI"             --  red, green, blue
            #"BdivG"            --  green, blue
            #"BCC"              --  red, green, blue
            #"CIVE"             --  red, green, blue
            #"COM1"             --  red, green, blue
            #"COM2"             --  red, green, blue
            #"ExG"              --  red, green, blue
            #"ExG2"             --  red, green, blue
            #"ExGR"             --  red, green, blue
            #"EXR"              --  red, green
            #"GmnB"             --  green, blue
            #"GmnR"             --  red, green
            #"GdivB"            --  green, blue
            #"GdivR"            --  red, green
            #"GCC"              --  red, green, blue
            #"MExG"             --  red, green, blue
            #"MGVRI"            --  red, green, blue
            #"NDI"              --  red, green
            #"NDRBI"            --  red, blue
            #"NGBDI"            --  green, blue
            #"RmnB"             --  red, blue
            #"RdiB"             --  red, blue
            #"RCC"              --  red, green, blue
            #"MRCCbyAlper"      --  red, green, blue
            #"RGBVI"            --  red, green, blue
            #"TGI"              --  red, green, blue
            #"VEG"              --  red, green, blue
            #"MyIndexi"         --  red, green, blue
            #"MSRGR"            --  red, green
            #"TNDGR"            --  red, green, blue
                       
           
    #list of all index implemented
        #allIndex = ["ARVI2","ATSAVI","BCC","BdivG","BI","BIM","BNDVI","BWDRVI",
        "CCCI","CI","CIgreen","CIrededge","CIVE","COM1","COM2","CTVI","CVI","DVI",
        "EVI","ExG","ExG2","ExGR","EXR","GBNDVI","GCC","GdivB","GdivR","GDVI","GEMI",
        "GLI","GmnB","GmnR","GNDVI","GOSAVI","GRNDVI","GSAVI","HI","HUE","I","IF","IPVI",
        "IVI","MExG","MGVRI","MRCCbyAlper","MRVI","MSAVI","MSRGR","MyIndexi","NDI","NDRBI",
        "NDRE","NDVI","NGBDI","NGRDI","NormG","NormNIR","NormR","PNDVI","RBNDVI","RCC","RdiB",
        "redEdgeNDVI","RGBVI","RI","RmnB","RVI","SI","TGI","TNDGR","TVI","VARI","VEG"]
    
    #list of index with not blue channel
        #notBlueIndex = ["ARVI2", "CCCI", "CVI", "NDVI", "redEdgeNDVI", "GNDVI",
                         "GRNDVI", "ATSAVI", "CIgreen", "CIrededge", "CTVI", "GDVI",
                         "GEMI", "GOSAVI", "GSAVI", "IVI", "IPVI", "RVI", "MRVI",
                         "MSAVI", "NormG", "NormNIR", "NormR", "NGRDI", "RI", "DVI",
                         "TVI", "NDRE"]
    #list of index just with RGB channels
        #RGBIndex = ["GLI", "CI", "HUE", "I", "NGRDI", "RI", "IF",
                    'BI','BIM','HI','NGRDI','SI','VARI','BdivG',
                    'BCC','CIVE','COM1','COM2','ExG','ExG2','ExGR','EXR',
                    'GmnB','GmnR','GdivB','GdivR','GCC','MExG','MGVRI',
                    'NDI','NDRBI','NGBDI','RmnB','RdiB','RCC','MRCCbyAlper',
                    'RGBVI','TGI','VEG','MyIndexi','MSRGR','TNDGR']
    """

    def __init__(self, red=None, green=None, blue=None, redEdge=None, nir=None):
        # print("Numpy version: " + np.__version__)
        self.setMatrices(red=red, green=green, blue=blue, redEdge=redEdge, nir=nir)

    def setMatrices(self, red=None, green=None, blue=None, redEdge=None, nir=None):
        if red is not None:
            self.red = red
        if green is not None:
            self.green = green
        if blue is not None:
            self.blue = blue
        if redEdge is not None:
            self.redEdge = redEdge
        if nir is not None:
            self.nir = nir
        return True

    def calculation(
        self, index="", red=None, green=None, blue=None, redEdge=None, nir=None
    ):
        """
        performs the calculation of the index with the values instantiated in the class
        :str index: abbreviation of index name to perform
        """
        self.setMatrices(red=red, green=green, blue=blue, redEdge=redEdge, nir=nir)
        funcs = {
            "ARVI2": self.ARVI2,
            "CCCI": self.CCCI,
            "CVI": self.CVI,
            "GLI": self.GLI,
            "NDVI": self.NDVI,
            "BNDVI": self.BNDVI,
            "redEdgeNDVI": self.redEdgeNDVI,
            "GNDVI": self.GNDVI,
            "GBNDVI": self.GBNDVI,
            "GRNDVI": self.GRNDVI,
            "RBNDVI": self.RBNDVI,
            "PNDVI": self.PNDVI,
            "ATSAVI": self.ATSAVI,
            "BWDRVI": self.BWDRVI,
            "CIgreen": self.CIgreen,
            "CIrededge": self.CIrededge,
            "CI": self.CI,
            "CTVI": self.CTVI,
            "GDVI": self.GDVI,
            "EVI": self.EVI,
            "GEMI": self.GEMI,
            "GOSAVI": self.GOSAVI,
            "GSAVI": self.GSAVI,
            "HUE": self.HUE,
            "IVI": self.IVI,
            "IPVI": self.IPVI,
            "I": self.I,
            "RVI": self.RVI,
            "MRVI": self.MRVI,
            "MSAVI": self.MSAVI,
            "NormG": self.NormG,
            "NormNIR": self.NormNIR,
            "NormR": self.NormR,
            "NGRDI": self.NGRDI,
            "RI": self.RI,
            "IF": self.IF,
            "DVI": self.DVI,
            "TVI": self.TVI,
            "NDRE": self.NDRE,
            "BI":self.BI,
            "BIM":self.BIM,
            "HI":self.HI,
            "SI":self.SI,
            "VARI":self.VARI,
            "BdivG":self.BdivG,
            "BCC":self.BCC,
            "CIVE":self.CIVE,
            "COM1":self.COM1,
            "COM2":self.COM2,
            "ExG":self.ExG,
            "ExG2":self.ExG2,
            "ExGR":self.ExGR,
            "EXR":self.EXR,
            "GmnB":self.GmnB,
            "GmnR":self.GmnR,
            "GdivB":self.GdivB,
            "GdivR":self.GdivR,
            "GCC":self.GCC,
            "MExG":self.MExG,
            "MGVRI":self.MGVRI,
            "NDI":self.NDI,
            "NDRBI":self.NDRBI,
            "NGBDI":self.NGBDI,
            "RmnB":self.RmnB,
            "RdiB":self.RdiB,
            "RCC":self.RCC,
            "MRCCbyAlper":self.MRCCbyAlper,
            "RGBVI":self.RGBVI,
            "TGI":self.TGI,
            "VEG":self.VEG,
            "MyIndexi":self.MyIndexi,
            "MSRGR":self.MSRGR,
            "TNDGR":self.TNDGR,

        }

        try:
            return funcs[index]()
        except KeyError:
            print("Index not in the list!")
            return False

    def ARVI2(self):
        """
        Atmospherically Resistant Vegetation Index 2
        https://www.indexdatabase.de/db/i-single.php?id=396
        :return: index
            −0.18+1.17*(self.nir−self.red)/(self.nir+self.red)
        """
        return -0.18 + (1.17 * ((self.nir - self.red) / (self.nir + self.red)))

    def CCCI(self):
        """
        Canopy Chlorophyll Content Index
        https://www.indexdatabase.de/db/i-single.php?id=224
        :return: index
        """
        return ((self.nir - self.redEdge) / (self.nir + self.redEdge)) / (
            (self.nir - self.red) / (self.nir + self.red)
        )

    def CVI(self):
        """
        Chlorophyll vegetation index
        https://www.indexdatabase.de/db/i-single.php?id=391
        :return: index
        """
        return self.nir * (self.red / (self.green**2))

    def GLI(self):
        """
        self.green leaf index
        https://www.indexdatabase.de/db/i-single.php?id=375
        :return: index
        """
        return (2 * self.green - self.red - self.blue) / (
            2 * self.green + self.red + self.blue
        )

    def NDVI(self):
        """
        Normalized Difference self.nir/self.red Normalized Difference Vegetation
        Index, Calibrated NDVI - CDVI
        https://www.indexdatabase.de/db/i-single.php?id=58
        :return: index
        """
        return (self.nir - self.red) / (self.nir + self.red)

    def BNDVI(self):
        """
            Normalized Difference self.nir/self.blue self.blue-normalized difference
        vegetation index
        https://www.indexdatabase.de/db/i-single.php?id=135
        :return: index
        """
        return (self.nir - self.blue) / (self.nir + self.blue)

    def redEdgeNDVI(self):
        """
        Normalized Difference self.rededge/self.red
        https://www.indexdatabase.de/db/i-single.php?id=235
        :return: index
        """
        return (self.redEdge - self.red) / (self.redEdge + self.red)

    def GNDVI(self):
        """
        Normalized Difference self.nir/self.green self.green NDVI
        https://www.indexdatabase.de/db/i-single.php?id=401
        :return: index
        """
        return (self.nir - self.green) / (self.nir + self.green)

    def GBNDVI(self):
        """
        self.green-self.blue NDVI
        https://www.indexdatabase.de/db/i-single.php?id=186
        :return: index
        """
        return (self.nir - (self.green + self.blue)) / (
            self.nir + (self.green + self.blue)
        )

    def GRNDVI(self):
        """
        self.green-self.red NDVI
        https://www.indexdatabase.de/db/i-single.php?id=185
        :return: index
        """
        return (self.nir - (self.green + self.red)) / (
            self.nir + (self.green + self.red)
        )

    def RBNDVI(self):
        """
        self.red-self.blue NDVI
        https://www.indexdatabase.de/db/i-single.php?id=187
        :return: index
        """
        return (self.nir - (self.blue + self.red)) / (self.nir + (self.blue + self.red))

    def PNDVI(self):
        """
        Pan NDVI
        https://www.indexdatabase.de/db/i-single.php?id=188
        :return: index
        """
        return (self.nir - (self.green + self.red + self.blue)) / (
            self.nir + (self.green + self.red + self.blue)
        )

    def ATSAVI(self, X=0.08, a=1.22, b=0.03):
        """
        Adjusted transformed soil-adjusted VI
        https://www.indexdatabase.de/db/i-single.php?id=209
        :return: index
        """
        return a * (
            (self.nir - a * self.red - b)
            / (a * self.nir + self.red - a * b + X * (1 + a**2))
        )

    def BWDRVI(self):
        """
        self.blue-wide dynamic range vegetation index
        https://www.indexdatabase.de/db/i-single.php?id=136
        :return: index
        """
        return (0.1 * self.nir - self.blue) / (0.1 * self.nir + self.blue)

    def CIgreen(self):
        """
        Chlorophyll Index self.green
        https://www.indexdatabase.de/db/i-single.php?id=128
        :return: index
        """
        return (self.nir / self.green) - 1

    def CIrededge(self):
        """
        Chlorophyll Index self.redEdge
        https://www.indexdatabase.de/db/i-single.php?id=131
        :return: index
        """
        return (self.nir / self.redEdge) - 1

    def CI(self):
        """
        Coloration Index
        https://www.indexdatabase.de/db/i-single.php?id=11
        :return: index
        """
        return (self.red - self.blue) / self.red

    def CTVI(self):
        """
        Corrected Transformed Vegetation Index
        https://www.indexdatabase.de/db/i-single.php?id=244
        :return: index
        """
        ndvi = self.NDVI()
        return ((ndvi + 0.5) / (abs(ndvi + 0.5))) * (abs(ndvi + 0.5) ** (1 / 2))

    def GDVI(self):
        """
        Difference self.nir/self.green self.green Difference Vegetation Index
        https://www.indexdatabase.de/db/i-single.php?id=27
        :return: index
        """
        return self.nir - self.green

    def EVI(self):
        """
        Enhanced Vegetation Index
        https://www.indexdatabase.de/db/i-single.php?id=16
        :return: index
        """
        return 2.5 * (
            (self.nir - self.red) / (self.nir + 6 * self.red - 7.5 * self.blue + 1)
        )

    def GEMI(self):
        """
        Global Environment Monitoring Index
        https://www.indexdatabase.de/db/i-single.php?id=25
        :return: index
        """
        n = (2 * (self.nir**2 - self.red**2) + 1.5 * self.nir + 0.5 * self.red) / (
            self.nir + self.red + 0.5
        )
        return n * (1 - 0.25 * n) - (self.red - 0.125) / (1 - self.red)

    def GOSAVI(self, Y=0.16):
        """
        self.green Optimized Soil Adjusted Vegetation Index
        https://www.indexdatabase.de/db/i-single.php?id=29
        mit Y = 0,16
        :return: index
        """
        return (self.nir - self.green) / (self.nir + self.green + Y)

    def GSAVI(self, L=0.5):
        """
        self.green Soil Adjusted Vegetation Index
        https://www.indexdatabase.de/db/i-single.php?id=31
        mit L = 0,5
        :return: index
        """
        return ((self.nir - self.green) / (self.nir + self.green + L)) * (1 + L)

    def HUE(self):
        """
        HUE
        https://www.indexdatabase.de/db/i-single.php?id=34
        :return: index
        """
        return np.arctan(2 *(self.blue - self.green - self.red) / 30.5 * (self.green - self.red))
            
       
    def IVI(self, a=None, b=None):
        """
        Ideal vegetation index
        https://www.indexdatabase.de/db/i-single.php?id=276
        b=intercept of vegetation line
        a=soil line slope
        :return: index
        """
        return (self.nir - b) / (a * self.red)

    def IPVI(self):
        """
        Infraself.red percentage vegetation index
        https://www.indexdatabase.de/db/i-single.php?id=35
        :return: index
        """
        return (self.nir / ((self.nir + self.red) / 2)) * (self.NDVI() + 1)

    def I(self):  # noqa: E741,E743
        """
        Intensity
        https://www.indexdatabase.de/db/i-single.php?id=36
        :return: index
        """
        return (self.red + self.green + self.blue) / 30.5

    def RVI(self):
        """
        Ratio-Vegetation-Index
        http://www.seos-project.eu/modules/remotesensing/remotesensing-c03-s01-p01.html
        :return: index
        """
        return self.nir / self.red

    def MRVI(self):
        """
        Modified Normalized Difference Vegetation Index RVI
        https://www.indexdatabase.de/db/i-single.php?id=275
        :return: index
        """
        return (self.RVI() - 1) / (self.RVI() + 1)

    def MSAVI(self):
        """
        Modified Soil Adjusted Vegetation Index
        https://www.indexdatabase.de/db/i-single.php?id=44
        :return: index
        """
        return (
            (2 * self.nir + 1)
            - ((2 * self.nir + 1) ** 2 - 8 * (self.nir - self.red)) ** (1 / 2)
        ) / 2

    def NormG(self):
        """
        Norm G
        https://www.indexdatabase.de/db/i-single.php?id=50
        :return: index
        """
        return self.green / (self.nir + self.red + self.green)

    def NormNIR(self):
        """
        Norm self.nir
        https://www.indexdatabase.de/db/i-single.php?id=51
        :return: index
        """
        return self.nir / (self.nir + self.red + self.green)

    def NormR(self):
        """
        Norm R
        https://www.indexdatabase.de/db/i-single.php?id=52
        :return: index
        """
        return self.red / (self.nir + self.red + self.green)

    def NGRDI(self):
        """
            Normalized Difference self.green/self.red Normalized self.green self.red
        difference index, Visible Atmospherically Resistant Indices self.green
        (VIself.green)
        https://www.indexdatabase.de/db/i-single.php?id=390
        :return: index
        """
        return (self.green - self.red) / (self.green + self.red)

    def RI(self):
        """
        Normalized Difference self.red/self.green self.redness Index
        https://www.indexdatabase.de/db/i-single.php?id=74
        :return: index
        """
        return (self.red - self.green) / (self.red + self.green)

    def IF(self):
        """
        Shape Index
        https://www.indexdatabase.de/db/i-single.php?id=79
        :return: index
        """
        return (2 * self.red - self.green - self.blue) / (self.green - self.blue)

    def DVI(self):
        """
        Simple Ratio self.nir/self.red Difference Vegetation Index, Vegetation Index
        Number (VIN)
        https://www.indexdatabase.de/db/i-single.php?id=12
        :return: index
        """
        return self.nir / self.red

    def TVI(self):
        """
        Transformed Vegetation Index
        https://www.indexdatabase.de/db/i-single.php?id=98
        :return: index
        """
        return (self.NDVI() + 0.5) ** (1 / 2)

    def NDRE(self):
        return (self.nir - self.redEdge) / (self.nir + self.redEdge)


    def BI(self):
        """
      
        """
        return ((self.red**2 + self.green**2 + self.blue**2)/3)**(1/2)
    
    
    def BIM(self):
        """
        
        """
        return ((self.red*2 + self.green*2 + self.blue*2)/3)**(1/2)
    
    def HI(self):
        """
      
        """
        return (self.red*2 - self.green - self.blue)/(self.green - self.blue)
    
    def SI(self):
        """
       
        """
        return (self.red - self.blue)/(self.red - self.blue)
    
    def VARI(self):
        """
        
        """
        return (self.green - self.red)/(self.green + self.red - self.blue) 
    
    def BdivG(self):
        """
        
        """
        return (self.blue)/(self.green)  
    
    def BCC(self):
        """
        
        """
        return (self.blue/(self.red + self.green + self.blue) )
    
    def CIVE(self):
        """
        
        """
        return ((0.441*self.red)-(0.811*self.green)+(0.385*self.blue)+18.78745)
    
    def COM1(self):
        """
     
        """
        return (((2*self.green)-self.red-self.blue)+((0.441*self.red)-(0.811*self.green)+(0.385*self.blue)+18.78745)+((3*self.green)-(2.4*self.red)-self.blue)+(self.green/(self.red**0.667*self.blue**0.334)))
    
    def COM2(self):
        """
     
        """
        return ((0.36*((2*self.green)-self.red-self.blue))+(0.47*((0.441*self.red)-(0.811*self.green)+(0.385*self.blue)+18.78745))+(0.17*(self.green/(self.red**0.667*self.blue**0.334))))
    
    def ExG(self):
        """
      
        """
        return ((2*self.green)-self.red-self.blue)
    
    def ExG2(self):
        """
     
        """
        return (((2*self.green)-self.red-self.blue)/(self.red+self.green+self.blue))
    
    def ExGR(self):
        """
       
        """
        return ((3*self.green)-(2.4*self.red)-self.blue)
    
    def EXR(self):
        """
     
        """
        return ((1.4*self.red)-self.green)
    
    def GmnB(self):
        """
      
        """
        return (self.green-self.blue)
    
    def GmnR(self):
        """
        
        """
        return (self.green-self.red)
    
    def GdivB(self):
        """
      
        """
        return (self.green/self.blue)
    
    def GdivR(self):
        """
        
        """
        return (self.green/self.red)
    
    def GCC(self):
        """
        
        """
        return (self.green/(self.red+self.green+self.blue))
    
    def MExG(self):
        """
       
        """
        return ((1.262*self.green)-(0.884*self.red)-(0.311*self.blue))
    
    def MGVRI(self):
        """
       
        """
        return ((self.green**2)-(self.red**2))/((self.green**2)+(self.red**2))
    
    def NDI(self):
        """
       
        """
        return 128*(((self.green-self.red)/(self.green+self.red))+1)
    
    def NDRBI(self):
        """
       
        """
        return ( self.red - self.blue)/( self.red+ self.blue)
    
    def NGBDI(self):
        """
       
        """
        return ( self.red - self.blue)/( self.green + self.blue)

    def RmnB(self):
        """
        
        """
        return ( self.red - self.blue)
    
    def RdiB(self):
        """
        
        """
        return ( self.red / self.blue)
    
    
    def RCC(self):
        """
        
        """
        return ( self.red/( self.red + self.green + self.blue))
    
    def MRCCbyAlper(self):
        """
        
        """
        return ( self.red**3/( self.red + self.green + self.blue))
    
    def RGBVI(self):
        """
        
        """
        return ((( self.green**2)-( self.red * self.blue))/(( self.green**2)+( self.red * self.blue)))
    
    def TGI(self):
        """
        
        """
        return ( self.green-((0.39 * self.red)-(0.69 * self.blue)))
    
    def VEG(self):
        """
        
        """
        return self.green/( self.red**0.667* self.blue**0.334)
    
    def MyIndexi(self):
        """
        
        """
        return (( self.red - self.blue) / self.green)
    
    def MSRGR(self):
        """
        
        """
        return ( self.green / self.red)**(1/2)
    
    def TNDGR(self):
        """
        
        """
        return ( ( self.green - self.red)/( self.green + self.red)+ 0.5)**(1/2)
 
    
def convert_bytes(num):
     """
     this function will convert bytes to MB.... GB... etc
     """
     for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
         if num < 1024.0:
             return "%3.1f %s" % (num, x)
         num /= 1024.0


def file_size(file_path):
     """
     this function will return the file size
     """
     if os.path.isfile(file_path):
         file_info = os.stat(file_path)
         return convert_bytes(file_info.st_size)

## Setting the folders directories

### Plots shapefiles
plots = gpd.read_file('C:\\Users\\leoag\\Michigan State University\\MSU Dry Bean Breeding Lab - General\\UAS_Beans\\2021\\MRC\\Data_VIs\\Shapefile\\Shapefile_MRC_02.shp')

### Open attribute table of data
plots_names=list(plots.columns.values)
print(plots_names)

### Open with `fiona` to better explore the plots features
shp_file=fiona.open('C:\\Users\\leoag\\Michigan State University\\MSU Dry Bean Breeding Lab - General\\UAS_Beans\\2021\\MRC\\Data_VIs\\Shapefile\\Shapefile_MRC_02.shp', "r") 
#first feature of the shapefile
shp_file_inf = shp_file.next()
print(shp_file_inf)

### Image diretory (.tif aka orthomosaic)   
img_dir = 'C:\\Users\\leoag\\Michigan State University\\MSU Dry Bean Breeding Lab - General\\UAS_Beans\\2021\\MRC\\Ortho'

#### Getting the images from the directory according to the image format 
img_list = os.listdir(img_dir)
img_list = [v for v in img_list if v.endswith('.tif')]
print(img_list)

### Creating a list with the image names to save in the final table
# Create an empty list to store the img names
img_list_names = []

for l in range(len(img_list)):
    #print(l)
    img_list_names.append(str(img_list[l]))

print(img_list_names)

### Checking the file size

file_path = os.path.join(img_dir, img_list[1])
print (file_size(file_path))

if os.stat(file_path).st_size > 1000000000 :
   print("This python code is not recommended to large files!!! Please, use the alternative python script avalaible at:[gitpage](website)")
else:
   print("This script will performing well. Continue")


## Runing the loop
### j img numbers
### v VIs
### p plots
### 'count' 'min' 'mean' 'max' 'median' 'std' function of pixel extraction
for j in range(len(img_list)):
    print(j)
    
    # Change the current working directory
    os.chdir('C:\\Users\\leoag\\Michigan State University\\MSU Dry Bean Breeding Lab - General\\UAS_Beans\\2021\\MRC\\Ortho')
   
    ## Reading the img and save the raster foe each band
    ortho = rasterio.open(img_list[j])
    print(ortho)
    crs = ortho.crs
    #dem_arr = dem.read()
    ortho_arr = ortho.read()
    #print(ortho_arr)
    #extract color bands
    # Add here more bands in case for MS
    red=ortho_arr[0]
    green=ortho_arr[1]
    blue=ortho_arr[2]
    #redEdge[]
    #nir[]
   
    # Obtaining the VI HUE to use in the soil extraction
    cl = IndexCalculation()            
    cl.setMatrices(red=red, green=green, blue=blue)
    
    HUE_msk=cl.HUE()
    #plt.imshow(HUE_msk)
    
    # Removing the R band using a mask obrained from the HUE and replacing to NaN values (no vlaues)
    msk_soil_R = ma.masked_where(HUE_msk > 0.7, ortho_arr[0])
    msk_soil_R = ma.filled(msk_soil_R.astype(float), np.nan)
    #plt.imshow(msk_soil_R)
    
    # Removing the G band using a mask obrained from the HUE and replacing to NaN values (no vlaues)
    msk_soil_G = ma.masked_where(HUE_msk > 0.7, ortho_arr[1])
    msk_soil_G = ma.filled(msk_soil_G.astype(float), np.nan)
    #plt.imshow(msk_soil_G)
    
    # Removing the B band using a mask obrained from the HUE and replacing to NaN values (no vlaues)
    msk_soil_B = ma.masked_where(HUE_msk > 0.7, ortho_arr[2])
    msk_soil_B = ma.filled(msk_soil_B.astype(float), np.nan)
    #plt.imshow(msk_soil_B)

    # Creating the new matrice with the adjusted bands (witout soil)
    cl.setMatrices(red=msk_soil_R, green=msk_soil_G, blue=msk_soil_B)
      
    # Getting the VIs (please, see the list with all the VIs available)
    R=msk_soil_R
    G=msk_soil_G
    B=msk_soil_B
    GLI=cl.GLI()
    CI=cl.CI()
    HUE=cl.HUE()
    I=cl.I()
    NGRDI=cl.NGRDI()
    RI=cl.RI()
    IF=cl.IF()
    BI=cl.BI()
    BIM=cl.BIM()
    HI=cl.HI()
    SI=cl.SI()
    VARI=cl.VARI()
    BdivG=cl.BdivG()
    BCC=cl.BCC()
    CIVE=cl.CIVE()
    COM1=cl.COM1()
    COM2=cl.COM2()
    ExG=cl.ExG()
    ExG2=cl.ExG2()
    ExGR=cl.ExGR()
    EXR=cl.EXR()
    GmnB=cl.GmnB()
    GmnR=cl.GmnR()
    GdivB=cl.GdivB()
    GdivR=cl.GdivR()
    GCC=cl.GCC()
    MExG=cl.MExG()
    MGVRI=cl.MGVRI()
    NDI=cl.NDI()
    NDRBI=cl.NDRBI()
    NGBDI=cl.NGBDI()
    RmnB=cl.RmnB()
    RdiB=cl.RdiB()
    RCC=cl.RCC()
    MRCCbyAlper=cl.MRCCbyAlper()
    RGBVI=cl.RGBVI()
    TGI=cl.TGI()
    VEG=cl.VEG()
    MyIndexi=cl.MyIndexi()
    MSRGR=cl.MSRGR()
    TNDGR=cl.TNDGR()

    # Getting a list to enter in the next lopp
    list_VIs = [R, G, B, GLI ,  CI ,  HUE ,  I ,  NGRDI ,  RI ,  IF , 
                BI , BIM , HI , NGRDI , 
                SI , VARI , BdivG , BCC , CIVE , COM1 , COM2 , ExG , ExG2 , ExGR , EXR , 
                GmnB , GmnR , GdivB , GdivR , GCC , MExG , MGVRI , NDI , NDRBI , NGBDI , 
                RmnB , RdiB , RCC , MRCCbyAlper ,RGBVI , TGI , VEG , MyIndexi , MSRGR , TNDGR ]  
    
    # Defining the name lista (same than the 'list_VIs')
    list_VIs_names = ['R', 'G', 'B',"GLI", "CI", "HUE", "I", "NGRDI", "RI", "IF",
                      'BI','BIM','HI','NGRDI','SI','VARI','BdivG',
                      'BCC','CIVE','COM1','COM2','ExG','ExG2','ExGR','EXR',
                      'GmnB','GmnR','GdivB','GdivR','GCC','MExG','MGVRI',
                      'NDI','NDRBI','NGBDI','RmnB','RdiB','RCC','MRCCbyAlper',
                      'RGBVI','TGI','VEG','MyIndexi','MSRGR','TNDGR']  
         
    
    # loop over features
    for  i in range(len(list_VIs)):
        print(i)
       # compute zonal statistics on each plot
        plot_zs = rs.zonal_stats(plots, 
                list_VIs[i], 
                nodata=0, 
                affine=ortho.transform, 
                geojson_out=True, 
                copy_properties=True, 
                stats="count min mean max median std")
    
            # build dataframe and display first 5 records
        results_df = gpd.GeoDataFrame.from_features(plot_zs)
        #display(results_df.head())
        results_df['VIs'] = list_VIs_names[i]
        
        if i == 0 :
            results_df2=results_df
            
        else :
            results_df2=pd.concat([results_df2,results_df], ignore_index=True)                     
            
    results_df2['Fligtg_date'] = img_list_names[j]
         
    if j == 0 :
         results_df3=results_df2
         
    else :
         results_df3=pd.concat([results_df3,results_df2], ignore_index=True)
        
    
results_df3.dtypes
results_df3.to_csv('VIs_RGB_MRC.csv')
    

    
    
    
    




