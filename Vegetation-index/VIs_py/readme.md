# *Python script to perform Vegetation indices (VIs) for phenotyping*
# Class Summary
> This algorithm consists in calculating vegetation indices looping through images files and plot numbers. These
indices can be used for precision agriculture for example (or remote sensing). There are functions to define the data and to calculate the
implemented indices.

# [Vegetation index](https://en.wikipedia.org/wiki/Vegetation_Index)

A Vegetation Index (VI) is a spectral transformation of two or more bands
designed to enhance the contribution of vegetation properties and allow
reliable spatial and temporal inter-comparisons of terrestrial
photosynthetic activity and canopy structural variations

# Information about channels (Wavelength range for each)

* [nir - near-infrared](https://www.malvernpanalytical.com/br/products/technology/near-infrared-spectroscopy)
    Wavelength Range 700 nm to 2500 nm
    
* [Red Edge](https://en.wikipedia.org/wiki/Red_edge)
    Wavelength Range 680 nm to 730 nm
    
* [red](https://en.wikipedia.org/wiki/Color)
    Wavelength Range 635 nm to 700 nm
    
* [blue](https://en.wikipedia.org/wiki/Color)
    Wavelength Range 450 nm to 490 nm
    
* [green](https://en.wikipedia.org/wiki/Color)
    Wavelength Range 520 nm to 560 nm
    
# Implemented index list 
> 73 in total from RGB and Multispectral (MS) sensors </p>

`abbreviationOfIndexName` -- list of channels used

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
            #"PSRI"             --  red, green, redEdge
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

## List of all index implemented 
        ARVI2","ATSAVI","BCC","BdivG","BI","BIM","BNDVI","BWDRVI",
        "CCCI","CI","CIgreen","CIrededge","CIVE","COM1","COM2","CTVI","CVI","DVI",
        "EVI","ExG","ExG2","ExGR","EXR","GBNDVI","GCC","GdivB","GdivR","GDVI","GEMI",
        "GLI","GmnB","GmnR","GNDVI","GOSAVI","GRNDVI","GSAVI","HI","HUE","I","IF","IPVI",
        "MExG","MGVRI","MRCCbyAlper","MRVI","MSAVI","MSRGR","MyIndexi","NDI","NDRBI",
        "NDRE","NDVI","NGBDI","NGRDI","NormG","NormNIR","NormR","PNDVI","RBNDVI","RCC","RdiB",
        "redEdgeNDVI","RGBVI","RI","RmnB","RVI","SI","TGI","TNDGR","TVI","VARI","VEG" , "PSRI"

## List of index for multispectral (MS) camera - 4 or 5 channels
 "ARVI2", "CCCI", "CVI", "NDVI", "redEdgeNDVI", "GNDVI",
 "GRNDVI", "ATSAVI", "CIgreen", "CIrededge", "CTVI", "GDVI",
 "GEMI", "GOSAVI", "GSAVI", "IPVI", "RVI", "MRVI",
 "MSAVI", "NormG", "NormNIR", "NormR", "NGRDI", "RI", "DVI",
 "TVI", "NDRE", "PSRI
 
## List of index just with RGB channels
 "GLI", "CI", "HUE", "I", "NGRDI", "RI", "IF", 'BI','BIM',
 'HI','NGRDI','SI','VARI','BdivG','BCC','CIVE','COM1','COM2',
 'ExG','ExG2','ExGR','EXR','GmnB','GmnR','GdivB','GdivR',
 'GCC','MExG','MGVRI','NDI','NDRBI','NGBDI','RmnB',
 'RdiB','RCC','MRCCbyAlper','RGBVI','TGI','VEG','MyIndexi','MSRGR','TNDGR'

# References
* Vegetation indices contributions from: 
[João Gustavo A. Amorim](https://github.com/TheAlgorithms/Python/blob/master/digital_image_processing/index_calculation.py)

____________________________________________________________________________________
____________________________________________________________________________________

# Import the packages to the python environment


```python
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
from pathlib import Path
warnings.filterwarnings('ignore') #don't display warnings
```

## Check RAM available


```python
from psutil import virtual_memory
ram_gb = virtual_memory().total / 1e9
print('Your runtime has {:.2f} gigabytes of available RAM\n'.format(ram_gb))

if ram_gb < 20:
  print('You should consider increasing the memory RAM before to start')
else:
  print('You are using a high-RAM runtime!')
```

    Your runtime has 34.16 gigabytes of available RAM
    
    You are using a high-RAM runtime!
    

## Function implemented to perform the index


```python
class IndexCalculation:
  
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
            "PSRI":self.PSRI,
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
        
    def PSRI(self):
        """
        Plant senescence reflectance index 
        """
        return (self.red - self.green) / (self.redEdge)


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
```

```python
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
```

# **RGB bands**

# Alternative #1:

## Setting the folders directories

### Plots shapefiles


```python
plots = gpd.read_file('Export_Output.shp')
```

### Open attribute table of data


```python
plots_names=list(plots.columns.values)
print(plots_names)
```
    

### Open with `fiona` to better explore the plots features


```python
shp_file=fiona.open('Export_Output.shp', "r") 
#first feature of the shapefile
shp_file_inf = shp_file.next()
print(shp_file_inf)
```
 

### Image diretory (.tif aka orthomosaic)  


```python
img_dir = Path("~\\path")
```

#### Getting the images from the directory according to the image format 


```python
img_list = os.listdir(img_dir)
img_list = [v for v in img_list if v.endswith('.tif')]
print(img_list)
```
    

### Creating a list with the image names to save in the final table


```python
# Create an empty list to store the img names
img_list_names = []

for l in range(len(img_list)):
    print(l)
    img_list_names.append(str(img_list[l]))

print(img_list_names)
```

# Checking the file size
```python
file_path = os.path.join(img_dir, img_list[1])
print (file_size(file_path))

if os.stat(file_path).st_size > 1000000000:
   print("This python code is not recommended for large files!!! Please, use the alternative python script available at "Alternative #2 below")
else:
   print("This script will performing well. Continue")

```

## Runing the loop
> **j**  img numbers </p>
> **i**  VIs</p>
* From the `rasterstats` library
> **p**  plots</p>
> **f**  function of pixel extraction: 'count' 'min' 'mean' 'max' 'median' 'std' 


```python
for j in range(len(img_list)):
    print(j)
    
    # Change the current working directory
    os.chdir(img_dir) #definy here the img directory
   
    
    ## Reading the img and save the raster foe each band
    ortho = rasterio.open(img_list[j])
    print(ortho)
    #crs = ortho.crs
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
    list_VIs_names = ['R', 'G', 'B', "GLI", "CI", "HUE", "I", "NGRDI", "RI", "IF",
                      'BI','BIM','HI','NGRDI','SI','VARI','BdivG',
                      'BCC','CIVE','COM1','COM2','ExG','ExG2','ExGR','EXR',
                      'GmnB','GmnR','GdivB','GdivR','GCC','MExG','MGVRI',
                      'NDI','NDRBI','NGBDI','RmnB','RdiB','RCC','MRCCbyAlper',
                      'RGBVI','TGI','VEG','MyIndexi','MSRGR','TNDGR']  
         
    
    # loop over features using `rasterstats` library
    for  i in range(len(list_VIs)):
        print(i)
       # compute zonal statistics on each plot and using 6 function of extraction
        plot_zs = rs.zonal_stats(plots, 
                list_VIs[i], 
                nodata=0, 
                affine=ortho.transform, 
                geojson_out=True, 
                copy_properties=True, 
                stats="count min mean max median std")
    
            # build dataframe and saving the VI name
        results_df = gpd.GeoDataFrame.from_features(plot_zs)
        results_df['VIs'] = list_VIs_names[i]
        
        if i == 0 :
            results_df2=results_df
            
        else :
            results_df2=pd.concat([results_df2,results_df], ignore_index=True)                     
    # Saving the image name and storing the results in a final file       
    results_df2['Fligtg_date'] = img_list_names[j]
         
    if j == 0 :
         results_df3=results_df2
         
    else :
         results_df3=pd.concat([results_df3,results_df2], ignore_index=True)
        
```


## Saving the results in a csv file


```python
results_df3.dtypes
results_df3.to_csv('VIs_names.csv')
```



# Alternative #2:
> Thanks to: [Jose Cleydson Ferreira Silva](https://www.linkedin.com/in/jose-cleydson-ferreira-silva-a790b617b/)
# Import the packages to the python environment

```python
import numpy as np
import numpy.ma as ma
#from matplotlib import pyplot as plt
import rasterio
import geopandas as gpd
#import rasterstats as rs
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
from rasterio.mask import mask
```
## Setting the folders directories

### Plots shapefiles


```python
shp_file = "Export_Output.shp"
```

### Open attribute table of data


```python
plots_names=list(plots.columns.values)
print(plots_names)
```
    

### Open with `fiona` to better explore the plots features


```python
shp_file_info=fiona.open('Export_Output.shp', "r") 
#first feature of the shapefile
shp_file_inf = shp_file_info.next()
print(shp_file_inf)
```
 

### Image diretory (.tif aka orthomosaic)  


```python
img_dir = Path("~\\path")
```

#### Getting the images from the directory according to the image format 


```python
img_list = os.listdir(img_dir)
img_list = [v for v in img_list if v.endswith('.tif')]
print(img_list)
```
    

### Creating a list with the image names to save in the final table


```python
# Create an empty list to store the img names
img_list_names = []

for l in range(len(img_list)):
    print(l)
    img_list_names.append(str(img_list[l]))

print(img_list_names)
```

## Checking the file size
```python
file_path = os.path.join(img_dir, img_list[1])
print (file_size(file_path))

if os.stat(file_path).st_size > 0:
   print("This python code is not recommended to large files!!! Please, use the alternative python script avalaible at:[gitpage](website)")
else:
   print("This script will performing well. Continue")
```

## Get the plot ID name
```python
plot_num_field_name = "MNPlot"
```
## Creating a array to saving the results
```python
array = list()
```

## Creating a name list to save the outputs
```python
colname = ['Flight_date','Plot_ID',
           "R_mean","R_median","R_count","R_std",
           "G_mean","G_median","G_count","G_std",
           "B_mean","B_median","B_count","B_std","BCC_mean","BCC_median","BCC_count",
           
           "BCC_std","BdivG_mean","BdivG_median","BdivG_count","BdivG_std","BI_mean",
           "BI_median","BI_count","BI_std","BIM_mean","BIM_median","BIM_count","BIM_std",
           "CI_mean","CI_median","CI_count","CI_std","CIVE_mean","CIVE_median","CIVE_count",
           "CIVE_std","COM1_mean","COM1_median","COM1_count","COM1_std","COM2_mean","COM2_median",
           "COM2_count","COM2_std","ExG_mean","ExG_median","ExG_count","ExG_std","ExG2_mean",
           "ExG2_median","ExG2_count","ExG2_std","ExGR_mean","ExGR_median","ExGR_count",
           "ExGR_std","EXR_mean","EXR_median","EXR_count","EXR_std","GCC_mean","GCC_median","GCC_count","GCC_std","GdivB_mean","GdivB_median",
           "GdivB_count","GdivB_std","GdivR_mean","GdivR_median","GdivR_count","GdivR_std",
           "GLI_mean","GLI_median","GLI_count","GLI_std","GmnB_mean","GmnB_median","GmnB_count",
           "GmnB_std","GmnR_mean","GmnR_median","GmnR_count","GmnR_std","HI_mean","HI_median",
           "HI_count","HI_std","HUE_mean","HUE_median","HUE_count","HUE_std","I_mean","I_median",
           "I_count","I_std","IF_mean","IF_median","IF_count","IF_std","MExG_mean","MExG_median",
           "MExG_count","MExG_std","MGVRI_mean","MGVRI_median","MGVRI_count","MGVRI_std","MRCCbyAlper_mean",
           "MRCCbyAlper_median","MRCCbyAlper_count","MRCCbyAlper_std","MSRGR_mean","MSRGR_median","MSRGR_count",
           "MSRGR_std","MyIndexi_mean","MyIndexi_median","MyIndexi_count","MyIndexi_std","NDI_mean","NDI_median",
           "NDI_count","NDI_std","NDRBI_mean","NDRBI_median","NDRBI_count","NDRBI_std","NGBDI_mean","NGBDI_median",
           "NGBDI_count","NGBDI_std","NGRDI_mean","NGRDI_median","NGRDI_count","NGRDI_std","RCC_mean","RCC_median","RCC_count","RCC_std","RdiB_mean","RdiB_median","RdiB_count","RdiB_std","RGBVI_mean",
           "RGBVI_median","RGBVI_count","RGBVI_std","RI_mean","RI_median","RI_count","RI_std","RmnB_mean","RmnB_median",
           "RmnB_count","RmnB_std","SI_mean","SI_median","SI_count","SI_std","TGI_mean","TGI_median","TGI_count","TGI_std",
           "TNDGR_mean","TNDGR_median","TNDGR_count","TNDGR_std","VARI_mean","VARI_median","VARI_count","VARI_std",
           "VEG_mean","VEG_median","VEG_count","VEG_std"]

```

## Runing the loop
> **j**  img numbers </p>
> **p**  plots</p>
> **i**  VIs</p>
* From the `numpy` library
> **f**  function of pixel extraction: 'mean' 'media' 'count' 'std

```python
for tiffile in img_list_names:
    print(tiffile)
    with fiona.open(shp_file, "r") as shapefile:
    
        for feature in shapefile:
    
            shape = [feature["geometry"]]
            
            Plot_ID = feature["properties"][plot_num_field_name]
            print(Plot_ID,"-",tiffile)
            
            # Change the current working directory
            os.chdir(img_dir)
            
            #open and mask the multispectral image to the plot
            
            with rasterio.open(tiffile, "r") as ras:
                out_image, out_transform = mask(ras, shape, crop=True, nodata=0)
                
                out_image = ma.masked_where(out_image == 0, out_image)
                out_image = ma.filled(out_image.astype(float), np.nan)

                #ep.plot_rgb(out_image, rgb=[0, 1, 2], title="Red Green Blue", stretch=True)
                      
                #extract color bands
                red = out_image[0,:,:].astype('float32')
                green = out_image[1,:,:].astype('float32')
                blue = out_image[2,:,:].astype('float32')
                
    
                #RGBIndex = ["GLI", "CI", "HUE", "I", "NGRDI", "RI", "S", "IF"]           
                cl = IndexCalculation()            
                cl.setMatrices(red=red, green=green, blue=blue)
                
                HUE_msk=cl.HUE()
                #plt.imshow(HUE_msk)
                
                # Removing the R band using a mask obrained from the HUE and replacing to NaN values (no vlaues)
                msk_soil_R = ma.masked_where(HUE_msk > 0.7, red)
                msk_soil_R = ma.filled(msk_soil_R.astype(float), np.nan)
                #plt.imshow(msk_soil_R)
                
                # Removing the G band using a mask obrained from the HUE and replacing to NaN values (no vlaues)
                msk_soil_G = ma.masked_where(HUE_msk > 0.7, green)
                msk_soil_G = ma.filled(msk_soil_G.astype(float), np.nan)
                #plt.imshow(msk_soil_G)
                
                # Removing the B band using a mask obrained from the HUE and replacing to NaN values (no vlaues)
                msk_soil_B = ma.masked_where(HUE_msk > 0.7, blue)
                msk_soil_B = ma.filled(msk_soil_B.astype(float), np.nan)
                #plt.imshow(msk_soil_B)

                # Creating the new matrice with the adjusted bands (witout soil)
                cl.setMatrices(red=msk_soil_R, green=msk_soil_G, blue=msk_soil_B)
               
                # print values
                array.append([tiffile,Plot_ID, 
                         np.nanmean(msk_soil_R.flatten()),np.nanmedian(msk_soil_R.flatten()),np.count_nonzero(~np.isnan(msk_soil_R.flatten())),np.nanstd(msk_soil_R.flatten()),
                         np.nanmean(msk_soil_G.flatten()),np.nanmedian(msk_soil_G.flatten()),np.count_nonzero(~np.isnan(msk_soil_G.flatten())),np.nanstd(msk_soil_G.flatten()),
                         np.nanmean(msk_soil_B.flatten()),np.nanmedian(msk_soil_B.flatten()),np.count_nonzero(~np.isnan(msk_soil_B.flatten())),np.nanstd(msk_soil_B.flatten()),
                         
                         np.nanmean(cl.BCC().flatten()),np.nanmedian(cl.BCC().flatten()),np.count_nonzero(~np.isnan(cl.BCC().flatten())),np.nanstd(cl.BCC().flatten()),
                         np.nanmean(cl.BdivG().flatten()),np.nanmedian(cl.BdivG().flatten()),np.count_nonzero(~np.isnan(cl.BdivG().flatten())),np.nanstd(cl.BdivG().flatten()),
                         np.nanmean(cl.BI().flatten()),np.nanmedian(cl.BI().flatten()),np.count_nonzero(~np.isnan(cl.BI().flatten())),np.nanstd(cl.BI().flatten()),
                         np.nanmean(cl.BIM().flatten()),np.nanmedian(cl.BIM().flatten()),np.count_nonzero(~np.isnan(cl.BIM().flatten())),np.nanstd(cl.BIM().flatten()),
                         np.nanmean(cl.CI().flatten()),np.nanmedian(cl.CI().flatten()),np.count_nonzero(~np.isnan(cl.CI().flatten())),np.nanstd(cl.CI().flatten()),
                         np.nanmean(cl.CIVE().flatten()),np.nanmedian(cl.CIVE().flatten()),np.count_nonzero(~np.isnan(cl.CIVE().flatten())),np.nanstd(cl.CIVE().flatten()),
                         np.nanmean(cl.COM1().flatten()),np.nanmedian(cl.COM1().flatten()),np.count_nonzero(~np.isnan(cl.COM1().flatten())),np.nanstd(cl.COM1().flatten()),
                         np.nanmean(cl.COM2().flatten()),np.nanmedian(cl.COM2().flatten()),np.count_nonzero(~np.isnan(cl.COM2().flatten())),np.nanstd(cl.COM2().flatten()), 
                         np.nanmean(cl.ExG().flatten()),np.nanmedian(cl.ExG().flatten()),np.count_nonzero(~np.isnan(cl.ExG().flatten())),np.nanstd(cl.ExG().flatten()),
                         np.nanmean(cl.ExG2().flatten()),np.nanmedian(cl.ExG2().flatten()),np.count_nonzero(~np.isnan(cl.ExG2().flatten())),np.nanstd(cl.ExG2().flatten()),
                         np.nanmean(cl.ExGR().flatten()),np.nanmedian(cl.ExGR().flatten()),np.count_nonzero(~np.isnan(cl.ExGR().flatten())),np.nanstd(cl.ExGR().flatten()),
                         np.nanmean(cl.EXR().flatten()),np.nanmedian(cl.EXR().flatten()),np.count_nonzero(~np.isnan(cl.EXR().flatten())),np.nanstd(cl.EXR().flatten()),
                         np.nanmean(cl.GCC().flatten()),np.nanmedian(cl.GCC().flatten()),np.count_nonzero(~np.isnan(cl.GCC().flatten())),np.nanstd(cl.GCC().flatten()),
                         np.nanmean(cl.GdivB().flatten()),np.nanmedian(cl.GdivB().flatten()),np.count_nonzero(~np.isnan(cl.GdivB().flatten())),np.nanstd(cl.GdivB().flatten()),
                         np.nanmean(cl.GdivR().flatten()),np.nanmedian(cl.GdivR().flatten()),np.count_nonzero(~np.isnan(cl.GdivR().flatten())),np.nanstd(cl.GdivR().flatten()),
                         np.nanmean(cl.GLI().flatten()),np.nanmedian(cl.GLI().flatten()),np.count_nonzero(~np.isnan(cl.GLI().flatten())),np.nanstd(cl.GLI().flatten()),
                         np.nanmean(cl.GmnB().flatten()),np.nanmedian(cl.GmnB().flatten()),np.count_nonzero(~np.isnan(cl.GmnB().flatten())),np.nanstd(cl.GmnB().flatten()),
                         np.nanmean(cl.GmnR().flatten()),np.nanmedian(cl.GmnR().flatten()),np.count_nonzero(~np.isnan(cl.GmnR().flatten())),np.nanstd(cl.GmnR().flatten()),
                         np.nanmean(cl.HI().flatten()),np.nanmedian(cl.HI().flatten()),np.count_nonzero(~np.isnan(cl.HI().flatten())),np.nanstd(cl.HI().flatten()),
                         np.nanmean(cl.HUE().flatten()),np.nanmedian(cl.HUE().flatten()),np.count_nonzero(~np.isnan(cl.HUE().flatten())),np.nanstd(cl.HUE().flatten()),
                         np.nanmean(cl.I().flatten()),np.nanmedian(cl.I().flatten()),np.count_nonzero(~np.isnan(cl.I().flatten())),np.nanstd(cl.I().flatten()),
                         np.nanmean(cl.IF().flatten()),np.nanmedian(cl.IF().flatten()),np.count_nonzero(~np.isnan(cl.IF().flatten())),np.nanstd(cl.IF().flatten()),
                         np.nanmean(cl.MExG().flatten()),np.nanmedian(cl.MExG().flatten()),np.count_nonzero(~np.isnan(cl.MExG().flatten())),np.nanstd(cl.MExG().flatten()),
                         np.nanmean(cl.MGVRI().flatten()),np.nanmedian(cl.MGVRI().flatten()),np.count_nonzero(~np.isnan(cl.MGVRI().flatten())),np.nanstd(cl.MGVRI().flatten()),
                         np.nanmean(cl.MRCCbyAlper().flatten()),np.nanmedian(cl.MRCCbyAlper().flatten()),np.count_nonzero(~np.isnan(cl.MRCCbyAlper().flatten())),np.nanstd(cl.MRCCbyAlper().flatten()),
                         np.nanmean(cl.MSRGR().flatten()),np.nanmedian(cl.MSRGR().flatten()),np.count_nonzero(~np.isnan(cl.MSRGR().flatten())),np.nanstd(cl.MSRGR().flatten()),
                         np.nanmean(cl.MyIndexi().flatten()),np.nanmedian(cl.MyIndexi().flatten()),np.count_nonzero(~np.isnan(cl.MyIndexi().flatten())),np.nanstd(cl.MyIndexi().flatten()),
                         np.nanmean(cl.NDI().flatten()),np.nanmedian(cl.NDI().flatten()),np.count_nonzero(~np.isnan(cl.NDI().flatten())),np.nanstd(cl.NDI().flatten()),
                         np.nanmean(cl.NDRBI().flatten()),np.nanmedian(cl.NDRBI().flatten()),np.count_nonzero(~np.isnan(cl.NDRBI().flatten())),np.nanstd(cl.NDRBI().flatten()) ,
                         np.nanmean(cl.NGBDI().flatten()),np.nanmedian(cl.NGBDI().flatten()),np.count_nonzero(~np.isnan(cl.NGBDI().flatten())),np.nanstd(cl.NGBDI().flatten()),
                         np.nanmean(cl.NGRDI().flatten()),np.nanmedian(cl.NGRDI().flatten()),np.count_nonzero(~np.isnan(cl.NGRDI().flatten())),np.nanstd(cl.NGRDI().flatten()),
                         np.nanmean(cl.RCC().flatten()),np.nanmedian(cl.RCC().flatten()),np.count_nonzero(~np.isnan(cl.RCC().flatten())),np.nanstd(cl.RCC().flatten()),
                         np.nanmean(cl.RdiB().flatten()),np.nanmedian(cl.RdiB().flatten()),np.count_nonzero(~np.isnan(cl.RdiB().flatten())),np.nanstd(cl.RdiB().flatten()),
                         np.nanmean(cl.RGBVI().flatten()),np.nanmedian(cl.RGBVI().flatten()),np.count_nonzero(~np.isnan(cl.RGBVI().flatten())),np.nanstd(cl.RGBVI().flatten()),
                         np.nanmean(cl.RI().flatten()),np.nanmedian(cl.RI().flatten()),np.count_nonzero(~np.isnan(cl.RI().flatten())),np.nanstd(cl.RI().flatten()),
                         np.nanmean(cl.RmnB().flatten()),np.nanmedian(cl.RmnB().flatten()),np.count_nonzero(~np.isnan(cl.RmnB().flatten())),np.nanstd(cl.RmnB().flatten()),
                         np.nanmean(cl.SI().flatten()),np.nanmedian(cl.SI().flatten()),np.count_nonzero(~np.isnan(cl.SI().flatten())),np.nanstd(cl.SI().flatten()),
                         np.nanmean(cl.TGI().flatten()),np.nanmedian(cl.TGI().flatten()),np.count_nonzero(~np.isnan(cl.TGI().flatten())),np.nanstd(cl.TGI().flatten()) ,
                         np.nanmean(cl.TNDGR().flatten()),np.nanmedian(cl.TNDGR().flatten()),np.count_nonzero(~np.isnan(cl.TNDGR().flatten())),np.nanstd(cl.TNDGR().flatten()),
                         np.nanmean(cl.VARI().flatten()),np.nanmedian(cl.VARI().flatten()),np.count_nonzero(~np.isnan(cl.VARI().flatten())),np.nanstd(cl.VARI().flatten()),
                         np.nanmean(cl.VEG().flatten()),np.nanmedian(cl.VEG().flatten()),np.count_nonzero(~np.isnan(cl.VEG().flatten())),np.nanstd(cl.VEG().flatten())
                        
                         ])
                                
```

## Saving the output into a csv table

```python
my_array = np.array(array)
dataframe = pd.DataFrame(my_array, columns = colname)
dataframe.to_csv('Output_VIs.csv')
```



# **Multispectral bands**


## Plots shapefiles
```python
plots_shp = gpd.read_file("C:\\Users\\leoag\\Michigan State University\\MSU Dry Bean Breeding Lab - General\\UAS_Beans\\2021\\MRC\\Data_VIs\\Shapefile\\Shapefile_MRC_02.shp")

#Open attribute table of data
plots_names=list(plots_shp.columns.values)
print(plots_names)
```

```python
## Open with `fiona` to better explore the plots features
shp_file = "C:\\Users\\leoag\\Michigan State University\\MSU Dry Bean Breeding Lab - General\\UAS_Beans\\2021\\MRC\\Data_VIs\\Shapefile\\Shapefile_MRC_02.shp"
```
```python
## Image diretory (.tif aka orthomosaic)   
img_dir = Path("~\\path")
```
```python
### Getting the images from the directory according to the image format 
img_list = os.listdir(img_dir)
img_list = [v for v in img_list if v.endswith('.tif')]
print(img_list) 
```
```python
## Creating a list with the image names to save in the final table

#Create an empty list to store the img names
img_list_names = []

for l in range(len(img_list)):
    #print(l)
    img_list_names.append(str(img_list[l]))

print(img_list_names)
```
## Plot ID
```python
plot_num_field_name = "name"
```
## Array to store the outputs
```python
array = list()
```
## Multispectral indices
> 28 VIs in total from MS camera + 5 reflectance (bands)
 
```python
colname = ['Flight_date','Plot_ID',
           "R_mean","R_median","R_count","R_std",
           "G_mean","G_median","G_count","G_std",
           "B_mean","B_median","B_count","B_std",
           "RE_mean","RE_median","RE_count","RE_std",
           "NIR_mean","NIR_median","NIR_count","NIR_std",
           "NDVI_mean","NDVI_median","NDVI_count","NDVI_std",
           "GNDVI_mean","GNDVI_median","GNDVI_count","GNDVI_std",
           "ARVI2_mean","ARVI2_median","ARVI2_count","ARVI2_std", 
           "CCCI_mean","CCCI_median","CCCI_count","CCCI_std", 
           "CVI_mean","CVI_median","CVI_count","CVI_std", 
           "redEdgeNDVI_mean","redEdgeNDVI_median","redEdgeNDVI_count","redEdgeNDVI_std",
           "GRNDVI_mean","GRNDVI_median","GRNDVI_count","GRNDVI_std", 
           "ATSAVI_mean","ATSAVI_median","ATSAVI_count","ATSAVI_std", 
           "CIgreen_mean","CIgreen_median","CIgreen_count","CIgreen_std",
           "CIrededge_mean","CIrededge_median","CIrededge_count","CIrededge_std",
           "CTVI_mean","CTVI_median","CTVI_count","CTVI_std", 
           "GDVI_mean","GDVI_median","GDVI_count","GDVI_std",
           "GEMI_mean","GEMI_median","GEMI_count","GEMI_std",
           "GOSAVI_mean","GOSAVI_median","GOSAVI_count","GOSAVI_std",
           "GSAVI_mean","GSAVI_median","GSAVI_count","GSAVI_std",
           "IPVI_mean","IPVI_median","IPVI_count","IPVI_std", 
           "RVI_mean","RVI_median","RVI_count","RVI_std", 
           "MRVI_mean","MRVI_median","MRVI_count","MRVI_std", 
           "MSAVI_mean","MSAVI_median","MSAVI_count","MSAVI_std",
           "NormG_mean","NormG_median","NormG_count","NormG_std", 
           "NormNIR_mean","NormNIR_median","NormNIR_count","NormNIR_std",
           "NormR_mean","NormR_median","NormR_count","NormR_std", 
           "NGRDI_mean","NGRDI_median","NGRDI_count","NGRDI_std", 
           "RI_mean","RI_median","RI_count","RI_std", 
           "DVI_mean","DVI_median","DVI_count","DVI_std",
           "TVI_mean","TVI_median","TVI_count","TVI_std",
           "NDRE_mean","NDRE_median","NDRE_count","NDRE_std",
           "PSRI_mean","PSRI_median","PSRI_count","PSRI_std"
           
           ]
           
```

## Setting the image parameter to start the loop
```python
nbands=5 #k
startnbands=0 #j

n_flights = [1,2,3,4,5]
```
## Runing the loop
> **j**  img numbers </p>
> **p**  plots</p>
> **i**  VIs</p>
* From the `numpy` library
> **f**  function of pixel extraction: 'mean' 'media' 'count' 'std

```python
for tiffile in n_flights:

#for tiffile in img_list_names[startnbands:nbands]:
        
    if (tiffile == 1):
        startnbands = startnbands + 0 
    else:
        startnbands = startnbands + 5
        
    ms_flight = img_list_names[(startnbands):(tiffile*nbands)]
        
        #print(ms_flight)
   
    with fiona.open(shp_file, "r") as shapefile:
    
        for feature in shapefile:
    
            shape = [feature["geometry"]]
            
            Plot_ID = feature["properties"][plot_num_field_name]
            print(Plot_ID,"--",ms_flight[1])
            
            # Change the current working directory
            os.chdir(img_dir)
            #open and mask the multispectral image to the plot
            
            ## blue band
            with rasterio.open(ms_flight[0], "r") as ras:
                out_image_blue, out_transform_blue = mask(ras, shape, crop=True, nodata=0)
                
                out_image_blue = ma.masked_where(out_image_blue == 0, out_image_blue)
                out_image_blue = ma.filled(out_image_blue.astype(float), np.nan)
                
             ## green band
            with rasterio.open(ms_flight[1], "r") as ras:
                 out_image_green, out_transform_green = mask(ras, shape, crop=True, nodata=0)
                 
                 out_image_green = ma.masked_where(out_image_green == 0, out_image_green)
                 out_image_green = ma.filled(out_image_green.astype(float), np.nan)
                 
            ## nir band
            with rasterio.open(ms_flight[2], "r") as ras:
                out_image_nir, out_transform_nir = mask(ras, shape, crop=True, nodata=0)
                
                out_image_nir = ma.masked_where(out_image_nir == 0, out_image_nir)
                out_image_nir = ma.filled(out_image_nir.astype(float), np.nan)
                
             ## red edge band
            with rasterio.open(ms_flight[3], "r") as ras:
                 out_image_redEdge, out_transform_redEdge = mask(ras, shape, crop=True, nodata=0)
                 
                 out_image_redEdge = ma.masked_where(out_image_redEdge == 0, out_image_redEdge)
                 out_image_redEdge = ma.filled(out_image_redEdge.astype(float), np.nan)
                 
            ## red band
            with rasterio.open(ms_flight[4], "r") as ras:
                out_image_red, out_transform_red = mask(ras, shape, crop=True, nodata=0)
                
                out_image_red = ma.masked_where(out_image_red == 0, out_image_red)
                out_image_red = ma.filled(out_image_red.astype(float), np.nan)
                              
           
                cl = IndexCalculation() 
                #cl.setMatrices(red=red, green=green, blue=blue)
                
                cl.setMatrices(red=out_image_red, green=out_image_green, blue=out_image_blue,redEdge=out_image_redEdge, nir=out_image_nir)
                             
                NDVI_msk=cl.NDVI()
                
                # Removing the R band using a mask obrained from the HI and replacing to NaN values (no vlaues)
                msk_soil_R = ma.masked_where(NDVI_msk < 0.7, out_image_red)
                msk_soil_R = ma.filled(msk_soil_R.astype(float), np.nan)
                #plt.imshow(msk_soil_R)
                
                # Removing the G band using a mask obrained from the HI and replacing to NaN values (no vlaues)
                msk_soil_G = ma.masked_where(NDVI_msk < 0.7, out_image_green)
                msk_soil_G = ma.filled(msk_soil_G.astype(float), np.nan)
                #plt.imshow(msk_soil_G)
                
                # Removing the B band using a mask obrained from the HI and replacing to NaN values (no vlaues)
                msk_soil_B = ma.masked_where(NDVI_msk < 0.7, out_image_blue)
                msk_soil_B = ma.filled(msk_soil_B.astype(float), np.nan)
                #plt.imshow(msk_soil_B)

                # Removing the redEdge band using a mask obrained from the HI and replacing to NaN values (no vlaues)
                msk_soil_RE = ma.masked_where(NDVI_msk < 0.7, out_image_redEdge)
                msk_soil_RE = ma.filled(msk_soil_RE.astype(float), np.nan)
                
                # Removing the B band using a mask obrained from the HI and replacing to NaN values (no vlaues)
                msk_soil_Nir = ma.masked_where(NDVI_msk < 0.7, out_image_nir)
                msk_soil_Nir = ma.filled(msk_soil_Nir.astype(float), np.nan)
                
                
                # Creating the new matrice with the adjusted bands (witout soil)
                cl.setMatrices(red=msk_soil_R, green=msk_soil_G, blue=msk_soil_B, redEdge= msk_soil_RE, nir= msk_soil_Nir)
               
                # print values
                array.append([ms_flight[0][0:6],Plot_ID, 
                         np.nanmean(msk_soil_R.flatten()),np.nanmedian(msk_soil_R.flatten()),np.count_nonzero(~np.isnan(msk_soil_R.flatten())),np.nanstd(msk_soil_R.flatten()),
                         np.nanmean(msk_soil_G.flatten()),np.nanmedian(msk_soil_G.flatten()),np.count_nonzero(~np.isnan(msk_soil_G.flatten())),np.nanstd(msk_soil_G.flatten()),
                         np.nanmean(msk_soil_B.flatten()),np.nanmedian(msk_soil_B.flatten()),np.count_nonzero(~np.isnan(msk_soil_B.flatten())),np.nanstd(msk_soil_B.flatten()),
                         np.nanmean(msk_soil_RE.flatten()),np.nanmedian(msk_soil_RE.flatten()),np.count_nonzero(~np.isnan(msk_soil_RE.flatten())),np.nanstd(msk_soil_RE.flatten()),
                         np.nanmean(msk_soil_Nir.flatten()),np.nanmedian(msk_soil_Nir.flatten()),np.count_nonzero(~np.isnan(msk_soil_Nir.flatten())),np.nanstd(msk_soil_Nir.flatten()),
                         np.nanmean(cl.NDVI().flatten()),np.nanmedian(cl.NDVI().flatten()),np.count_nonzero(~np.isnan(cl.NDVI().flatten())),np.nanstd(cl.NDVI().flatten()),
                         np.nanmean(cl.GNDVI().flatten()),np.nanmedian(cl.GNDVI().flatten()),np.count_nonzero(~np.isnan(cl.GNDVI().flatten())),np.nanstd(cl.GNDVI().flatten()),
                         np.nanmean(cl.ARVI2().flatten()),np.nanmedian(cl.ARVI2().flatten()),np.count_nonzero(~np.isnan(cl.ARVI2().flatten())),np.nanstd(cl.ARVI2().flatten()),
                         np.nanmean(cl.CCCI().flatten()),np.nanmedian(cl.CCCI().flatten()),np.count_nonzero(~np.isnan(cl.CCCI().flatten())),np.nanstd(cl.CCCI().flatten()),
                         np.nanmean(cl.CVI().flatten()),np.nanmedian(cl.CVI().flatten()),np.count_nonzero(~np.isnan(cl.CVI().flatten())),np.nanstd(cl.CVI().flatten()),
                         np.nanmean(cl.redEdgeNDVI().flatten()),np.nanmedian(cl.redEdgeNDVI().flatten()),np.count_nonzero(~np.isnan(cl.redEdgeNDVI().flatten())),np.nanstd(cl.redEdgeNDVI().flatten()),
                         np.nanmean(cl.GRNDVI().flatten()),np.nanmedian(cl.GRNDVI().flatten()),np.count_nonzero(~np.isnan(cl.GRNDVI().flatten())),np.nanstd(cl.GRNDVI().flatten()),
                         np.nanmean(cl.ATSAVI().flatten()),np.nanmedian(cl.ATSAVI().flatten()),np.count_nonzero(~np.isnan(cl.ATSAVI().flatten())),np.nanstd(cl.ATSAVI().flatten()),
                         np.nanmean(cl.CIgreen().flatten()),np.nanmedian(cl.CIgreen().flatten()),np.count_nonzero(~np.isnan(cl.CIgreen().flatten())),np.nanstd(cl.CIgreen().flatten()),
                         np.nanmean(cl.CIrededge().flatten()),np.nanmedian(cl.CIrededge().flatten()),np.count_nonzero(~np.isnan(cl.CIrededge().flatten())),np.nanstd(cl.CIrededge().flatten()),
                         np.nanmean(cl.CTVI().flatten()),np.nanmedian(cl.CTVI().flatten()),np.count_nonzero(~np.isnan(cl.CTVI().flatten())),np.nanstd(cl.CTVI().flatten()),
                         np.nanmean(cl.GDVI().flatten()),np.nanmedian(cl.GDVI().flatten()),np.count_nonzero(~np.isnan(cl.GDVI().flatten())),np.nanstd(cl.GDVI().flatten()),
                         np.nanmean(cl.GEMI().flatten()),np.nanmedian(cl.GEMI().flatten()),np.count_nonzero(~np.isnan(cl.GEMI().flatten())),np.nanstd(cl.GEMI().flatten()),
                         np.nanmean(cl.GOSAVI().flatten()),np.nanmedian(cl.GOSAVI().flatten()),np.count_nonzero(~np.isnan(cl.GOSAVI().flatten())),np.nanstd(cl.GOSAVI().flatten()),
                         np.nanmean(cl.GSAVI().flatten()),np.nanmedian(cl.GSAVI().flatten()),np.count_nonzero(~np.isnan(cl.GSAVI().flatten())),np.nanstd(cl.GSAVI().flatten()),
                         np.nanmean(cl.IPVI().flatten()),np.nanmedian(cl.IPVI().flatten()),np.count_nonzero(~np.isnan(cl.IPVI().flatten())),np.nanstd(cl.IPVI().flatten()),
                         np.nanmean(cl.RVI().flatten()),np.nanmedian(cl.RVI().flatten()),np.count_nonzero(~np.isnan(cl.RVI().flatten())),np.nanstd(cl.RVI().flatten()),
                         np.nanmean(cl.MRVI().flatten()),np.nanmedian(cl.MRVI().flatten()),np.count_nonzero(~np.isnan(cl.MRVI().flatten())),np.nanstd(cl.MRVI().flatten()),
                         
                         np.nanmean(cl.MSAVI().flatten()),np.nanmedian(cl.MSAVI().flatten()),np.count_nonzero(~np.isnan(cl.MSAVI().flatten())),np.nanstd(cl.MSAVI().flatten()),
                         np.nanmean(cl.NormG().flatten()),np.nanmedian(cl.NormG().flatten()),np.count_nonzero(~np.isnan(cl.NormG().flatten())),np.nanstd(cl.NormG().flatten()),
                         np.nanmean(cl.NormNIR().flatten()),np.nanmedian(cl.NormNIR().flatten()),np.count_nonzero(~np.isnan(cl.NormNIR().flatten())),np.nanstd(cl.NormNIR().flatten()),
                         np.nanmean(cl.NormR().flatten()),np.nanmedian(cl.NormR().flatten()),np.count_nonzero(~np.isnan(cl.NormR().flatten())),np.nanstd(cl.NormR().flatten()),
                         np.nanmean(cl.NGRDI().flatten()),np.nanmedian(cl.NGRDI().flatten()),np.count_nonzero(~np.isnan(cl.NGRDI().flatten())),np.nanstd(cl.NGRDI().flatten()),
                         np.nanmean(cl.RI().flatten()),np.nanmedian(cl.RI().flatten()),np.count_nonzero(~np.isnan(cl.RI().flatten())),np.nanstd(cl.RI().flatten()),
                         np.nanmean(cl.DVI().flatten()),np.nanmedian(cl.DVI().flatten()),np.count_nonzero(~np.isnan(cl.DVI().flatten())),np.nanstd(cl.DVI().flatten()),
                         np.nanmean(cl.TVI().flatten()),np.nanmedian(cl.TVI().flatten()),np.count_nonzero(~np.isnan(cl.TVI().flatten())),np.nanstd(cl.TVI().flatten()),
                         np.nanmean(cl.NDRE().flatten()),np.nanmedian(cl.NDRE().flatten()),np.count_nonzero(~np.isnan(cl.NDRE().flatten())),np.nanstd(cl.NDRE().flatten()),
                         np.nanmean(cl.PSRI().flatten()),np.nanmedian(cl.PSRI().flatten()),np.count_nonzero(~np.isnan(cl.PSRI().flatten())),np.nanstd(cl.PSRI().flatten())
                         ])
```                                       

# Saving the output into a csv table

```python
my_array = np.array(array)
dataframe = pd.DataFrame(my_array, columns = colname)
dataframe.to_csv('python_test3.csv')
```





